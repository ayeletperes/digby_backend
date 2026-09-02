from app import vdjbase_dbs, app, genomic_dbs
from flask import request
from flask_restx import Resource, reqparse
from api.restx import api
from api.system.system import digby_protected
from sqlalchemy import or_, func, cast, Float, distinct
from copy import deepcopy
from functools import lru_cache

from receptor_utils import sequence_alignment

from api.vdjbase.vdjbase import get_vdjbase_species, find_datasets as get_vdjbase_datasets
from api.genomic.genomic import get_genomic_species, get_genomic_datasets
from db.vdjbase_model import Gene as VDJbaseGene, Allele as VDJbaseAllele, AllelesSample as VDJbaseAllelesSample
from db.genomic_db import Gene as GenomicGene, Sequence as GenomicSequence, SampleSequence as GenomicSampleSequence
from db.genomic_airr_model import Sample as GenomicSample, Study as GenomicStudy
from db.vdjbase_airr_model import Sample as VDJbaseSample, Study as VDJbaseStudy

ns = api.namespace('refbook', description='Refbook related operations')


_species_and_loci_cache = None


def species_and_loci():
    """ The species/loci catalogue, merged across the genomic and AIRR-seq databases.

    Every endpoint calls this to validate its arguments, and it walks each open
    dataset, so it previously re-ran the whole scan on every request.
    """
    global _species_and_loci_cache

    if _species_and_loci_cache is not None:
        return deepcopy(_species_and_loci_cache)

    ret = {'species': [], 'loci': {}}

    for get_species, get_datasets in ((get_genomic_species, get_genomic_datasets),
                                      (get_vdjbase_species, get_vdjbase_datasets)):
        for sp in get_species():
            if sp not in ret['species']:
                ret['species'].append(sp)
            for rec in get_datasets(sp):
                loc = rec['dataset']
                # a locus held in both databases must still appear once
                if loc not in ret['loci'].setdefault(sp, []):
                    ret['loci'][sp].append(loc)

    # the open dataset set is fixed when the app starts, so this cannot go stale
    _species_and_loci_cache = ret
    return deepcopy(ret)


def check_species_locus(species, locus):
    """ Returns a 404 response if the species/locus pair is unknown, otherwise None. """
    catalogue = species_and_loci()

    if species not in catalogue['species']:
        return {'message': 'Species not found'}, 404
    if locus not in catalogue['loci'].get(species, []):
        return {'message': 'Locus not found'}, 404

    return None


SOURCES = ('genomic', 'airrseq')


def requested_sources():
    """ The databases this request should read, from a `sources` query parameter.

    Defaults to both, so an omitted parameter behaves as before. Unknown names are
    ignored rather than rejected, and an empty result falls back to both: a request
    that selects nothing is a caller error, not a reason to report no data.
    """
    raw = request.args.get('sources')
    if not raw:
        return set(SOURCES)

    wanted = {s.strip().lower() for s in raw.split(',')} & set(SOURCES)
    return wanted or set(SOURCES)


def requested_list(name):
    """ A comma-separated query parameter as a list, or None when absent.

    None and an empty list mean different things here: no parameter means no
    filtering, whereas an explicit empty value selects nothing.
    """
    raw = request.args.get(name)
    if raw is None:
        return None
    return [value.strip() for value in raw.split(',') if value.strip()]


def applicable(session, Study, Sample, projects, samples):
    """ The part of a project/sample selection that exists in one database.

    The two databases name entirely different studies - AIRR-seq holds P2..P34 and
    genomic holds P25/P28 - so a selection made in one of them matches nothing in
    the other. Filtering both by the whole selection therefore reported zero for
    whichever database the user had not picked from, which reads as "this allele is
    absent" rather than "you did not ask about this database".

    A database with none of the selection is left unfiltered, and `scoped` says so
    for callers that need to describe what a figure covers.
    """
    known_projects = []
    if projects:
        known_projects = [name for (name,) in
                          session.query(Study.study_name)
                          .filter(Study.study_name.in_(projects)).all()]

    known_samples = []
    if samples:
        known_samples = [name for (name,) in
                         session.query(Sample.sample_name)
                         .filter(Sample.sample_name.in_(samples)).all()]

    return known_projects, known_samples, bool(known_projects or known_samples)


def sample_filter(query, projects, samples, session=None):
    """ Narrow a query already joined to Sample by project and by sample name. """
    if session is not None:
        projects, samples, _ = applicable(session, VDJbaseStudy, VDJbaseSample, projects, samples)

    if projects:
        query = query.join(VDJbaseStudy, VDJbaseStudy.id == VDJbaseSample.study_id) \
                     .filter(VDJbaseStudy.study_name.in_(projects))
    if samples:
        query = query.filter(VDJbaseSample.sample_name.in_(samples))
    return query


def genomic_sample_filter(query, projects, samples, session=None):
    """ Narrow a genomic query already joined to Sample by project and sample. """
    if session is not None:
        projects, samples, _ = applicable(session, GenomicStudy, GenomicSample, projects, samples)

    if projects:
        query = query.join(GenomicStudy, GenomicStudy.id == GenomicSample.study_id) \
                     .filter(GenomicStudy.study_name.in_(projects))
    if samples:
        query = query.filter(GenomicSample.sample_name.in_(samples))
    return query


def dataset_session(dbs, species, locus, sources=None, source=None):
    """ Session for one dataset, or None if this database holds nothing for that species/locus.

    `dbs[species]` is absent entirely for a species held only in the other database,
    so indexing it directly raises KeyError rather than missing gracefully.

    Passing `sources` and this database's `source` name also returns None when the
    caller asked not to read it, so selecting genomic-only filters the data and not
    just the set of panels offered.
    """
    if sources is not None and source is not None and source not in sources:
        return None

    dataset = dbs.get(species, {}).get(locus)
    return dataset.session if dataset is not None else None


@ns.route('/species_and_loci')
@api.response(404, 'No species available!')
class SpeciesApi(Resource):
    @digby_protected()
    def get(self):
        """ Returns the list of species and loci for which information is held """

        return species_and_loci()


@ns.route('/projects/<string:species>/<string:locus>')
@api.response(404, 'Species or locus not found')
class ProjectsApi(Resource):
    @digby_protected()
    def get(self, species, locus):
        """ Returns the projects held for a species and locus, per database """

        error = check_species_locus(species, locus)
        if error:
            return error

        sources = requested_sources()
        found = {}

        # the two databases hold different studies entirely - the genomic side has
        # P25/P28 where AIRR-seq has P1..P34 - so which projects exist depends on
        # which databases the caller asked for
        for source, dbs, Study, Sample in (
                ('airrseq', vdjbase_dbs, VDJbaseStudy, VDJbaseSample),
                ('genomic', genomic_dbs, GenomicStudy, GenomicSample)):
            session = dataset_session(dbs, species, locus, sources, source)
            if session is None:
                continue

            rows = session.query(Study.study_name, Study.study_id, func.count(Sample.id)) \
                .join(Sample, Sample.study_id == Study.id) \
                .group_by(Study.study_name, Study.study_id) \
                .all()

            for name, accession, count in rows:
                entry = found.setdefault(name, {'name': name, 'accession': accession,
                                                'samples': 0, 'sources': [],
                                                'by_source': {}})
                entry['samples'] += count
                entry['sources'].append(source)
                # kept split as well as totalled: a project in both databases is
                # one bar per database, and `samples` alone cannot be unstacked
                entry['by_source'][source] = count

        projects = sorted(found.values(), key=lambda p: _project_order(p['name']))
        return {'projects': projects}


def _project_order(name):
    """ P2 before P10: the names are a letter and a number, not plain strings. """
    digits = ''.join(c for c in name if c.isdigit())
    return (name[:1], int(digits) if digits else 0, name)


@ns.route('/samples/<string:species>/<string:locus>')
@api.response(404, 'Species or locus not found')
class SamplesApi(Resource):
    @digby_protected()
    def get(self, species, locus):
        """ Returns the samples held for a species and locus, per database """

        error = check_species_locus(species, locus)
        if error:
            return error

        sources = requested_sources()
        projects = requested_list('projects')
        found = {}

        for source, dbs, Study, Sample in (
                ('airrseq', vdjbase_dbs, VDJbaseStudy, VDJbaseSample),
                ('genomic', genomic_dbs, GenomicStudy, GenomicSample)):
            session = dataset_session(dbs, species, locus, sources, source)
            if session is None:
                continue

            query = session.query(Sample.sample_name, Study.study_name) \
                .join(Study, Study.id == Sample.study_id)
            if projects:
                query = query.filter(Study.study_name.in_(projects))

            for name, project in query.all():
                entry = found.setdefault(name, {'name': name, 'project': project, 'sources': []})
                entry['sources'].append(source)

        samples = sorted(found.values(), key=lambda s: s['name'])
        return {'samples': samples}


@ns.route('/summary')
class RefbookSummary(Resource):
    @digby_protected()
    def get(self):
        """ Headline counts for the landing page

        Samples are counted as distinct names per species+locus, unioned across
        the two databases. That matters in both directions: the human cohorts are
        disjoint, so genomic and AIRR-seq add up, while the rhesus macaque cohort
        was sequenced both ways and its 106 animals would otherwise be counted
        twice.
        """
        catalogue = species_and_loci()
        samples, projects, datasets = 0, set(), 0

        for species in catalogue['species']:
            for locus in catalogue['loci'][species]:
                names = set()
                # joined through Sample, not read off Study directly: a study row
                # with no samples in this locus is not a project anyone can look at
                session = dataset_session(vdjbase_dbs, species, locus, None, 'airrseq')
                if session is not None:
                    rows = (session.query(VDJbaseSample.sample_name, VDJbaseStudy.study_name)
                            .join(VDJbaseStudy, VDJbaseStudy.id == VDJbaseSample.study_id).all())
                    names.update(n for n, _ in rows)
                    projects.update(p for _, p in rows)

                session = dataset_session(genomic_dbs, species, locus, None, 'genomic')
                if session is not None:
                    rows = (session.query(GenomicSample.sample_name, GenomicStudy.study_name)
                            .join(GenomicStudy, GenomicStudy.id == GenomicSample.study_id).all())
                    names.update(n for n, _ in rows)
                    projects.update(p for _, p in rows)

                if names:
                    datasets += 1
                samples += len(names)

        return {'samples': samples, 'projects': len(projects), 'datasets': datasets,
                'species': len(catalogue['species'])}


@ns.route('/ascs_in_locus/<string:species>/<string:locus>')
@api.response(404, 'Species or locus not found')
class AscsInLocusApi(Resource):
    @digby_protected()
    def get(self, species, locus):
        """ Returns the list of ASCs in a given locus for a given species """

        error = check_species_locus(species, locus)
        if error:
            return error

        sources = requested_sources()

        ascs = []
        genomic = False
        airr_seq = False

        # segment comes from the stored type, not the name: IGHA1 and IGHG1 do not
        # carry it, and IGHD is both the delta constant gene and a D-segment prefix
        segments = {}

        # a gene with no sequence has nothing for any panel to draw. The IGHC
        # dataset carries 78 V, 37 D and 6 J gene rows with zero sequences
        # between them, and listing those offers 98 genes that open empty.
        session = dataset_session(vdjbase_dbs, species, locus, sources, 'airrseq')
        if session is not None:
            genes = session.query(VDJbaseGene.name, VDJbaseGene.type) \
                .join(VDJbaseAllele, VDJbaseAllele.gene_id == VDJbaseGene.id) \
                .filter(VDJbaseGene.pseudo_gene == 0).distinct().all()
            ascs.extend([g[0] for g in genes])
            segments.update({g[0]: segment_of_type(g[1]) for g in genes})
            airr_seq = True

        session = dataset_session(genomic_dbs, species, locus, sources, 'genomic')
        if session is not None:
            genes = session.query(GenomicGene.name, GenomicGene.type) \
                .join(GenomicSequence, GenomicSequence.gene_id == GenomicGene.id) \
                .filter(GenomicGene.pseudo_gene == 0).distinct().all()
            ascs.extend([g[0] for g in genes])
            segments.update({g[0]: segment_of_type(g[1]) for g in genes})
            genomic = True

        # both databases contribute: re-deriving this from `genes` dropped whichever
        # ran first, and raised NameError when neither held the locus
        ascs = sorted(set(ascs))
        ascs = [g for g in ascs if '/OR' not in g] # filter orphons
        return {'ascs': ascs, 'segments': {g: segments[g] for g in ascs},
                'genomic': genomic, 'airr_seq': airr_seq}


@ns.route('/ascs_overview/<string:species>/<string:locus>/<path:asc>')
@api.response(404, 'Species or locus not found')
class AscsOverview(Resource):
    @digby_protected()
    def get(self, species, locus, asc):
        """ Returns data for the overview refbook component """

        error = check_species_locus(species, locus)
        if error:
            return error

        sources = requested_sources()
        projects = requested_list('projects')
        samples = requested_list('samples')

        ret = {}
        alleles = {}

        session = dataset_session(vdjbase_dbs, species, locus, sources, 'airrseq')
        if session is not None:
            vdjbase_alleles = session.query(VDJbaseAllele.name, VDJbaseAllele.novel, VDJbaseAllele.appears) \
                .join(VDJbaseGene) \
                .filter(VDJbaseGene.name == asc) \
                .all()
            
            for a, novel, appearances in vdjbase_alleles:
                alleles[a] = {'VDJbase': appearances, 'Genomic': 0, 'novel': novel}

        session = dataset_session(genomic_dbs, species, locus, sources, 'genomic')
        if session is not None:
            genomic_alleles = session.query(GenomicSequence.name, GenomicSequence.novel, GenomicSequence.appearances) \
                .join(GenomicGene) \
                .filter(GenomicGene.name == asc, or_(GenomicSequence.functional == 'Functional', GenomicSequence.functional == 'ORF')) \
                .all()
            
            for a, novel, appearances in genomic_alleles:
                if a not in alleles:
                    alleles[a] = {'VDJbase': 0, 'Genomic': appearances, 'novel': novel}
                else:
                    alleles[a]['Genomic'] = appearances

        wanted = requested_list('alleles')
        if wanted:
            alleles = {name: rec for name, rec in alleles.items() if name in wanted}

        # Allele.appears and Sequence.appearances are whole-dataset totals written
        # when the database was built, so they ignore the project and sample
        # filters and had drifted (248 recorded against 250 actual). Count the
        # samples directly instead, so these agree with the carrier figures.
        #
        # Sample *names*, not counts, because "both" has to be an intersection.
        #
        # Whether that intersection can be non-empty depends on the species. The
        # rhesus macaque data is one cohort sequenced both ways - 106 of 106
        # samples carry the same name in both databases - so an allele really can
        # be seen in the same animal by both. The human databases hold disjoint
        # cohorts (genomic P25/P28, AIRR-seq P2..P34; no shared sample, no shared
        # project), so there the intersection is empty and saying so is the
        # honest answer. `min(genomic, airrseq)`, which this used to report, is
        # neither: it invented 299 shared samples for IGHV1-18*01 out of two
        # cohorts with nobody in common.
        #
        # Whether the per-sample rows are needed at all. Pulling one row per
        # (allele, sample) so the Both column can be an intersection costs ~300ms
        # for a big gene against ~40ms for plain counts, and it is wasted whenever
        # the two cohorts are disjoint - which they are for every human locus,
        # where Both is necessarily zero. The names are fetched only when a sample
        # can actually be in both databases.
        carriers_airrseq = {}
        carriers_genomic = {}
        observed_airrseq = {}
        observed_genomic = {}
        airrseq_scoped = False
        genomic_scoped = False
        airrseq_cohort = set()
        genomic_cohort = set()
        airrseq_session = genomic_session = None

        # Pass one: who is in each cohort. One column per sample, so it is cheap,
        # and it decides whether the per-allele rows below need names at all.
        airrseq_session = dataset_session(vdjbase_dbs, species, locus, sources, 'airrseq')
        if airrseq_session is not None:
            in_scope = sample_filter(airrseq_session.query(VDJbaseSample.sample_name),
                                     projects, samples, airrseq_session).all()
            airrseq_cohort = {name for (name,) in in_scope}
            _, _, airrseq_scoped = applicable(airrseq_session, VDJbaseStudy, VDJbaseSample,
                                              projects, samples)

        genomic_session = dataset_session(genomic_dbs, species, locus, sources, 'genomic')
        if genomic_session is not None:
            in_scope = genomic_sample_filter(genomic_session.query(GenomicSample.sample_name),
                                             projects, samples, genomic_session).all()
            genomic_cohort = {name for (name,) in in_scope}
            _, _, genomic_scoped = applicable(genomic_session, GenomicStudy, GenomicSample,
                                              projects, samples)

        # Pass two. Names are only needed for the intersection, so where no sample
        # is in both databases - every human locus - count in SQL instead.
        shared_cohort = airrseq_cohort & genomic_cohort
        need_names = bool(shared_cohort)

        if airrseq_session is not None:
            base = (
                airrseq_session.query(VDJbaseAllele.name, VDJbaseSample.sample_name)
                .select_from(VDJbaseAllelesSample)
                .join(VDJbaseAllele, VDJbaseAllele.id == VDJbaseAllelesSample.allele_id)
                .join(VDJbaseGene, VDJbaseGene.id == VDJbaseAllele.gene_id)
                .join(VDJbaseSample, VDJbaseSample.id == VDJbaseAllelesSample.sample_id)
                .filter(VDJbaseGene.name == asc)
            )
            if need_names:
                for allele_name, sample_name in sample_filter(base, projects, samples,
                                                              airrseq_session).all():
                    carriers_airrseq.setdefault(allele_name, set()).add(sample_name)
            else:
                counted = base.with_entities(
                    VDJbaseAllele.name, func.count(distinct(VDJbaseAllelesSample.sample_id))
                ).group_by(VDJbaseAllele.name)
                observed_airrseq = dict(sample_filter(counted, projects, samples,
                                                      airrseq_session).all())

        if genomic_session is not None:
            base = (
                genomic_session.query(GenomicSequence.name, GenomicSample.sample_name)
                .select_from(GenomicSampleSequence)
                .join(GenomicSequence, GenomicSequence.id == GenomicSampleSequence.sequence_id)
                .join(GenomicGene, GenomicGene.id == GenomicSequence.gene_id)
                .join(GenomicSample, GenomicSample.id == GenomicSampleSequence.sample_id)
                .filter(GenomicGene.name == asc)
            )
            if need_names:
                for allele_name, sample_name in genomic_sample_filter(base, projects, samples,
                                                                      genomic_session).all():
                    carriers_genomic.setdefault(allele_name, set()).add(sample_name)
            else:
                counted = base.with_entities(
                    GenomicSequence.name, func.count(distinct(GenomicSampleSequence.sample_id))
                ).group_by(GenomicSequence.name)
                observed_genomic = dict(genomic_sample_filter(counted, projects, samples,
                                                              genomic_session).all())

        if need_names:
            observed_airrseq = {name: len(rows) for name, rows in carriers_airrseq.items()}
            observed_genomic = {name: len(rows) for name, rows in carriers_genomic.items()}

        ret['total'] = len(alleles)
        ret['novel'] = sum(rec['novel'] for rec in alleles.values())
        ret['baseline'] = ret['total'] - ret['novel']
        alleles = dict(sorted(alleles.items()))
        ret['alleles'] = list(alleles.keys())
        # these counts may not be exactly what we want, I am not sure what to do if there are samples
        # for which we don't have both genomic and airr-seq results
        # Three exclusive buckets over the same samples, so they add up to the
        # number of distinct samples carrying the allele in either database.
        ret['genomic_only_counts'] = [
            len(carriers_genomic.get(name, set()) - carriers_airrseq.get(name, set()))
            for name in alleles]
        ret['vdjbase_only_counts'] = [
            len(carriers_airrseq.get(name, set()) - carriers_genomic.get(name, set()))
            for name in alleles]
        ret['both_counts'] = [
            len(carriers_genomic.get(name, set()) & carriers_airrseq.get(name, set()))
            for name in alleles]

        # The per-database totals, which the buckets above no longer hide.
        ret['genomic_counts'] = [observed_genomic.get(name, 0) for name in alleles]
        ret['vdjbase_counts'] = [observed_airrseq.get(name, 0) for name in alleles]
        ret['scoped'] = {'genomic': genomic_scoped, 'airrseq': airrseq_scoped}

        # How many samples the figures are drawn from, and how many of them are
        # the same animal in both databases. A UI needs the last number to tell
        # "no allele happens to be shared" apart from "these cohorts have nobody
        # in common, so Both cannot be anything but zero".
        ret['cohort'] = {
            'genomic': len(genomic_cohort),
            'airrseq': len(airrseq_cohort),
            'shared': len(genomic_cohort & airrseq_cohort),
        }

        return ret


def collect_asc_sequences(species, locus, asc, sources=None, allele_names=None,
                          keep_reference=False):
    """ Every allele sequence for one ASC, merged across the requested databases.

    `allele_names` narrows the result to particular alleles, which is what a
    drill-down from another panel asks for.

    Returns a list of {name, seq_gapped, seq}. AIRR-seq is read first so that an
    allele held in both keeps the AIRR-seq record, as before.
    """
    alleles = []
    recs = []

    session = dataset_session(vdjbase_dbs, species, locus, sources, 'airrseq')
    if session is not None:
        vdjbase_alleles = session.query(VDJbaseAllele.name, VDJbaseAllele.seq) \
            .join(VDJbaseGene) \
            .filter(VDJbaseGene.name == asc, ~VDJbaseAllele.name.contains('Del')) \
            .all()

        for a, seq_gapped in vdjbase_alleles:
            if not seq_gapped:      # a row with no sequence would crash on .upper()
                continue
            recs.append({'name': a, 'seq_gapped': seq_gapped.upper(),
                         'seq': seq_gapped.upper().replace('.', '')})
            alleles.append(a)

    session = dataset_session(genomic_dbs, species, locus, sources, 'genomic')
    if session is not None:
        genomic_alleles = session.query(GenomicSequence.name, GenomicSequence.gapped_sequence, GenomicSequence.sequence) \
            .join(GenomicGene) \
            .filter(GenomicGene.name == asc, or_(GenomicSequence.functional == 'Functional', GenomicSequence.functional == 'ORF')) \
            .all()

        for a, gapped, ungapped in genomic_alleles:
            if a in alleles or not gapped:
                continue
            recs.append({'name': a, 'seq_gapped': gapped.upper(),
                         'seq': (ungapped or gapped.replace('.', '')).upper()})

    if allele_names:
        wanted = [rec for rec in recs if rec['name'] in allele_names]

        # An alignment of a single sequence shows nothing, since every row is read
        # against a reference. Keep the gene's first allele as that baseline when
        # the filter would otherwise leave nothing to compare against.
        if keep_reference and len(wanted) < 2 and recs:
            reference = min(recs, key=lambda rec: rec['name'])
            if reference not in wanted:
                wanted = [reference] + wanted

        recs = wanted

    return recs


SEGMENTS = ('V', 'D', 'J', 'C')


def segment_of_type(gene_type):
    """ The segment a gene's stored type says it is: IGHV -> V, IGHC -> C.

    The type is what the database records, so this works for a constant gene,
    whose name does not carry the segment where a V, D or J name does. Prefer it
    to segment_of.
    """
    code = (gene_type or '')[3:4].upper()
    return code if code in SEGMENTS else '?'


def segment_of(asc):
    """ The segment read off an ASC name (IGHV1-2 -> V), where nothing better is at hand.

    Only V, D and J names carry it. IGHA1 and IGHG1 do not, and IGHD is the delta
    constant gene as well as a D-segment prefix, so a name is not enough for a
    constant locus - use segment_of_type against Gene.type there. This returned
    'V' for anything it did not recognise, which made a constant gene a V gene
    silently; it says so now.
    """
    segment = (asc or '')[3:4].upper()
    return segment if segment in SEGMENTS else '?'


def dataset_stamp(species, locus):
    """ A build stamp for the databases behind one species/locus.

    Part of the alignment cache key, so rebuilding a database invalidates its
    entries rather than serving sequences that no longer exist.
    """
    stamp = []
    for dbs in (vdjbase_dbs, genomic_dbs):
        dataset = dbs.get(species, {}).get(locus)
        stamp.append(str(dataset.created) if dataset is not None else None)
    return tuple(stamp)


# Sized to hold every ASC in every dataset at once (~1,140 today, ~7MB of text).
# A smaller cache is actively worse than none for the larger loci: sweeping the 524
# Rhesus IGH ASCs through a 256-entry LRU evicts each entry before it is reused and
# returns a 0% hit rate, paying the bookkeeping for nothing.
NAME_LIMIT = 26


def _fnv1a(text):
    """FNV-1a 32-bit, mirrored byte-for-byte in gene-naming.ts."""
    h = 0x811C9DC5
    for ch in text:
        h ^= ord(ch) & 0xFFFF
        h = (h * 0x01000193) & 0xFFFFFFFF
    return h


def _suffix_token(suffix):
    """Three base36 characters of that hash - stable for an allele, whatever else is loaded."""
    b36 = '0123456789abcdefghijklmnopqrstuvwxyz'
    h = _fnv1a(suffix)
    out = ''
    for _ in range(3):
        out += b36[h % 36]
        h //= 36
    return out


def abbreviate_names(names):
    """ Map allele names to short display labels, and back.

    create_alignment indents every sequence past the longest name, and VDJbase
    encodes each SNP in a novel allele's name: the longest here reaches 193
    characters against a median of 23, which would push the sequences ~200
    columns to the right and make the alignment unreadable.

    A name over the limit is shortened to its stem plus the number of suffixes it
    carried, so IGHV1-2*02_t211c_t213c_g225a becomes IGHV1-2*02+3.
    """
    labels = {}
    used = set()

    # Sorted, so a `#2` lands on the same allele as it does in the client's
    # shortenAlleleNames. Assignment depends on iteration order, and the two
    # sides receive their names in different orders.
    for name in sorted(names):
        if len(name) <= NAME_LIMIT:
            label = name
        else:
            # split at the first _ AFTER the allele: gene names themselves carry
            # underscores (IGHV4-NL_1*01_a157g), and partitioning on the first
            # one dropped the allele and miscounted the substitutions - 458 of
            # allele_server's 3,449 suffixed IG names are affected
            star = name.find('*')
            cut = name.find('_', star + 1) if star >= 0 else name.find('_')
            if cut >= 0:
                suffix = name[cut + 1:]
                # the count alone is not distinguishing: 231 of the 793 names
                # over the limit in the full HUSA set share a stem and a count.
                # The token comes from the allele's own suffix, so it survives
                # more data being loaded, which a positional #N does not.
                label = f'{name[:cut]}+{len(suffix.split("_"))}~{_suffix_token(suffix)}'
            else:
                label = name[:NAME_LIMIT - 1] + '~'

        # labels index the alignment rows, so they have to stay distinct
        candidate, n = label, 2
        while candidate in used:
            candidate = f'{label}#{n}'
            n += 1

        used.add(candidate)
        labels[candidate] = name

    return labels


@lru_cache(maxsize=4096)
def _render_alignment(species, locus, asc, stamp, codon_wrap, sources, allele_names):
    """ Rendered alignment for one ASC. Memoised: see dataset_stamp for invalidation.

    `stamp` is unused in the body and present only to key the cache. `sources` is a
    frozenset so that it stays hashable and genomic-only renders cache separately.
    """
    by_name = {r['name']: r['seq_gapped']
               for r in collect_asc_sequences(species, locus, asc, sources,
                                              list(allele_names) if allele_names else None,
                                              keep_reference=True)}
    if not by_name:
        return None

    labels = abbreviate_names(by_name)
    sequences = {label: by_name[name] for label, name in labels.items()}

    alignment = sequence_alignment.create_alignment(
        sequences, sequence_type=segment_of(asc), codon_wrap=codon_wrap)

    # only the abbreviated ones are worth listing
    legend = {label: name for label, name in labels.items() if label != name}
    return alignment, legend


@ns.route('/asc_alignment/<string:species>/<string:locus>/<path:asc>')
@api.response(404, 'Species or locus not found')
class AscAlignment(Resource):
    @digby_protected()
    def get(self, species, locus, asc):
        """ Returns a formatted alignment of every allele in an ASC """

        error = check_species_locus(species, locus)
        if error:
            return error

        sources = requested_sources()

        try:
            codon_wrap = max(5, min(60, int(request.args.get('wrap', 20))))
        except ValueError:
            return {'message': 'wrap must be a number'}, 400

        alleles = requested_list('alleles')
        rendered = _render_alignment(species, locus, asc, dataset_stamp(species, locus),
                                     codon_wrap, frozenset(sources),
                                     frozenset(alleles) if alleles else None)

        if rendered is None:
            return {'message': f'No sequences for {asc}'}, 404

        alignment, legend = rendered
        info = _render_alignment.cache_info()
        return {'asc': asc, 'segment': segment_of(asc), 'codon_wrap': codon_wrap,
                'alignment': alignment, 'legend': legend,
                'cache': {'hits': info.hits, 'misses': info.misses, 'size': info.currsize}}


@ns.route('/asc_seqs/<string:species>/<string:locus>/<path:asc>')
@api.response(404, 'Species or locus not found')
class AscSeqs(Resource):
    @digby_protected()
    def get(self, species, locus, asc):
        """ Returns sequences of all alleles in an ASC """

        error = check_species_locus(species, locus)
        if error:
            return error

        sources = requested_sources()

        return {'alleles': collect_asc_sequences(species, locus, asc, sources, requested_list('alleles'))}

@ns.route('/asc_usage/<string:species>/<string:locus>/<path:asc>')
@api.response(404, 'Species or locus not found')
class AscUsage(Resource):
    @digby_protected()
    def get(self, species, locus, asc):
        """ Returns usage statistics for all alleles in an ASC """

        error = check_species_locus(species, locus)
        if error:
            return error

        sources = requested_sources()
        projects = requested_list('projects')
        samples = requested_list('samples')
        alleles = requested_list('alleles')

        vdjbase_alleles = []

        session = dataset_session(vdjbase_dbs, species, locus, sources, 'airrseq')
        if session is not None:

            totals = (
                session.query(
                    VDJbaseAllelesSample.patient_id,
                    VDJbaseAllelesSample.sample_id,
                    func.sum(VDJbaseAllelesSample.total_count).label("total")
                )
                .group_by(VDJbaseAllelesSample.patient_id, VDJbaseAllelesSample.sample_id)
                .subquery()
            )
            
            fraction_expr = (
                cast(VDJbaseAllelesSample.count, Float) /
                func.nullif(cast(totals.c.total, Float), 0.0)
            )

            usage_query = (
                session.query(
                    VDJbaseAllele.name,
                    func.group_concat(VDJbaseSample.sample_name).label("samples"),
                    func.group_concat(func.coalesce(fraction_expr, 0.0)).label("fraction"),
                )
                .join(VDJbaseAllele, VDJbaseAllele.id == VDJbaseAllelesSample.allele_id)
                .join(VDJbaseSample, VDJbaseSample.id == VDJbaseAllelesSample.sample_id)
                .join(
                    totals,
                    (VDJbaseAllelesSample.patient_id == totals.c.patient_id)
                    & (VDJbaseAllelesSample.sample_id == totals.c.sample_id),
                )
                .join(VDJbaseGene, VDJbaseGene.id == VDJbaseAllele.gene_id)
                .filter(
                    VDJbaseGene.name == asc,
                    VDJbaseAllelesSample.count.isnot(None)
                )
            )

            if alleles:
                usage_query = usage_query.filter(VDJbaseAllele.name.in_(alleles))

            vdjbase_alleles = (
                sample_filter(usage_query, projects, samples, session)
                .group_by(VDJbaseAllele.name)
                .all()
            )

        recs = [{'name': allele, 'usage': list(usages.split(',') if usages else []), 'samples': list(samples.split(',') if samples else [])} for allele, samples, usages in vdjbase_alleles]

        return {'alleles': recs}

@ns.route('/asc_zygousity/<string:species>/<string:locus>/<path:asc>')
@api.response(404, 'Species or locus not found')
class AscZygosity(Resource):
    @digby_protected()
    def get(self, species, locus, asc):
        """ Returns zygosity statistics for all subjects in a given ASC """

        error = check_species_locus(species, locus)
        if error:
            return error

        sources = requested_sources()
        projects = requested_list('projects')
        samples = requested_list('samples')
        alleles = requested_list('alleles')

        # Keyed on sample name and unioned across the databases, because the unit
        # of zygosity is the subject, not the record. Where the same subject was
        # sequenced both ways - the rhesus macaque cohort is, 106 of 106 names
        # match - the two databases describe one animal and must not become two
        # rows. Where the cohorts are disjoint, as in the human data, the union is
        # simply a concatenation.
        carried = {}

        session = dataset_session(vdjbase_dbs, species, locus, sources, 'airrseq')
        if session is not None:

            zygosity_query = (
                session.query(
                    VDJbaseSample.sample_name,
                    func.group_concat(func.distinct(VDJbaseAllele.name)).label("alleles"),
                )
                # without this the first selected column decides the FROM, and the
                # explicit join below then adds `sample` a second time
                .select_from(VDJbaseAllelesSample)
                .join(VDJbaseAllele, VDJbaseAllele.id == VDJbaseAllelesSample.allele_id)
                .join(VDJbaseGene, VDJbaseGene.id == VDJbaseAllele.gene_id)
                .join(VDJbaseSample, VDJbaseSample.id == VDJbaseAllelesSample.sample_id)
                .filter(
                    VDJbaseGene.name == asc,
                )
            )

            if alleles:
                carriers = (
                    session.query(VDJbaseAllelesSample.sample_id)
                    .join(VDJbaseAllele, VDJbaseAllele.id == VDJbaseAllelesSample.allele_id)
                    .filter(VDJbaseAllele.name.in_(alleles))
                    .subquery()
                )
                zygosity_query = zygosity_query.filter(VDJbaseSample.id.in_(carriers))

            alleles_per_sample = (
                sample_filter(zygosity_query, projects, samples, session)
                .group_by(VDJbaseSample.id)
                .all()
            )

            # not `alleles`: that name holds the requested allele filter, and
            # rebinding it here left the genomic branch below filtering on a
            # row's comma-joined string instead of the caller's list
            for sample_name, carried_names in alleles_per_sample:
                if carried_names:
                    carried.setdefault(sample_name, set()).update(carried_names.split(','))

        # Genomic carries the same information: which alleles of the gene a
        # sample holds. Only *usage* is genuinely AIRR-seq-only, so the panel
        # used to be gated on a database it did not actually need.
        session = dataset_session(genomic_dbs, species, locus, sources, 'genomic')
        if session is not None:

            genomic_query = (
                session.query(
                    GenomicSample.sample_name,
                    func.group_concat(func.distinct(GenomicSequence.name)).label("alleles"),
                )
                # as above: let the association table decide the FROM, so the
                # explicit join does not add `sample` a second time
                .select_from(GenomicSampleSequence)
                .join(GenomicSequence, GenomicSequence.id == GenomicSampleSequence.sequence_id)
                .join(GenomicGene, GenomicGene.id == GenomicSequence.gene_id)
                .join(GenomicSample, GenomicSample.id == GenomicSampleSequence.sample_id)
                .filter(
                    GenomicGene.name == asc,
                    or_(GenomicSequence.functional == 'Functional',
                        GenomicSequence.functional == 'ORF'),
                )
            )

            if alleles:
                carriers = (
                    session.query(GenomicSampleSequence.sample_id)
                    .join(GenomicSequence,
                          GenomicSequence.id == GenomicSampleSequence.sequence_id)
                    .filter(GenomicSequence.name.in_(alleles))
                    .subquery()
                )
                genomic_query = genomic_query.filter(GenomicSample.id.in_(carriers))

            genomic_per_sample = (
                genomic_sample_filter(genomic_query, projects, samples, session)
                .group_by(GenomicSample.id)
                .all()
            )

            for sample_name, carried_names in genomic_per_sample:
                if carried_names:
                    carried.setdefault(sample_name, set()).update(carried_names.split(','))

        # The subqueries above pick the subjects who carry one of the requested
        # alleles; they do not say which alleles to draw, so every other allele
        # those subjects hold came back too and the "seen in" thresholds had no
        # visible effect here. Intersect, and drop a subject left holding none.
        if alleles:
            wanted = set(alleles)
            carried = {name: sets & wanted for name, sets in carried.items()}
            carried = {name: sets for name, sets in carried.items() if sets}

        recs = [{'name': name, 'sets': sorted(sets)}
                for name, sets in sorted(carried.items())]

        return {'samples': recs}
