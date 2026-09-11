"""Services for gene-usage QTL results.

Synchronous JSON, the same shape as the refbook namespace: these are single
indexed reads over a local sqlite file, so there is nothing for Celery to do.

The Manhattan view is deliberately an aggregate. The association table holds one
row per variant *and* ASC - 658,140 of them for IGH - but a Manhattan plot has
one point per variant, so when no ASC is named the strongest signal per variant
is returned. That is 9,402 points rather than 658,140, and it is the plot people
actually mean, so no thinning or sampling is needed.
"""

import math

from flask import request
from flask_restx import Resource
from sqlalchemy import Integer, cast, func, select

from api.restx import api
from api.system.system import digby_protected
from api.qtl import dbsnp
from app import qtl_dbs
from db.qtl_model import (
    Asc, AscUsage, Dosage, Subject, Threshold, UsageAssociation, Variant,
)

ns = api.namespace('qtl', description='Gene-usage QTL results')

# Loci offered by the API. IGH was held back while its cohort had no linked
# repertoire in VDJbase; that gate only ever mattered to analyses that join the
# two, and none of the guQTL views do, so every built locus is offered. Kept as a
# set rather than deleted: it is the one place to hold a locus back again.
HIDDEN_LOCI: set = set()

# The smallest genotype class a fit needs before the run calls it well powered.
# Read off the run's own published leads rather than assumed: across all 1,228
# rows carrying the flag in IGK and IGL, every one at 4 or below is flagged not
# well powered and every one at 5 or above is flagged well powered, with no
# exceptions in either locus. Applying it to the rest of the table says the same
# thing about the same numbers, rather than leaving the caveat off most rows.
WELL_POWERED_MIN = 5


# What each built database says it holds, keyed by species. Read once: the set of
# databases is fixed at startup.
_catalogue = None


def _describe(key, provider):
    """The locus and project a database records about itself.

    Asked of the database, never of the directory it sits in. A directory has to
    be uniquely named so two projects can both hold an IGH, which makes its name
    a filing label; taking a fact back out of it is how `IGKC` ended up filed as
    a V gene next door.
    """
    session = provider.session
    try:
        locus = session.execute('select locus from details limit 1').scalar()
    except Exception:
        locus = None
    try:
        projects = [r[0] for r in session.execute('select project from qtl_run order by id')]
    except Exception:
        # built before the table existed; unknown, which is not the same as none
        projects = []
    return (locus or key), (projects or [None])


def datasets(species=None):
    """Every offered database, as {species: [{locus, project, provider}, ...]}."""
    global _catalogue
    if _catalogue is None:
        _catalogue = {}
        for name in sorted(qtl_dbs):
            found = []
            for key in sorted(qtl_dbs[name]):
                provider = qtl_dbs[name][key]
                locus, projects = _describe(key, provider)
                if locus in HIDDEN_LOCI:
                    continue
                # one entry per study the database holds, all sharing its provider
                for project in projects:
                    found.append({'locus': locus, 'project': project,
                                  'provider': provider, 'studies': len(projects)})
            if found:
                _catalogue[name] = found
    return _catalogue.get(species, []) if species else _catalogue


def loci_for(species, project=None):
    """The loci one project holds, for the views that sweep every locus at once.

    Without this those sweeps ran over every locus any project holds, and each
    one that resolved to two databases was skipped in silence - a search that
    quietly stopped covering a locus as soon as a second study was loaded.
    """
    if project is None:
        project = current_project()
    return sorted({d['locus'] for d in datasets(species)
                   if project is None or d['project'] == project})


def qtl_session(species, locus, project=None):
    """Session for one guQTL dataset, or None if it is absent or ambiguous.

    A scan is computed within one cohort and never pooled across them, so a
    study qualifies every query. The schema carries that as qtl_run, and the
    dimension tables carry run_id, but the queries here do not filter on it yet.

    So a database holding more than one study is refused rather than answered
    from: every view would silently mix two cohorts whose p-values are on
    different scales. Serving the wrong cohort's numbers under the right locus
    name is the one outcome worth a hard failure. Until the run filter reaches
    the queries, one study per database is the supported case.
    """
    if project is None:
        project = current_project()

    found = [d for d in datasets(species) if d['locus'] == locus]
    if any(d.get('studies', 1) > 1 for d in found):
        return None

    if project is not None:
        found = [d for d in found if d['project'] == project]

    return found[0]['provider'].session if len(found) == 1 else None


def current_project():
    """The project this request asked for, if it named one.

    Read here rather than declared on twelve routes: it qualifies every one of
    them the same way, and a route that forgot to accept it would quietly answer
    from whichever database sorted first.
    """
    try:
        value = request.args.get('project')
    except RuntimeError:      # outside a request context
        return None
    return value or None


def available():
    """Species, projects and loci that have guQTL results and are offered."""
    ret = {'species': [], 'loci': {}, 'projects': {}, 'datasets': []}

    for species, found in datasets().items():
        ret['species'].append(species)
        ret['loci'][species] = sorted({d['locus'] for d in found})
        # '' is a database built before projects were recorded, which is a
        # different statement from a database whose project is known
        ret['projects'][species] = sorted({d['project'] or '' for d in found})
        ret['datasets'] += [{'species': species, 'locus': d['locus'],
                             'project': d['project']}
                            for d in sorted(found, key=lambda d: (d['project'] or '',
                                                                  d['locus']))]

    return ret


def display_name(locus, asc):
    """An ASC under the name someone would type, without doubling the locus.

    IGH's D clusters are stored with the locus already on them - `IGHD5-12` -
    while V and J are stored bare - `V5-10-1`. Prefixing either way round without
    looking produced `IGHIGHD5-12`, which matches nothing.
    """
    asc = str(asc)
    return asc if asc.startswith(locus) else f'{locus}{asc}'


def gene_names(locus, asc):
    """Every gene name an ASC answers to, for searching.

    ASC names drop the locus, so IGL's `V9-49` is the gene IGLV9-49. A name
    carrying a slash is a merged cluster of genes that cannot be told apart -
    `V1-13/1D-13` is IGKV1-13 together with IGKV1D-13 - and the members after the
    first are written without their segment letter, so it is put back. Someone
    who types either member should find the cluster that contains it.
    """
    parts = str(asc).split('/')
    # IGH's D clusters are stored with the locus already on them - `IGHD5-12`,
    # and `IGHD4-11/IGHD4-4` for a merged one - while V and J are not. Prefixing
    # blindly made `IGHIGHD5-12`, which matches nothing anyone would type. The
    # stored string is left alone either way; this only builds search aliases.
    def qualify(name):
        return display_name(locus, name)

    segment = parts[0].replace(locus, '', 1)[:1]

    names = {asc, qualify(parts[0])}
    for member in parts[1:]:
        if member.startswith(locus):
            names.add(member)
        else:
            names.add(qualify(member if member.startswith(segment)
                              else f'{segment}{member}'))
    return names


def check_species_locus(species, locus, project=None):
    """A 404 or 400 response if this request names no single dataset, else None."""
    catalogue = available()

    if species not in catalogue['species']:
        return {'message': 'Species not found'}, 404
    if locus not in catalogue['loci'].get(species, []):
        return {'message': 'Locus not found'}, 404

    if project is None:
        project = current_project()
    found = [d for d in datasets(species) if d['locus'] == locus]
    if project is not None:
        named = [d for d in found if d['project'] == project]
        if not named:
            return {'message': f'No {locus} results for project {project}'}, 404
        found = named

    if len(found) > 1:
        return {'message': f'{len(found)} projects hold {locus} results; name one '
                           'with ?project=. These analyses are per project and '
                           'are not pooled across them.'}, 400

    return None


@ns.route('/species_and_loci')
class QtlSpeciesApi(Resource):
    @digby_protected()
    def get(self):
        """ Returns the species and loci for which guQTL results are held """
        return available()


@ns.route('/ascs/<string:species>/<string:locus>')
@api.response(404, 'Species or locus not found')
class QtlAscsApi(Resource):
    @digby_protected()
    def get(self, species, locus):
        """ Returns the ASCs tested in a locus, with how strong their best hit was """

        error = check_species_locus(species, locus)
        if error:
            return error

        session = qtl_session(species, locus)

        # Three cheap queries instead of one expensive one. Asking for the count,
        # the significant total and the best p in a single GROUP BY costs 558 ms
        # on IGH's 658,140 associations, because summing `significant` has to
        # read every row; split out, the other two ride `ix_usage_asc_p` as
        # covering scans and the significant count reads only the 11,122 rows
        # that are significant. 558 ms -> 191 ms, same numbers.
        tested = dict(session.query(UsageAssociation.asc_id, func.count())
                      .group_by(UsageAssociation.asc_id).all())

        # cast first: summing a Boolean column runs the total back through the
        # Boolean result processor, so 206 arrives as True and counts as 1.
        # Counting a filtered set avoids the sum, and the cast with it.
        significant = dict(session.query(UsageAssociation.asc_id, func.count())
                           .filter(UsageAssociation.significant == True)
                           .group_by(UsageAssociation.asc_id).all())

        # `GROUP BY asc_id` walks all 658,140 index entries to find 70 maxima.
        # Asked one ASC at a time it is 70 seeks instead: with `ix_usage_asc_p`
        # on (asc_id, neglog10_p) each max is the last entry of its own range,
        # which SQLite reaches directly. 94ms -> 9ms, and no new index for it.
        best_of = (select([func.max(UsageAssociation.neglog10_p)])
                   .where(UsageAssociation.asc_id == Asc.id)
                   .correlate(Asc)
                   .as_scalar())
        best = dict(session.query(Asc.id, best_of).all())

        ascs = [{'asc': row.asc, 'segment': row.segment, 'n_member': row.n_member,
                 'n_variants': int(tested.get(row.id, 0)),
                 'n_significant': int(significant.get(row.id, 0)),
                 'best_neglog10_p': best.get(row.id)}
                for row in session.query(Asc).all()
                # an ASC with no associations was not scanned, and was not in the
                # joined form of this query either
                if row.id in tested]
        ascs.sort(key=lambda a: (a['segment'] or '', a['asc']))

        return {'ascs': ascs, 'thresholds': _thresholds(session)}


_summary_columns = {}


def _has_summary(session):
    """Whether this database carries the precomputed best-per-variant columns.

    A schema question, asked once per database, not a data question: the answer
    chooses between two queries that return the same numbers, so an older
    database is slower here and never wrong.
    """
    key = id(session.bind)
    if key not in _summary_columns:
        try:
            session.query(Variant.best_neglog10_p).limit(1).all()
            _summary_columns[key] = True
        except Exception:
            session.rollback()
            _summary_columns[key] = False
    return _summary_columns[key]


def _thresholds(session):
    return [{'analysis': t.analysis, 'conditional': t.conditional or None,
             'grouped_by': t.grouped_by or None, 'threshold': t.threshold,
             'neglog10_threshold': -math.log10(t.threshold) if t.threshold else None,
             'n_subjects': t.n_subjects, 'n_variants': t.n_variants,
             'n_independent': t.n_independent, 'n_asc': t.n_asc,
             'n_significant_variants': t.n_significant_variants}
            for t in session.query(Threshold).all()]


@ns.route('/manhattan/<string:species>/<string:locus>')
@api.response(404, 'Species or locus not found')
class QtlManhattanApi(Resource):
    @digby_protected()
    def get(self, species, locus):
        """ Returns one point per variant: its position and how strong the signal is

        With `asc`, the association with that ASC. Without, the strongest signal
        the variant showed against any ASC, which is the usual whole-locus view.
        """

        error = check_species_locus(species, locus)
        if error:
            return error

        session = qtl_session(species, locus)
        asc = request.args.get('asc')

        # Each ASC is its own scan, so a Manhattan is properly per-ASC. Without
        # one this returns each variant's strongest signal across all of them,
        # which is a summary rather than a test - so it also reports which ASC
        # produced that maximum, or a click on a point would not know whose usage
        # to plot.
        if asc:
            # One row per (variant, ASC) - the table is unique on the pair - so
            # there is nothing to group and no maximum to take. Grouping anyway
            # cost a temp b-tree over the filtered set for a result that was
            # already one row per variant.
            query = (
                session.query(Variant.variant, Variant.pos, Variant.maf,
                              Variant.gene, Variant.feature,
                              UsageAssociation.neglog10_p,
                              cast(UsageAssociation.significant, Integer),
                              Asc.asc)
                .join(UsageAssociation, UsageAssociation.variant_id == Variant.id)
                .join(Asc, Asc.id == UsageAssociation.asc_id)
                .filter(Asc.asc == asc)
            )
        elif _has_summary(session):
            # Read the answer rather than recompute it. The maxima are stored on
            # the variant at build time; see Variant.best_neglog10_p.
            query = (
                session.query(Variant.variant, Variant.pos, Variant.maf,
                              Variant.gene, Variant.feature,
                              Variant.best_neglog10_p,
                              cast(Variant.best_significant, Integer),
                              Asc.asc)
                .outerjoin(Asc, Asc.id == Variant.best_asc_id)
                .filter(Variant.best_neglog10_p.isnot(None))
            )
        else:
            # A database built before the summary columns existed. Same numbers,
            # several seconds slower - the aggregate this view was named for.
            query = (
                session.query(Variant.variant, Variant.pos, Variant.maf,
                              Variant.gene, Variant.feature,
                              func.max(UsageAssociation.neglog10_p),
                              cast(UsageAssociation.significant, Integer),
                              Asc.asc)
                .join(UsageAssociation, UsageAssociation.variant_id == Variant.id)
                .join(Asc, Asc.id == UsageAssociation.asc_id)
                .group_by(Variant.id)
            )

        points = [{'variant': variant, 'pos': pos, 'maf': maf, 'gene': gene,
                   'feature': feature, 'neglog10_p': neglog10_p,
                   'significant': bool(significant), 'asc': best_asc}
                  for variant, pos, maf, gene, feature, neglog10_p, significant,
                      best_asc in query.all()]
        points.sort(key=lambda p: p['pos'] if p['pos'] is not None else 0)

        # the contig the positions are on. IGH sits on a locus-relative contig
        # literally named `igh` while the light chains use chr2 / chr22, so the
        # axis has to say which frame it is drawn in rather than name the locus
        contig = session.query(Variant.contig).filter(Variant.contig.isnot(None)).first()

        return {'locus': locus, 'asc': asc, 'points': points,
                'contig': contig[0] if contig else None,
                'thresholds': _thresholds(session),
                'leads': _leads(session, asc)}


def _leads(session, asc=None, limit=10):
    """The strongest independent signals, for labelling the plot."""
    query = (
        session.query(Variant.variant, Variant.pos, Asc.asc, UsageAssociation.neglog10_p,
                      UsageAssociation.beta, Variant.min_genotype_group,
                      Variant.well_powered)
        .join(UsageAssociation, UsageAssociation.variant_id == Variant.id)
        .join(Asc, Asc.id == UsageAssociation.asc_id)
        .filter(UsageAssociation.is_lead == True)      # noqa: E712 - SQL, not Python
    )
    if asc:
        query = query.filter(Asc.asc == asc)

    rows = query.order_by(UsageAssociation.neglog10_p.desc()).limit(limit).all()
    return [{'variant': variant, 'pos': pos, 'asc': asc_name, 'neglog10_p': p,
             'beta': beta, 'min_genotype_group': min_group,
             'well_powered': bool(powered) if powered is not None else None}
            for variant, pos, asc_name, p, beta, min_group, powered in rows]


def _genotype_counts(session, record):
    """How many subjects carry each genotype at a variant.

    A property of the variant, not of the scan: the cohort is the same for every
    ASC, and the only subjects a scan drops are the ones with no call here, which
    is what `genotype is not null` already excludes.
    """
    rows = (
        session.query(Dosage.genotype, func.count(Dosage.id))
        .filter(Dosage.variant_id == record.id, Dosage.genotype.isnot(None))
        .group_by(Dosage.genotype)
        .all()
    )
    return {int(genotype): count for genotype, count in rows}


def _variant_payload(session, record):
    """A variant's identity, and every ASC scan it appeared in.

    Shared by the two ways in: from a Manhattan point, where the locus is already
    known, and from a bare variant id, where it is not.

    `min_genotype_group` is filled in here. The run publishes it only for lead
    variants - 988 of IGK's 7,255 significant rows, 240 of IGL's 1,573 - which
    would leave the column that qualifies every p-value blank on most of the
    table. It is recovered from the genotype counts instead, which reproduces the
    run's own value exactly wherever the run published one. `well_powered`
    follows from it by the run's own rule - see WELL_POWERED_MIN.
    """
    counts = _genotype_counts(session, record)
    smallest = min(counts.values()) if counts else None

    rows = (
        # n, min_genotype_group and well_powered are the variant's, not each
        # association's, so they come off `record` rather than a join
        session.query(Asc.asc, Asc.segment, UsageAssociation.beta, UsageAssociation.se,
                      UsageAssociation.p_value, UsageAssociation.neglog10_p,
                      UsageAssociation.significant, UsageAssociation.is_lead)
        .join(UsageAssociation, UsageAssociation.asc_id == Asc.id)
        .filter(UsageAssociation.variant_id == record.id)
        .order_by(UsageAssociation.neglog10_p.desc())
        .all()
    )

    associations = [
        {'asc': asc, 'segment': segment, 'beta': beta, 'se': se, 'p_value': p,
         'neglog10_p': neglog10_p, 'significant': bool(significant), 'n': record.n,
         'min_genotype_group': (record.min_genotype_group
                                if record.min_genotype_group is not None else smallest),
         'well_powered': (bool(record.well_powered) if record.well_powered is not None
                          else None if smallest is None
                          else smallest >= WELL_POWERED_MIN),
         'is_lead': bool(lead)}
        for asc, segment, beta, se, p, neglog10_p, significant, lead in rows]

    return {
        'variant': {'variant': record.variant, 'contig': record.contig, 'pos': record.pos,
                    'maf': record.maf, 'gene': record.gene, 'feature': record.feature,
                    'sub_feature': record.sub_feature,
                    'distance_to_gene': record.distance_to_gene},
        'associations': associations,
        # a variant is tested against every ASC in the locus, so "how many genes
        # does it actually drive" is the number the table is really asked for
        'n_tested': len(associations),
        'n_significant': sum(1 for a in associations if a['significant']),
        # the same three numbers the boxplot is drawn over, stated once: they
        # qualify every row of the table rather than any one of them
        'genotype_counts': counts,
        'min_genotype_group': smallest,
        # what the world outside VDJbase calls this variant, when it can be said.
        # IGH is named for a locus-relative contig, so it takes a map.
        'dbsnp': dbsnp.lookup(record.contig, record.pos),
        'thresholds': _thresholds(session),
        'has_genotypes': session.query(Dosage.id)
            .filter(Dosage.variant_id == record.id).first() is not None,
    }


@ns.route('/variant/<string:species>/<string:locus>/<path:variant>')
@api.response(404, 'Species, locus or variant not found')
class QtlVariantApi(Resource):
    @digby_protected()
    def get(self, species, locus, variant):
        """ Returns what a variant is, and every ASC it was tested against """

        error = check_species_locus(species, locus)
        if error:
            return error

        session = qtl_session(species, locus)
        record = session.query(Variant).filter(Variant.variant == variant).one_or_none()
        if record is None:
            return {'message': f'No such variant: {variant}'}, 404

        payload = _variant_payload(session, record)
        payload['locus'] = locus
        return payload


@ns.route('/variant_lookup/<string:species>/<path:variant>')
@api.response(404, 'Species or variant not found')
class QtlVariantLookupApi(Resource):
    @digby_protected()
    def get(self, species, variant):
        """ Resolves a bare variant id to its locus and every ASC it was tested against

        The lookup people actually arrive with. Someone holding an id from a GWAS
        hit knows the id and nothing else - not which locus it falls in, and not
        which genes were scanned against it - so every offered locus is searched
        for an exact match rather than the caller being asked to guess.

        Ids carry their contig (`chr22_22756855`, `igh_...`) and the loci sit on
        different contigs, so a match is unambiguous and the first one is the one.
        """
        catalogue = available()
        if species not in catalogue['species']:
            return {'message': 'Species not found'}, 404

        for locus in loci_for(species):
            session = qtl_session(species, locus)
            if session is None:
                continue

            record = (session.query(Variant)
                      .filter(Variant.variant == variant).one_or_none())
            if record is None:
                continue

            payload = _variant_payload(session, record)
            payload['locus'] = locus
            return payload

        return {'message': f'No such variant: {variant}'}, 404


@ns.route('/variant_usage/<string:species>/<string:locus>/<path:variant>')
@api.response(404, 'Species, locus or variant not found')
class QtlVariantUsageApi(Resource):
    @digby_protected()
    def get(self, species, locus, variant):
        """ Returns each subject's genotype at a variant against their ASC usage

        This is the plot behind a Manhattan hit: the association says a variant
        explains usage, and this shows the usage it explains, per subject.
        """

        error = check_species_locus(species, locus)
        if error:
            return error

        session = qtl_session(species, locus)
        asc = request.args.get('asc')
        if not asc:
            return {'message': 'An asc is required'}, 400

        record = session.query(Variant).filter(Variant.variant == variant).one_or_none()
        if record is None:
            return {'message': f'No such variant: {variant}'}, 404

        asc_record = session.query(Asc).filter(Asc.asc == asc).one_or_none()
        if asc_record is None:
            return {'message': f'No such ASC: {asc}'}, 404

        # one row per subject: their genotype here, and their usage of that ASC
        rows = (
            session.query(Subject.subject, Dosage.dosage, Dosage.genotype,
                          AscUsage.usage, AscUsage.logit_usage, AscUsage.count,
                          AscUsage.total)
            .select_from(Dosage)
            .join(Subject, Subject.id == Dosage.subject_id)
            .join(AscUsage, (AscUsage.subject_id == Dosage.subject_id) &
                            (AscUsage.asc_id == asc_record.id))
            .filter(Dosage.variant_id == record.id)
            .all()
        )

        association = (
            session.query(UsageAssociation)
            .filter(UsageAssociation.variant_id == record.id,
                    UsageAssociation.asc_id == asc_record.id)
            .one_or_none()
        )

        subjects = [{'subject': subject, 'dosage': dosage, 'genotype': genotype,
                     'usage': usage, 'logit_usage': logit, 'count': count, 'total': total}
                    for subject, dosage, genotype, usage, logit, count, total in rows]

        return {
            'variant': variant,
            'asc': asc,
            'segment': asc_record.segment,
            'subjects': subjects,
            'association': None if association is None else {
                'beta': association.beta, 'se': association.se,
                'p_value': association.p_value, 'neglog10_p': association.neglog10_p,
                'n': association.n, 'significant': bool(association.significant),
                # the pipeline documents the extreme tail as anti-conservative, so
                # the smallest genotype class travels with the p-value
                'min_genotype_group': association.min_genotype_group,
                'well_powered': (bool(association.well_powered)
                                 if association.well_powered is not None else None),
            },
        }


@ns.route('/search/<string:species>')
@api.response(404, 'Species not found')
class QtlSearchApi(Resource):
    @digby_protected()
    def get(self, species):
        """ Resolves a typed query to variants and genes across every locus

        The two questions people arrive with are "I have a variant, is it a QTL
        here" and "I have a gene, does anything drive its usage" - neither of
        which starts from a Manhattan plot. Both are answered from the same box,
        so the caller does not have to know which locus to look in first.
        """
        query = (request.args.get('q') or '').strip()
        # each side of the dashboard asks only about its own kind, so the other
        # scan is skipped rather than run and discarded
        kind = request.args.get('kind') or 'both'
        want_variants = kind in ('both', 'variant')
        want_genes = kind in ('both', 'gene')

        if len(query) < 2:
            return {'query': query, 'variants': [], 'genes': []}

        catalogue = available()
        if species not in catalogue['species']:
            return {'message': 'Species not found'}, 404

        needle = f'%{query}%'
        variants, genes = [], []

        for locus in loci_for(species):
            session = qtl_session(species, locus)
            if session is None:
                continue

            # a variant: how strong its best hit was, and against which gene
            rows = [] if not want_variants else (
                session.query(Variant.variant, Variant.pos, Variant.gene, Variant.feature,
                              func.max(UsageAssociation.neglog10_p), Asc.asc,
                              cast(UsageAssociation.significant, Integer))
                .join(UsageAssociation, UsageAssociation.variant_id == Variant.id)
                .join(Asc, Asc.id == UsageAssociation.asc_id)
                .filter(Variant.variant.like(needle))
                .group_by(Variant.id)
                .order_by(func.max(UsageAssociation.neglog10_p).desc())
                .limit(25)
                .all()
            )
            variants.extend(
                {'locus': locus, 'variant': variant, 'pos': pos, 'gene': gene,
                 'feature': feature, 'neglog10_p': best, 'asc': asc,
                 'significant': bool(significant)}
                for variant, pos, gene, feature, best, asc, significant in rows)

            # ASC names drop the locus - IGL's V9-49 is the gene IGLV9-49 - and
            # people type the full name, so match against both. There are only a
            # few dozen per locus, so filtering here beats string-building in SQL.
            rows = [] if not want_genes else (
                session.query(Asc.asc, Asc.segment,
                              func.count(UsageAssociation.id),
                              func.sum(cast(UsageAssociation.significant, Integer)),
                              func.max(UsageAssociation.neglog10_p))
                .join(UsageAssociation, UsageAssociation.asc_id == Asc.id)
                .group_by(Asc.id)
                .all()
            )
            wanted = query.upper()
            for asc, segment, tested, significant, best in rows:
                if not any(wanted in name.upper() for name in gene_names(locus, asc)):
                    continue
                genes.append(
                    {'locus': locus, 'asc': asc,
                     'gene': display_name(locus, asc), 'segment': segment,
                     'n_variants': tested, 'n_significant': int(significant or 0),
                     'best_neglog10_p': best})

        variants.sort(key=lambda v: v['neglog10_p'], reverse=True)
        genes.sort(key=lambda g: g['n_significant'], reverse=True)

        return {'query': query, 'variants': variants[:25], 'genes': genes[:25]}
