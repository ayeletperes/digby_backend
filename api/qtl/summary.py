"""Where the significant gene-usage variants are, and what they sit in.

The Manhattan answers this one scan at a time. This answers it for the whole
run: how many variants reached significance, in which segment's genes, and
whether they are coding changes, leader, RSS, UTR or intergenic.

COUNTING. A variant is counted once per (segment, feature), not once per
association. A variant significant for six V ASCs is one V variant, not six -
`count(distinct variant)`, because the question is how many *sites* matter, and
counting associations would rank a locus by how many genes it happens to have.
A variant significant against both a V and a J ASC is counted under both, since
it genuinely did both; the segment totals therefore need not sum to the locus
total, and the response says so.

DENOMINATORS. Every count travels with what it was drawn from: the locus totals
here, and `n_tested` per gene below, so a count can be read as a rate rather than
a raw height. There is deliberately no per-segment denominator - in all three
loci every gene is tested against every variant, so it was the locus total under
another name and cost a 658,000-row distinct pass to say so.

UTR. The manuscript figure folds `utr` into `intergenic`; the database keeps them
apart, so this does too. The two IGL variants involved are the whole difference
between the counts here and the counts in the figure, and folding them silently
would leave that unexplained.
"""

from flask import request
from flask_restx import Resource
from sqlalchemy import distinct, func

from api.restx import api
from api.system.system import digby_protected
from api.qtl.qtl import (available, check_species_locus, current_project,
                          loci_for, qtl_session)
from db.qtl_model import Asc, UsageAssociation, Variant

ns = api.namespace('qtl_summary', path='/qtl',
                   description='Where the significant guQTL variants are')

# The order these read in, which is the order the locus is laid out in.
SEGMENT_ORDER = ['V', 'D', 'J', 'C']
FEATURE_ORDER = ['coding', 'leader', 'rss', 'utr', 'intergenic']


def _rank(values, order):
    """Known things in their conventional order, unknown ones after, sorted."""
    known = [v for v in order if v in values]
    return known + sorted(v for v in values if v not in order)


def _locus_counts(session):
    """Distinct significant variants by segment and feature."""
    return (session.query(Asc.segment, Variant.feature,
                          func.count(distinct(Variant.id)))
            .select_from(UsageAssociation)
            .join(Variant, Variant.id == UsageAssociation.variant_id)
            .join(Asc, Asc.id == UsageAssociation.asc_id)
            .filter(UsageAssociation.significant == True)
            .group_by(Asc.segment, Variant.feature).all())


def _tested_per_asc(session):
    """How many variants each gene was tested against.

    Two things keep this off a full pass over 658,000 rows. `count(*)` rather
    than `count(distinct variant_id)`: the table is unique on (variant_id,
    asc_id), so within one gene the two are the same number - checked, and they
    agree on every row. And grouped by `asc_id` rather than by the gene's name,
    which lets `ix_usage_asc_p` do the grouping instead of joining every
    association to `qtl_asc` first; the 70 names are looked up afterwards.
    """
    counts = dict(session.query(UsageAssociation.asc_id, func.count())
                  .group_by(UsageAssociation.asc_id).all())
    return {name: int(counts.get(asc_id, 0))
            for asc_id, name in session.query(Asc.id, Asc.asc).all()}


@ns.route('/usage_summary/<string:species>')
@api.response(404, 'Species not found')
class QtlUsageSummaryApi(Resource):
    @digby_protected()
    def get(self, species):
        """ Returns significant gene-usage variants by locus, segment and location """

        catalogue = available()
        if species not in catalogue['species']:
            return {'message': 'Species not found'}, 404

        # one project's loci, not every locus any project holds: this rolls three
        # scans into one figure, and a figure that mixed cohorts would be reading
        # across studies the analysis never pooled
        project = current_project()
        loci = loci_for(species, project)

        rows = []
        totals = {}
        segments, features = set(), set()

        for locus in loci:
            session = qtl_session(species, locus)
            if session is None:
                continue

            for segment, feature, n in _locus_counts(session):
                rows.append({'locus': locus, 'segment': segment,
                             'feature': feature, 'n': int(n)})
                segments.add(segment)
                features.add(feature)

            significant = (session.query(func.count(distinct(Variant.id)))
                           .select_from(UsageAssociation)
                           .join(Variant, Variant.id == UsageAssociation.variant_id)
                           .filter(UsageAssociation.significant == True).scalar())
            totals[locus] = {
                'n_variants': int(session.query(func.count(Variant.id)).scalar() or 0),
                'n_significant': int(significant or 0),
            }

        return {
            'species': species,
            'project': project,
            'loci': loci,
            'segments': _rank(segments, SEGMENT_ORDER),
            'features': _rank(features, FEATURE_ORDER),
            'rows': rows,
            'totals': totals,
            'counting': 'distinct variants; a variant significant for several ASCs '
                        'of one segment counts once, and one significant in two '
                        'segments counts in both, so segment totals need not sum '
                        'to the locus total',
        }


@ns.route('/gene_summary/<string:species>/<string:locus>')
@api.response(404, 'Species or locus not found')
class QtlGeneSummaryApi(Resource):
    @digby_protected()
    def get(self, species, locus):
        """ Returns per-gene significant variant counts, by where the variants sit

        One row per ASC. `by_feature` counts that gene's own significant
        variants, so unlike the locus summary these do sum: a variant significant
        for two genes appears under each of them.
        """

        error = check_species_locus(species, locus)
        if error:
            return error

        session = qtl_session(species, locus)

        counts = (session.query(Asc.asc, Asc.segment, Variant.feature,
                                func.count(distinct(Variant.id)))
                  .select_from(UsageAssociation)
                  .join(Variant, Variant.id == UsageAssociation.variant_id)
                  .join(Asc, Asc.id == UsageAssociation.asc_id)
                  .filter(UsageAssociation.significant == True)
                  .group_by(Asc.asc, Asc.segment, Variant.feature).all())

        tested = _tested_per_asc(session)

        # where the gene's strongest significant variant sits, so a row can be
        # placed on the locus as well as counted
        leads = (session.query(Asc.asc, Variant.variant, Variant.pos, Variant.contig,
                               func.max(UsageAssociation.neglog10_p))
                 .select_from(UsageAssociation)
                 .join(Variant, Variant.id == UsageAssociation.variant_id)
                 .join(Asc, Asc.id == UsageAssociation.asc_id)
                 .filter(UsageAssociation.significant == True)
                 .group_by(Asc.asc).all())
        lead = {asc: {'variant': v, 'pos': pos, 'contig': contig, 'neglog10_p': p}
                for asc, v, pos, contig, p in leads}

        genes = {}
        features = set()
        for asc, segment, feature, n in counts:
            entry = genes.setdefault(asc, {'asc': asc, 'segment': segment,
                                           'by_feature': {}, 'n_significant': 0})
            entry['by_feature'][feature] = int(n)
            entry['n_significant'] += int(n)
            features.add(feature)

        # every scanned gene, including the ones with nothing: a list of only the
        # genes that hit reads as though the rest were not looked at
        for row in session.query(Asc).all():
            entry = genes.setdefault(row.asc, {'asc': row.asc, 'segment': row.segment,
                                               'by_feature': {}, 'n_significant': 0})
            entry['n_tested'] = int(tested.get(row.asc, 0))
            entry['lead'] = lead.get(row.asc)

        return {
            'species': species,
            'locus': locus,
            'features': _rank(features, FEATURE_ORDER),
            'segments': _rank({g['segment'] for g in genes.values()}, SEGMENT_ORDER),
            'genes': sorted(genes.values(), key=lambda g: -g['n_significant']),
            'n_genes_scanned': len(genes),
            'n_genes_with_signal': sum(1 for g in genes.values() if g['n_significant']),
        }
