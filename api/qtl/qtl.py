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
from sqlalchemy import Integer, cast, func

from api.restx import api
from api.system.system import digby_protected
from app import qtl_dbs
from db.qtl_model import (
    Asc, AscUsage, Dosage, Subject, Threshold, UsageAssociation, Variant,
)

ns = api.namespace('qtl', description='Gene-usage QTL results')

# Loci offered by the API. IGH is held back for now: its cohort has no linked
# repertoire in VDJbase, which the analyses that join the two will need. The data
# is built and on disk, so this is a one-line change when that lands.
HIDDEN_LOCI = {'IGH'}


def qtl_session(species, locus):
    """Session for one guQTL dataset, or None if it is absent or held back."""
    if locus in HIDDEN_LOCI:
        return None

    dataset = qtl_dbs.get(species, {}).get(locus)
    return dataset.session if dataset is not None else None


def available():
    """Species and loci that have guQTL results and are offered."""
    ret = {'species': [], 'loci': {}}

    for species in sorted(qtl_dbs):
        loci = sorted(l for l in qtl_dbs[species] if l not in HIDDEN_LOCI)
        if loci:
            ret['species'].append(species)
            ret['loci'][species] = loci

    return ret


def check_species_locus(species, locus):
    """A 404 response if this species/locus has no guQTL results, else None."""
    catalogue = available()

    if species not in catalogue['species']:
        return {'message': 'Species not found'}, 404
    if locus not in catalogue['loci'].get(species, []):
        return {'message': 'Locus not found'}, 404

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

        rows = (
            session.query(Asc.asc, Asc.segment, Asc.n_member,
                          func.count(UsageAssociation.id),
                          # cast first: summing a Boolean column runs the total
                          # back through the Boolean result processor, so 206
                          # arrives as True and counts as 1
                          func.sum(cast(UsageAssociation.significant, Integer)),
                          func.max(UsageAssociation.neglog10_p))
            .join(UsageAssociation, UsageAssociation.asc_id == Asc.id)
            .group_by(Asc.id)
            .all()
        )

        ascs = [{'asc': asc, 'segment': segment, 'n_member': n_member,
                 'n_variants': tested, 'n_significant': int(significant or 0),
                 'best_neglog10_p': best}
                for asc, segment, n_member, tested, significant, best in rows]
        ascs.sort(key=lambda a: (a['segment'] or '', a['asc']))

        return {'ascs': ascs, 'thresholds': _thresholds(session)}


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

        query = (
            session.query(Variant.variant, Variant.pos, Variant.maf, Variant.gene,
                          Variant.feature,
                          func.max(UsageAssociation.neglog10_p),
                          func.max(cast(UsageAssociation.significant, Integer)))
            .join(UsageAssociation, UsageAssociation.variant_id == Variant.id)
            .group_by(Variant.id)
        )

        if asc:
            query = query.join(Asc, Asc.id == UsageAssociation.asc_id).filter(Asc.asc == asc)

        points = [{'variant': variant, 'pos': pos, 'maf': maf, 'gene': gene,
                   'feature': feature, 'neglog10_p': neglog10_p,
                   'significant': bool(significant)}
                  for variant, pos, maf, gene, feature, neglog10_p, significant in query.all()]
        points.sort(key=lambda p: p['pos'] if p['pos'] is not None else 0)

        return {'locus': locus, 'asc': asc, 'points': points,
                'thresholds': _thresholds(session),
                'leads': _leads(session, asc)}


def _leads(session, asc=None, limit=10):
    """The strongest independent signals, for labelling the plot."""
    query = (
        session.query(Variant.variant, Variant.pos, Asc.asc, UsageAssociation.neglog10_p,
                      UsageAssociation.beta, UsageAssociation.min_genotype_group,
                      UsageAssociation.well_powered)
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

        rows = (
            session.query(Asc.asc, Asc.segment, UsageAssociation.beta, UsageAssociation.se,
                          UsageAssociation.p_value, UsageAssociation.neglog10_p,
                          UsageAssociation.significant, UsageAssociation.n,
                          UsageAssociation.min_genotype_group, UsageAssociation.well_powered,
                          UsageAssociation.is_lead)
            .join(UsageAssociation, UsageAssociation.asc_id == Asc.id)
            .filter(UsageAssociation.variant_id == record.id)
            .order_by(UsageAssociation.neglog10_p.desc())
            .all()
        )

        return {
            'variant': {'variant': record.variant, 'contig': record.contig, 'pos': record.pos,
                        'maf': record.maf, 'gene': record.gene, 'feature': record.feature,
                        'sub_feature': record.sub_feature,
                        'distance_to_gene': record.distance_to_gene},
            'associations': [
                {'asc': asc, 'segment': segment, 'beta': beta, 'se': se, 'p_value': p,
                 'neglog10_p': neglog10_p, 'significant': bool(significant), 'n': n,
                 'min_genotype_group': min_group,
                 'well_powered': bool(powered) if powered is not None else None,
                 'is_lead': bool(lead)}
                for asc, segment, beta, se, p, neglog10_p, significant, n,
                    min_group, powered, lead in rows],
            'has_genotypes': session.query(Dosage.id)
                .filter(Dosage.variant_id == record.id).first() is not None,
        }


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
