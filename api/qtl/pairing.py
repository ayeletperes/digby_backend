"""What a variant does to which partner a gene recombines with.

The usage scan asks whether a variant changes how much a gene is used. This asks
something else: whether it changes the *company that gene keeps* - given a D, how
the J it joins to is distributed, and the other way round. A variant can leave
every marginal usage untouched and still move the pairing, so the two scans are
not restatements of each other.

CONDITIONALS. `P(J|D)` and `P(D|J)` are two separate scans over the same
variants, and the anchor is whichever side is being conditioned on: for `P(J|D)`
the anchor is the J and the partners are the D genes, for `P(D|J)` it is the
other way round. They are never mixed - `conditional` filters everything here,
because a variant appears in both and pooling them double-reports it.

NAMES. IGH's D clusters are stored with the locus on them (`IGHD5-12`, and
`IGHD4-11/IGHD4-4` when merged) while J is stored bare (`J1`). Those strings are
the join keys between `qtl_dj_enrichment`, `qtl_cell_test` and `qtl_asc`, so they
travel verbatim; shortening for display is the client's business.

STARS. `marked` and `marked_strict` are the pipeline's own call on a cell, kept
rather than re-derived so the dashboard and the manuscript figures agree. A cell
is only meaningful under an omnibus the scan already called significant, which
`omnibus_significant` records; that is passed through rather than recomputed.
"""

import numpy as np
from flask import request
from flask_restx import Resource
from sqlalchemy import case, func

from api.restx import api
from api.system.system import digby_protected
from api.qtl.qtl import _thresholds, check_species_locus, qtl_session
from db.qtl_model import (
    Asc, CellTest, DjEnrichment, Dosage, PairingAssociation, UsageAssociation, Variant,
)

# Its own module, mounted on the guQTL path: one more question about a guQTL
# variant, so it belongs beside the endpoints that answer the others.
ns = api.namespace('qtl_pairing', path='/qtl',
                   description='D/J and V/J pairing under a guQTL variant')

CONDITIONALS = ('P(J|D)', 'P(D|J)')

# Genotype classes, in the order a reader expects to see them.
GENOTYPES = (0, 1, 2)


def _conditional():
    """The requested conditional, or None if it is not one this data holds.

    No default: the two scans answer different questions over the same variants,
    and silently picking one would make the answer depend on a spelling mistake.
    """
    asked = request.args.get('conditional')
    return asked if asked in CONDITIONALS else None


def _box(values):
    """Five-number summary Plotly can draw without being handed every point.

    Quartiles by linear interpolation, which is numpy's default and R's type 7 -
    the same convention the manuscript figures were drawn with, so a box here and
    a box there are the same box. Fences are the data's own extremes rather than
    1.5 IQR: these are proportions from a bounded scale and a subject outside the
    whiskers is a real subject, not a candidate for hiding.
    """
    if not values:
        return None
    a = np.asarray(values, dtype=float)
    q1, med, q3 = (float(x) for x in np.quantile(a, [0.25, 0.5, 0.75]))
    return {'n': int(a.size), 'min': float(a.min()), 'q1': q1, 'median': med,
            'q3': q3, 'max': float(a.max()), 'mean': float(a.mean())}


@ns.route('/pairing_variants/<string:species>/<string:locus>')
@api.response(404, 'Species or locus not found')
class QtlPairingVariantsApi(Resource):
    @digby_protected()
    def get(self, species, locus):
        """ Returns the variants that have a pairing scan, strongest first """

        error = check_species_locus(species, locus)
        if error:
            return error

        session = qtl_session(species, locus)
        conditional = _conditional()
        if conditional is None:
            return {'message': f'conditional must be one of {", ".join(CONDITIONALS)}'}, 400

        # aggregated in SQL rather than by capping the rows and folding in
        # Python: a row cap drops whole variants, and which ones it drops depends
        # on how many anchors each happens to have
        limit = min(int(request.args.get('limit') or 200), 1000)
        summary = (session.query(
                       PairingAssociation.variant_id.label('vid'),
                       func.max(PairingAssociation.neglog10_p).label('best'),
                       # list-of-whens form: this SQLAlchemy predates the 2.0
                       # positional signature
                       func.sum(case([(PairingAssociation.significant == True, 1)],
                                     else_=0)).label('n_sig'),
                       func.count().label('n_anchors'))
                   .filter(PairingAssociation.conditional == conditional)
                   .group_by(PairingAssociation.variant_id)
                   .subquery())

        rows = (session.query(Variant, summary.c.best, summary.c.n_sig,
                              summary.c.n_anchors, PairingAssociation.anchor_gene)
                .join(summary, summary.c.vid == Variant.id)
                # the anchor named is the one that produced the best p, which is
                # the only anchor the summary line can honestly claim
                .join(PairingAssociation,
                      (PairingAssociation.variant_id == Variant.id)
                      & (PairingAssociation.conditional == conditional)
                      & (PairingAssociation.neglog10_p == summary.c.best))
                .order_by(summary.c.best.desc())
                .limit(limit).all())

        total = session.query(func.count(func.distinct(PairingAssociation.variant_id))) \
            .filter(PairingAssociation.conditional == conditional).scalar()

        variants = [{'variant': v.variant, 'pos': v.pos, 'contig': v.contig,
                     'gene': v.gene, 'feature': v.feature, 'maf': v.maf,
                     'anchor_gene': anchor, 'neglog10_p': best,
                     'n_significant_anchors': int(n_sig or 0),
                     'n_anchors': int(n_anchors or 0),
                     'significant': bool(n_sig)}
                    for v, best, n_sig, n_anchors, anchor in rows]

        return {'species': species, 'locus': locus, 'conditional': conditional,
                'conditionals': list(CONDITIONALS),
                # zero scanned and zero significant are different answers, and an
                # empty list cannot tell them apart: only IGH has a pairing scan,
                # and a light chain returning nothing must not read as a scan that
                # found nothing
                'scanned': bool(total),
                # every count it took to get here, so a shortened list cannot be
                # mistaken for the whole scan
                'n_variants_scanned': int(total or 0),
                'n_returned': len(variants),
                'limit': limit,
                'variants': variants}


@ns.route('/pairing/<string:species>/<string:locus>/<path:variant>')
@api.response(404, 'Species, locus or variant not found')
class QtlPairingApi(Resource):
    @digby_protected()
    def get(self, species, locus, variant):
        """ Returns the partner distribution of every gene, split by genotype

        One cell is one (D, J) pair under one genotype class: the distribution
        across subjects of P(J|D) or P(D|J), summarised as a box. The anchor is
        whichever side the conditional conditions on.
        """

        error = check_species_locus(species, locus)
        if error:
            return error

        conditional = _conditional()
        if conditional is None:
            return {'message': f'conditional must be one of {", ".join(CONDITIONALS)}'}, 400

        session = qtl_session(species, locus)
        record = session.query(Variant).filter(Variant.variant == variant).one_or_none()
        if record is None:
            return {'message': f'No such variant: {variant}'}, 404

        # genotype per subject. A subject with no call is not a fourth class - it
        # is absent from every box, and the counts below say how many that is.
        dosages = (session.query(Dosage.subject_id, Dosage.dosage)
                   .filter(Dosage.variant_id == record.id).all())
        genotype = {sid: int(d) for sid, d in dosages if d is not None}
        if not genotype:
            return {'message': f'{variant} has no genotypes, so it has no pairing'}, 404

        enrichment = session.query(DjEnrichment).all()
        if not enrichment:
            return {'message': f'{locus} holds no pairing data'}, 404

        by_j = conditional == 'P(J|D)'
        anchor_of = (lambda e: e.j_gene) if by_j else (lambda e: e.d_gene)
        partner_of = (lambda e: e.d_gene) if by_j else (lambda e: e.j_gene)
        value_of = (lambda e: e.p_j_given_d) if by_j else (lambda e: e.p_d_given_j)
        # the anchor's own marginal: P(J) when conditioning on D, P(D) otherwise
        anchor_marginal = (lambda e: e.p_j) if by_j else (lambda e: e.p_d)
        partner_marginal = (lambda e: e.p_d) if by_j else (lambda e: e.p_j)

        cells = {}
        anchor_values = {}
        partner_values = {}
        skipped = 0

        for e in enrichment:
            gt = genotype.get(e.subject_id)
            if gt is None:
                skipped += 1
                continue
            a, p = anchor_of(e), partner_of(e)
            value = value_of(e)
            if value is not None:
                cells.setdefault((a, p, gt), []).append(value)
            # a marginal is one number per subject per gene, repeated across every
            # partner row, so it is collected against the subject and not appended
            if anchor_marginal(e) is not None:
                anchor_values.setdefault((a, gt), {})[e.subject_id] = anchor_marginal(e)
            if partner_marginal(e) is not None:
                partner_values.setdefault((p, gt), {})[e.subject_id] = partner_marginal(e)

        anchors = sorted({a for a, _, _ in cells})
        partners = sorted({p for _, p, _ in cells})

        # the pipeline's own stars, and the omnibus they sit under
        tests = (session.query(CellTest)
                 .filter(CellTest.conditional == conditional,
                         CellTest.variant_id == record.id).all())
        marks = {}
        for t in tests:
            a, p = (t.j_gene, t.d_gene) if by_j else (t.d_gene, t.j_gene)
            marks[(a, p)] = {
                'marked': bool(t.marked), 'marked_strict': bool(t.marked_strict),
                'p_value': t.p_value, 'delta_mean': t.delta_mean, 'beta': t.beta,
                'n_low': t.n_low, 'n_high': t.n_high,
                'omnibus_significant': bool(t.omnibus_significant),
            }

        omnibus = {a.anchor_gene: {
                       'n': a.n, 'pillai': a.pillai, 'f_stat': a.f_stat,
                       'p_value': a.p_value, 'neglog10_p': a.neglog10_p,
                       'min_genotype_group': a.min_genotype_group,
                       'significant': bool(a.significant)}
                   for a in session.query(PairingAssociation)
                   .filter(PairingAssociation.conditional == conditional,
                           PairingAssociation.variant_id == record.id).all()}

        segment = {row.asc: row.segment for row in session.query(Asc).all()}

        # The marginal panels are the usage scan's question, not this one's: does
        # the variant change how MUCH a gene is used, as against who it pairs
        # with. Passed through so the two can be read against each other, which
        # is the whole reason the marginals are drawn beside the grid.
        usage = {
            asc: {'beta': beta, 'p_value': p_value, 'neglog10_p': neglog10_p,
                  'significant': bool(sig), 'n': record.n,
                  'min_genotype_group': record.min_genotype_group}
            for asc, beta, p_value, neglog10_p, sig in
            session.query(Asc.asc, UsageAssociation.beta, UsageAssociation.p_value,
                          UsageAssociation.neglog10_p, UsageAssociation.significant)
            .select_from(UsageAssociation)
            .join(Asc, Asc.id == UsageAssociation.asc_id)
            .filter(UsageAssociation.variant_id == record.id).all()}

        return {
            'species': species,
            'locus': locus,
            'conditional': conditional,
            'conditionals': list(CONDITIONALS),
            'anchor_side': 'J' if by_j else 'D',
            'variant': {'variant': record.variant, 'contig': record.contig,
                        'pos': record.pos, 'maf': record.maf, 'gene': record.gene,
                        'feature': record.feature, 'sub_feature': record.sub_feature},
            'genotypes': [{'genotype': g,
                           'n': sum(1 for v in genotype.values() if v == g)}
                          for g in GENOTYPES],
            # every filter says what it dropped: enrichment rows for subjects this
            # variant was not called in are not in any box, and a silent drop here
            # would look like a smaller cohort rather than a missing genotype
            'subjects': {'genotyped': len(genotype),
                         'enrichment_rows_without_genotype': skipped},
            'anchors': anchors,
            'partners': partners,
            'segments': segment,
            # keyed by the stored ASC name, the same key the anchors and partners
            # use, so no name has to be rebuilt to line the two scans up
            'usage': usage,
            # the run's own thresholds. The usage flag above is set against the
            # study-wide corrected one, not against 0.05, and the panel has to be
            # able to say so rather than let a cell-level star key stand for it.
            'thresholds': _thresholds(session),
            'omnibus': omnibus,
            'cells': [{'anchor': a, 'partner': p, 'genotype': g, 'box': _box(values)}
                      for (a, p, g), values in sorted(cells.items())],
            # a star belongs to the pair, not to one of its three genotype boxes;
            # kept apart so a client cannot draw it three times
            'marks': [{'anchor': a, 'partner': p, **mark}
                      for (a, p), mark in sorted(marks.items())],
            'anchor_marginal': [{'gene': a, 'genotype': g,
                                 'box': _box(list(per_subject.values()))}
                                for (a, g), per_subject in sorted(anchor_values.items())],
            'partner_marginal': [{'gene': p, 'genotype': g,
                                  'box': _box(list(per_subject.values()))}
                                 for (p, g), per_subject in sorted(partner_values.items())],
        }
