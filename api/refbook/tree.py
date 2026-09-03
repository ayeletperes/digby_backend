# Hierarchical clustering of the alleles of one ASC. Live and memoised: the worst

import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from flask_restx import Resource
from functools import lru_cache

from api.restx import api
from api.system.system import digby_protected
from api.refbook.refbook import (check_species_locus, collect_asc_sequences, dataset_stamp,
                                 requested_sources, requested_list, segment_of)

ns = api.namespace('refbook_tree', description='Hierarchical clustering of an ASC',
                   path='/refbook')

GAP = ord('.')
PAD = ord(' ')


def _matrix(seqs):
    """ Sequences as a (n, width) uint8 matrix, right-padded with spaces. """
    width = max(len(s) for s in seqs)
    return np.stack([np.frombuffer(s.encode().ljust(width, b' '), dtype=np.uint8)
                     for s in seqs])


def _extent(M):
    """ First and last real column per row. Terminal gaps are padding, not """
    real = (M != GAP) & (M != PAD)
    cols = np.arange(M.shape[1])
    start = np.where(real.any(1), np.argmax(real, 1), 0)
    end = np.where(real.any(1), M.shape[1] - 1 - np.argmax(real[:, ::-1], 1), -1)
    return start, end, cols


def regap(seqs):
    """ Put every sequence back into the gene's IMGT column frame. """
    widths = [len(s) for s in seqs]
    frame_width = max(set(widths), key=widths.count)
    if len(set(widths)) == 1:
        return seqs, [0] * len(seqs)

    from Bio import Align
    aligner = Align.PairwiseAligner(mode='global', match_score=1, mismatch_score=-1,
                                    open_gap_score=-2, extend_gap_score=-1)
    frame = seqs[widths.index(frame_width)]

    out, dropped = [], []
    for seq, width in zip(seqs, widths):
        if width == frame_width:
            out.append(seq)
            dropped.append(0)
            continue

        row = ['.'] * frame_width
        blocks = aligner.align(frame, seq)[0].aligned
        for (t0, t1), (q0, q1) in zip(*blocks):
            row[t0:t1] = seq[q0:q1]
        out.append(''.join(row))
        dropped.append(width - sum(int(q1 - q0) for q0, q1 in blocks[1]))

    return out, dropped


def gapped_distances(seqs):
    """ Difference counts over IMGT-gapped sequences, plus the informative columns. """
    seqs, dropped = regap(seqs)

    M = _matrix(seqs)
    start, end, cols = _extent(M)

    lo = np.maximum(start[:, None], start[None, :])[:, :, None]
    hi = np.minimum(end[:, None], end[None, :])[:, :, None]
    valid = (cols >= lo) & (cols <= hi)

    dist = ((M[:, None, :] != M[None, :, :]) & valid).sum(2).astype(float)

    # an insertion has no column in the frame, so it is counted by its length: two
    drop = np.array(dropped, dtype=float)
    dist += np.abs(drop[:, None] - drop[None, :])

    # a column is informative if the alleles that cover it do not all agree
    covered = (cols >= start[:, None]) & (cols <= end[:, None])
    informative = [int(c) for c in cols
                   if len(np.unique(M[covered[:, c], c])) > 1]

    return dist, informative, int(M.shape[1])


def nw_distances(seqs):
    """ Edit distance for D and J, which have no gapped form. """
    from Bio import Align

    aligner = Align.PairwiseAligner(mode='global', match_score=0, mismatch_score=-1,
                                    open_gap_score=-1, extend_gap_score=-1)

    n = len(seqs)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dist[i, j] = dist[j, i] = -aligner.score(seqs[i], seqs[j])
    return dist, [], max(len(s) for s in seqs)


def linkage(dist):
    """ Complete-linkage clustering, scipy's format: [i, j, height, size]. """
    return [[float(a), float(b), float(h), float(n)]
            for a, b, h, n in hierarchy.linkage(squareform(dist, checks=False),
                                                method='complete')]


def leaf_order(merges, n):
    """ Leaves left to right, as the dendrogram draws them. """
    if not merges:
        return list(range(n))
    return [int(i) for i in hierarchy.leaves_list(np.array(merges, dtype=float))]


def duplicate_groups(names, seqs):
    """ Alleles sharing an identical sequence. Reported, not collapsed: they stay """
    groups = {}
    for name, seq in zip(names, seqs):
        groups.setdefault(seq, []).append(name)
    return [g for g in groups.values() if len(g) > 1]


# large enough that sweeping one locus does not evict before reuse
@lru_cache(maxsize=4096)
def _tree(species, locus, asc, stamp, sources, allele_names):
    """ The tree for one ASC. `stamp` only keys the cache: see dataset_stamp. """
    recs = collect_asc_sequences(species, locus, asc, set(sources),
                                 list(allele_names) if allele_names else None)
    if not recs:
        return None

    recs.sort(key=lambda r: r['name'])
    names = [r['name'] for r in recs]
    segment = segment_of(asc)

    # V is stored IMGT-gapped and therefore already column-aligned; D and J are not
    gapped = segment == 'V' and all(r['seq_gapped'] for r in recs)
    seqs = [r['seq_gapped'] if gapped else r['seq'] for r in recs]

    if len(recs) == 1:
        return {'labels': names, 'merges': [], 'order': [0], 'informative_columns': [],
                'duplicate_groups': [], 'columns': len(seqs[0]), 'gapped': gapped,
                'metric': 'hamming_gapped' if gapped else 'edit'}

    dist, informative, columns = (gapped_distances(seqs) if gapped else nw_distances(seqs))
    merges = linkage(dist)

    return {'labels': names,
            'merges': [[int(a), int(b), h, int(s)] for a, b, h, s in merges],
            'order': leaf_order(merges, len(names)),
            'informative_columns': informative,
            'duplicate_groups': duplicate_groups(names, seqs),
            'columns': columns,
            'gapped': gapped,
            'metric': 'hamming_gapped' if gapped else 'edit'}


@ns.route('/asc_tree/<string:species>/<string:locus>/<path:asc>')
@api.response(404, 'Species or locus not found')
class AscTree(Resource):
    @digby_protected()
    def get(self, species, locus, asc):
        """ Hierarchical clustering of every allele in an ASC """

        error = check_species_locus(species, locus)
        if error:
            return error

        alleles = requested_list('alleles')
        tree = _tree(species, locus, asc, dataset_stamp(species, locus),
                     frozenset(requested_sources()),
                     frozenset(alleles) if alleles else None)

        if tree is None:
            return {'message': f'No sequences for {asc}'}, 404

        info = _tree.cache_info()
        return {'asc': asc, 'segment': segment_of(asc), 'linkage': 'complete', **tree,
                'cache': {'hits': info.hits, 'misses': info.misses, 'size': info.currsize}}
