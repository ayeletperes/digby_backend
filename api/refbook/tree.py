# Hierarchical clustering of the alleles of one ASC, for the refbook dashboard.
#
# Computed live and memoised, not precomputed: the worst gene in the dataset
# (Human/IGL/IGLV3-1, 133 alleles) costs ~60ms, and the panel accepts a `sources`
# and an `alleles` subset, which a precomputed table cannot enumerate.

import numpy as np
from flask import request
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

LINKAGES = ('complete', 'single', 'average')


def _matrix(seqs):
    """ Sequences as a (n, width) uint8 matrix, right-padded with spaces. """
    width = max(len(s) for s in seqs)
    return np.stack([np.frombuffer(s.encode().ljust(width, b' '), dtype=np.uint8)
                     for s in seqs])


def _extent(M):
    """ First and last column of real sequence in each row.

    Terminal gaps are alignment padding, not a difference: an allele recorded
    only from position 20 onwards must not read as 19 substitutions against a
    full-length one. Internal gaps are kept, since those are real indels.
    """
    real = (M != GAP) & (M != PAD)
    cols = np.arange(M.shape[1])
    start = np.where(real.any(1), np.argmax(real, 1), 0)
    end = np.where(real.any(1), M.shape[1] - 1 - np.argmax(real[:, ::-1], 1), -1)
    return start, end, cols


def regap(seqs):
    """ Put every sequence back into the gene's own IMGT column frame.

    Almost all gapped sequences for a gene are the same width and their columns
    already correspond. The exceptions are novel alleles carrying an insertion or a
    deletion, which shift every column after it: 62 of the 589 multi-allele V genes
    have at least one, and in IGHV4-39 (12 of 42 alleles) it put the root of the
    tree at 188 differences instead of ~25. Those rows are realigned to the widest
    common frame and rebuilt on its coordinates, so the fast column comparison holds
    for the whole gene.

    Returns the reframed sequences and, per sequence, the number of inserted columns
    dropped, which the frame cannot represent and which the caller adds back.
    """
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
    """ Pairwise difference counts over IMGT-gapped sequences, and the informative columns.

    The database stores V alleles already IMGT-gapped, so the columns correspond by
    construction and comparing them costs ~40ms where re-aligning every pair of the
    same gene with a pairwise aligner costs ~40s.
    """
    seqs, dropped = regap(seqs)

    M = _matrix(seqs)
    start, end, cols = _extent(M)

    lo = np.maximum(start[:, None], start[None, :])[:, :, None]
    hi = np.minimum(end[:, None], end[None, :])[:, :, None]
    valid = (cols >= lo) & (cols <= hi)

    dist = ((M[:, None, :] != M[None, :, :]) & valid).sum(2).astype(float)

    # an insertion has no column in the frame, so it is counted by its length: two
    # alleles sharing one are still identical, and one against an allele without it
    # differs by its length
    drop = np.array(dropped, dtype=float)
    dist += np.abs(drop[:, None] - drop[None, :])

    # a column is informative if the alleles that cover it do not all agree
    covered = (cols >= start[:, None]) & (cols <= end[:, None])
    informative = [int(c) for c in cols
                   if len(np.unique(M[covered[:, c], c])) > 1]

    return dist, informative, int(M.shape[1])


def nw_distances(seqs):
    """ Pairwise distances for D and J, which have no gapped form and ragged ends.

    Needleman-Wunsch scored so that minus the alignment score is the edit distance.
    n is at most ~10 here, so all pairs cost about a millisecond.

    End gaps are penalised, unlike the V path: there the IMGT columns say a terminal
    gap is sequence that was never recorded, whereas here nothing distinguishes that
    from a genuinely shorter allele - and free end gaps make the empty overlap free,
    which scores every pair as identical.
    """
    from Bio import Align

    aligner = Align.PairwiseAligner(mode='global', match_score=0, mismatch_score=-1,
                                    open_gap_score=-1, extend_gap_score=-1)

    n = len(seqs)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dist[i, j] = dist[j, i] = -aligner.score(seqs[i], seqs[j])
    return dist, [], max(len(s) for s in seqs)


def linkage(dist, method='complete'):
    """ Agglomerative clustering, in scipy's linkage format.

    scipy is not in the backend environment and this is the only thing it would be
    used for, so the ~25 lines are here instead. Rows are [i, j, height, size],
    leaves are 0..n-1 and merge k creates cluster n+k, as scipy does it.
    """
    n = dist.shape[0]
    d = dist.astype(float).copy()
    np.fill_diagonal(d, np.inf)

    cid = list(range(n))            # cluster id currently held in each row
    size = [1] * n
    merges = []

    for k in range(n - 1):
        i, j = np.unravel_index(np.argmin(d), d.shape)
        height = float(d[i, j])

        if method == 'single':
            row = np.minimum(d[i], d[j])
        elif method == 'average':
            row = (size[i] * d[i] + size[j] * d[j]) / (size[i] + size[j])
        else:
            row = np.maximum(d[i], d[j])

        merges.append([min(cid[i], cid[j]), max(cid[i], cid[j]), height, size[i] + size[j]])

        d[i] = d[:, i] = row        # row i becomes the merged cluster
        d[i, i] = np.inf
        d[j] = d[:, j] = np.inf     # row j leaves the pool
        cid[i] = n + k
        size[i] += size[j]

    return merges


def leaf_order(merges, n):
    """ Leaves left to right, by walking the tree from its root. """
    if not merges:
        return list(range(n))

    order, stack = [], [n + len(merges) - 1]
    while stack:
        node = stack.pop()
        if node < n:
            order.append(node)
        else:
            left, right, _, _ = merges[node - n]
            stack.extend([right, left])         # left comes off the stack first
    return order


def duplicate_groups(names, seqs):
    """ Alleles that share an identical sequence, reported rather than collapsed.

    26 of the 672 multi-allele genes have at least two alleles with the same gapped
    sequence. They stay as separate leaves joined at height 0, so no allele silently
    disappears from the tree.
    """
    groups = {}
    for name, seq in zip(names, seqs):
        groups.setdefault(seq, []).append(name)
    return [g for g in groups.values() if len(g) > 1]


# 1,138 genes today, and a gene can be asked for under several source or allele
# filters. Sized like the alignment cache: sweeping one large locus through a
# smaller LRU evicts each entry before it is reused. A tree is a few KB.
@lru_cache(maxsize=4096)
def _tree(species, locus, asc, stamp, sources, allele_names, method):
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
    merges = linkage(dist, method)

    return {'labels': names,
            'merges': [[int(a), int(b), h, int(s)] for a, b, h, s in merges],
            'order': leaf_order(merges, len(names)),
            'informative_columns': informative,
            'duplicate_groups': duplicate_groups(names, seqs),
            'columns': columns,
            'gapped': gapped,
            # two different measures, and the caller has to say which: over the
            # IMGT columns it is a count of differing positions, and off them it
            # is an edit distance that can also count a gap
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

        method = request.args.get('linkage', 'complete').lower()
        if method not in LINKAGES:
            return {'message': f'linkage must be one of {", ".join(LINKAGES)}'}, 400

        alleles = requested_list('alleles')
        tree = _tree(species, locus, asc, dataset_stamp(species, locus),
                     frozenset(requested_sources()),
                     frozenset(alleles) if alleles else None, method)

        if tree is None:
            return {'message': f'No sequences for {asc}'}, 404

        info = _tree.cache_info()
        return {'asc': asc, 'segment': segment_of(asc), 'linkage': method, **tree,
                'cache': {'hits': info.hits, 'misses': info.misses, 'size': info.currsize}}


def demo():
    """ Self-check: python -c 'import app; from api.refbook.tree import demo; demo()'

    (imported through `app` because api.restx imports it back, so this module
    cannot be the one that starts the chain)
    """
    # terminal gaps are padding, so 0 and 1 are the same sequence
    dist, informative, width = gapped_distances(['..ACGT..', 'AAACGT..', '..ACTT..', '..ATTA..'])
    assert width == 8
    assert dist[0, 1] == 0, dist
    assert dist[0, 2] == 1 and dist[0, 3] == 3 and dist[2, 3] == 2, dist
    assert informative == [3, 4, 5], informative    # 0-1 are covered by one allele only

    # complete linkage: {0,1} at 0, 2 joins at max(1,1), 3 at max(3,3,2)
    merges = linkage(dist)
    assert merges == [[0, 1, 0.0, 2], [2, 4, 1.0, 3], [3, 5, 3.0, 4]], merges
    assert linkage(dist, 'single')[-1][2] == 2.0    # single takes the nearest instead
    assert leaf_order(merges, 4) == [3, 2, 0, 1]    # the outlier ends up on one edge
    assert duplicate_groups(['a', 'b', 'c'], ['AC', 'AC', 'AG']) == [['a', 'b']]

    # an insertion shifts every later column, so the row is put back in the frame
    framed, dropped = regap(['ACGTACGTAC', 'ACGTTTACGTAC', 'ACGTACGTAC'])
    assert framed[1] == 'ACGTACGTAC' and dropped == [0, 2, 0], (framed, dropped)
    d2, _, _ = gapped_distances(['ACGTACGTAC', 'ACGTTTACGTAC', 'ACGTACGTAG'])
    assert d2[0, 1] == 2 and d2[0, 2] == 1 and d2[1, 2] == 3, d2

    # D and J: edit distance, so a shift costs two indels and a substitution costs one
    d, _, _ = nw_distances(['GGTATAAC', 'GTATAACT', 'GGTAAAAC'])
    assert d[0, 1] == 2 and d[0, 2] == 1, d
    print('ok')


if __name__ == '__main__':
    demo()
