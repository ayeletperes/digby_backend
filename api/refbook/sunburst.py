""" The sunburst panel's data: one locus as depth-ordered parallel arrays.

The panel draws the whole locus at once and never fetches again, so this is a
single request per species+locus. The hierarchy is

    chain -> gene_type -> subgroup -> asc -> allele

and it is returned flattened, level by level, as parallel arrays. Because the
nodes are emitted in level order, `parent[i] < i` holds for every node, which is
what lets the client derive arc sweeps, depths and subtree membership in linear
forward/backward passes instead of walking parent chains per node.
"""

import re

from flask_restx import Namespace, Resource
from sqlalchemy import or_

from app import vdjbase_dbs, genomic_dbs
from api.restx import api
from api.system.system import digby_protected
from api.refbook.refbook import check_species_locus, dataset_session, requested_sources
from db.vdjbase_model import Gene as VDJbaseGene, Allele as VDJbaseAllele
from db.genomic_db import Gene as GenomicGene, Sequence as GenomicSequence

# a namespace of its own so app.py can register it separately; the path puts it
# alongside the rest of the refbook API
ns = Namespace('refbook_sunburst', description='Locus hierarchy for the sunburst panel',
               path='/refbook')

LEVELS = ['chain', 'gene_type', 'subgroup', 'asc', 'allele']

SEGMENT_ORDER = 'VDJC'

_FAMILY = re.compile(r'^([A-Z]{2,3}[VDJC]\d*)')
_DIGITS = re.compile(r'(\d+)')


_ISOTYPE = re.compile(r'^(IG[HKL][ADEGM])')


def family_of(gene, segment=None):
    """ IGHV1-18 -> IGHV1, IGHG1 -> IGHG.

    Derived from the gene name rather than read off Gene.family, because the two
    databases fill that column differently: the AIRR-seq side holds 'IGHV1' where
    the genomic side holds '1'. Merging on it would split every subgroup in two.

    Constant genes are named for their isotype rather than a numbered subgroup,
    so they group by it: IGHG1 to IGHG4 and IGHG4D are the IgG subclasses and
    belong under IGHG. The segment has to be passed in, because IGHD is the delta
    constant gene and also the prefix every D gene shares.
    """
    if segment == 'C':
        match = _ISOTYPE.match(gene or '')
        return match.group(1) if match else (gene or '?')

    match = _FAMILY.match(gene or '')
    return match.group(1) if match else (gene or '?')


def segment_of(gene_type):
    """ IGHV -> V. The locus prefix is three characters, the segment is what follows. """
    if gene_type and len(gene_type) > 3:
        return gene_type[3:]
    return gene_type or '?'


def _natural(text):
    """ Sort key that reads runs of digits as numbers, so IGHV2 precedes IGHV10. """
    return [int(part) if part.isdigit() else part for part in _DIGITS.split(text or '')]


def _order(record):
    segment, subgroup, asc, allele = record[0], record[1], record[2], record[3]
    rank = SEGMENT_ORDER.find(segment[:1])
    return (rank if rank >= 0 else len(SEGMENT_ORDER),
            _natural(subgroup), _natural(asc), _natural(allele))


def collect(species, locus, sources):
    """ Every allele of the locus, unioned across the databases the caller asked for.

    Returns a list of (segment, subgroup, asc, allele, novel, in_genomic, in_airrseq).

    Pseudogenes and orphons are excluded, matching /ascs_in_locus, and on the
    genomic side only Functional and ORF sequences are taken, matching
    /ascs_overview. An allele held in both databases is one node, flagged as
    present in both.
    """
    found = {}

    def record(allele, gene, gene_type, novel, source):
        if '/OR' in gene:
            return
        entry = found.get(allele)
        if entry is None:
            segment = segment_of(gene_type)
            entry = found[allele] = {'gene': gene, 'segment': segment,
                                     'subgroup': family_of(gene, segment), 'novel': False,
                                     'genomic': 0, 'airrseq': 0}
        entry[source] = 1
        # a novel call in either database is enough to mark the allele novel; the
        # databases are built separately and only one of them may have seen it
        entry['novel'] = entry['novel'] or bool(novel)

    session = dataset_session(vdjbase_dbs, species, locus, sources, 'airrseq')
    if session is not None:
        rows = session.query(VDJbaseAllele.name, VDJbaseAllele.novel,
                             VDJbaseGene.name, VDJbaseGene.type) \
            .join(VDJbaseGene, VDJbaseGene.id == VDJbaseAllele.gene_id) \
            .filter(VDJbaseGene.pseudo_gene == 0).all()
        for allele, novel, gene, gene_type in rows:
            record(allele, gene, gene_type, novel, 'airrseq')

    session = dataset_session(genomic_dbs, species, locus, sources, 'genomic')
    if session is not None:
        rows = session.query(GenomicSequence.name, GenomicSequence.novel,
                             GenomicGene.name, GenomicGene.type) \
            .join(GenomicGene, GenomicGene.id == GenomicSequence.gene_id) \
            .filter(GenomicGene.pseudo_gene == 0,
                    or_(GenomicSequence.functional == 'Functional',
                        GenomicSequence.functional == 'ORF')).all()
        for allele, novel, gene, gene_type in rows:
            record(allele, gene, gene_type, novel, 'genomic')

    return sorted(((rec['segment'], rec['subgroup'], rec['gene'], allele,
                    int(rec['novel']), rec['genomic'], rec['airrseq'])
                   for allele, rec in found.items()), key=_order)


def build(chain, alleles):
    """ Flatten the hierarchy into depth-ordered parallel arrays.

    `alleles` is what collect() returns, already sorted, so each level comes out
    in that order too and the arcs are drawn in gene order rather than at random.

    novel/nG/nA are counts of alleles in the subtree, summed up the tree in one
    backward pass. They are counts of alleles, not of samples: samples cannot be
    added across alleles without counting the same subject many times, whereas an
    allele belongs to exactly one node at every level.
    """
    label = [chain]
    parent = [-1]
    novel, n_genomic, n_airrseq = [0], [0], [0]
    level_start = [0]
    index = {(): 0}

    for depth in range(1, len(LEVELS)):
        level_start.append(len(label))
        for record in alleles:
            key = tuple(record[:depth])
            if key in index:
                continue
            index[key] = len(label)
            label.append(record[depth - 1])
            parent.append(index[key[:-1]])
            leaf = depth == len(LEVELS) - 1
            novel.append(record[4] if leaf else 0)
            n_genomic.append(record[5] if leaf else 0)
            n_airrseq.append(record[6] if leaf else 0)

    # parent[i] < i, so one backward pass totals every subtree
    for i in range(len(label) - 1, 0, -1):
        p = parent[i]
        novel[p] += novel[i]
        n_genomic[p] += n_genomic[i]
        n_airrseq[p] += n_airrseq[i]

    return {'levels': LEVELS, 'levelStart': level_start, 'label': label,
            'parent': parent, 'novel': novel, 'nG': n_genomic, 'nA': n_airrseq}


@ns.route('/sunburst/<string:species>/<string:locus>')
@api.response(404, 'Species or locus not found')
class SunburstApi(Resource):
    @digby_protected()
    def get(self, species, locus):
        """ The locus hierarchy for the sunburst panel, as depth-ordered parallel arrays """

        error = check_species_locus(species, locus)
        if error:
            return error

        return build(locus, collect(species, locus, requested_sources()))


def _selfcheck():
    """ The invariants the client relies on, on a hand-built locus. """
    rows = [
        ('V', 'IGHV1', 'IGHV1-2', 'IGHV1-2*02', 0, 1, 1),
        ('V', 'IGHV1', 'IGHV1-2', 'IGHV1-2*04', 1, 0, 1),
        ('V', 'IGHV1', 'IGHV1-18', 'IGHV1-18*01', 0, 1, 0),
        ('D', 'IGHD1', 'IGHD1-7', 'IGHD1-7*01', 0, 1, 1),
    ]
    # D before V alphabetically, and *18 before *2 as strings: neither is what we want
    assert sorted(reversed(rows), key=_order) == rows

    out = build('IGH', rows)
    n = len(out['label'])

    assert n == 1 + 2 + 2 + 3 + 4, n                      # chain, segments, subgroups, ascs, alleles
    assert out['parent'][0] == -1
    assert all(out['parent'][i] < i for i in range(1, n))  # the whole point of depth ordering
    assert out['levelStart'] == [0, 1, 3, 5, 8], out['levelStart']
    assert all(out['levelStart'][i] < out['levelStart'][i + 1] for i in range(len(LEVELS) - 1))
    assert out['label'][0] == 'IGH'
    assert out['label'][1:3] == ['V', 'D']                 # V before D, not alphabetical
    assert out['label'][5:8] == ['IGHV1-2', 'IGHV1-18', 'IGHD1-7']   # IGHV1-2 before IGHV1-18

    assert out['novel'][0] == 1 and out['nG'][0] == 3 and out['nA'][0] == 3   # rolled up to the root
    v = out['label'].index('V')
    assert out['nG'][v] == 2 and out['nA'][v] == 2 and out['novel'][v] == 1

    # sorting must not change what the tree contains
    assert sum(1 for p in out['parent'] if p == 0) == 2
    print('ok: %d nodes' % n)


if __name__ == '__main__':
    _selfcheck()
