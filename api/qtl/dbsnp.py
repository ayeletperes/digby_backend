"""dbSNP identifiers for guQTL variants, and the coordinates that reach a browser.

A guQTL variant on IGH is named for where it sits on a locus-relative contig -
`igh_963055` - which is meaningless to dbSNP, Ensembl, gnomAD or UCSC. This
carries the supplied map from that contig to GRCh38 chr14 and to an rsID, so a
variant can be looked up outside VDJbase.

COORDINATE FRAME. Both position columns in the file are 0-based, and nothing in
the file says so. It was established two ways before anything relied on it:

  * Against Ensembl, four of four rsIDs sit one base higher than the file says
    (rs1778 106606198 -> 106606199, and the same for rs1780, rs453878, rs690537).
  * Against `qtl_variant.pos`, which is 1-based, matching on `igh_pos - 1`
    resolves 73.3% of IGH variants where 20.6% is what chance gives at this
    dbSNP density. Matching the columns as they stand resolves 26.5%, which *is*
    chance - so a naive join does not merely miss things, it attaches a wrong
    rsID to most of what it appears to resolve.

Hence `pos_1based = igh_pos + 1` on the way in and `GRCh38_pos + 1` on the way
out, and both are done here so no caller has to remember.

ONE POSITION, SEVERAL rsIDs. 7,107 mapped positions carry more than one, up to
five. They are all returned. Collapsing to one would invent an identity, and
which one you would pick is not a question the data answers.

IGH ONLY. The file is chr14. IGK and IGL variants are already named for their
GRCh38 chromosome and position (`chr2_88872893`), so they need no map to be
linked out by position - but their frame has not been checked the way this one
has, so this module does not claim them.
"""

import gzip
import os

from app import app

# The map lives outside the repository, beside the annotation release it shares a
# reference with. `QTL_DBSNP_PATH` overrides the location.
FILENAME = 'GRCh38_igh_dbSNP_id.tsv.gz'

# The locus-relative contig the map is written against. A variant on any other
# contig is not something this file can speak about.
CONTIG = 'igh'


def _directory():
    return (app.config.get('QTL_DBSNP_PATH')
            or os.environ.get('QTL_DBSNP_PATH')
            or os.path.join(app.config['STATIC_PATH'], 'study_data/QTL/dbsnp'))


_cache = {}


def _load():
    """`{pos_1based: (grch38_contig, grch38_pos_1based, [rsid, ...])}`.

    Read once and kept: a static release of about 245,000 rows. Absent, this
    returns None rather than an empty map - a variant with no rsID and a locus
    with no map are different answers, and only one of them should be shown as
    "no identifier".
    """
    directory = _directory()
    path = os.path.join(directory, FILENAME)
    stamp = os.path.getmtime(path) if os.path.exists(path) else None

    hit = _cache.get(path)
    if hit is not None and hit[0] == stamp:
        return hit[1]

    table = None
    if stamp is not None:
        table = {}
        with gzip.open(path, 'rt') as handle:
            header = handle.readline().rstrip('\n').split('\t')
            try:
                c_chrom = header.index('GRCh38')
                c_gpos = header.index('GRCh38_pos')
                c_ipos = header.index('igh_pos')
                c_rs = header.index('db_SNP_id')
            except ValueError:
                raise ValueError(
                    f'{path} does not have the expected columns '
                    'GRCh38, GRCh38_pos, igh_pos, db_SNP_id')

            for line in handle:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(c_chrom, c_gpos, c_ipos, c_rs):
                    continue
                # 0-based -> 1-based on both columns; see the module docstring
                pos = int(parts[c_ipos]) + 1
                entry = table.get(pos)
                if entry is None:
                    table[pos] = (parts[c_chrom], int(parts[c_gpos]) + 1,
                                  [parts[c_rs]])
                else:
                    entry[2].append(parts[c_rs])

    _cache[path] = (stamp, table)
    return table


def available():
    """Whether the map is configured at all."""
    return _load() is not None


def release():
    """What the map says about itself, for provenance beside the identifiers."""
    try:
        with open(os.path.join(_directory(), 'RELEASE.txt')) as handle:
            return handle.readline().strip() or None
    except OSError:
        return None


def lookup(contig, pos):
    """What dbSNP calls the variant at this 1-based locus position.

    Returns None when there is no map, so a caller can tell "not configured"
    from "configured and this variant is not in dbSNP", which is `mapped: False`.
    """
    table = _load()
    if table is None:
        return None
    if contig != CONTIG or pos is None:
        # said plainly rather than returning an empty result: the light chains
        # are already on chromosome coordinates and this map has nothing to add
        return {'mapped': False, 'applies': False, 'rsids': [],
                'grch38': None, 'source': release()}

    entry = table.get(pos)
    if entry is None:
        return {'mapped': False, 'applies': True, 'rsids': [],
                'grch38': None, 'source': release()}

    chrom, gpos, rsids = entry
    return {
        'mapped': True,
        'applies': True,
        'rsids': rsids,
        # 1-based, the frame dbSNP, Ensembl and UCSC all use
        'grch38': {'contig': chrom, 'pos': gpos},
        'source': release(),
    }
