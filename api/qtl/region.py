"""The genomic neighbourhood a guQTL variant sits in.

An association tells you a variant explains a gene's usage. It does not tell you
what the variant *is* - whether it changes the coding sequence, sits in the
recombination signal, or is 8 kb of intergenic sequence that happens to be in LD
with something that does. That is a question about coordinates, so it is answered
from the same annotation the analysis itself was built against.

COORDINATE FRAME. The BED annotation and the guQTL databases share contigs: `igh`
(locus-relative), `chr2` for IGK, `chr22` for IGL. BED is 0-based half-open;
`qtl_variant.pos` is 1-based. **Everything this module emits is converted to
1-based inclusive**, the frame `pos` is already in, so the caller draws one frame
and never has to know a BED file was involved.

That mapping is not assumed, it is checked: taking the BED interval as
`[start + 1, end]` and testing `start <= pos <= end` reproduces the feature the
database independently assigns for 651 of the 652 variants that fall in one,
across all three loci - the exception being IGHD1-1, a gene these BED files do
not contain at all. Nearest-gene agrees on all 15,689 intergenic variants. The
half-open reading (`start <= pos < end`) instead misses six variants that sit on
an interval boundary, which is what pins the frame.

The annotation lives outside the repository; see `_annotation_dir`.
"""

import os

from flask import request
from flask_restx import Resource
from sqlalchemy import Integer, cast, func

from api.restx import api
from api.system.system import digby_protected
from api.qtl.qtl import _thresholds, check_species_locus, qtl_session
from app import app
from db.qtl_model import Asc, UsageAssociation, Variant

# Its own namespace so it can be registered separately, but mounted on the guQTL
# path: this is one more question about a guQTL variant, so it belongs beside the
# endpoints that answer the others rather than on a prefix of its own.
ns = api.namespace('qtl_region', path='/qtl',
                   description='The genomic neighbourhood of a guQTL variant')

DEFAULT_WINDOW = 20000
MAX_WINDOW = 2000000
# below this a window holds no context at all, and a stray drag can produce it
MIN_WINDOW = 200

# BED file -> the name the guQTL database gives that feature, so the track and
# the variant panel use one vocabulary rather than two. `exon_2` is deliberately
# absent: it is `l-part2` followed by `region`, both of which are drawn, so
# including it would lay a second rectangle over the top of those two.
BED_KINDS = {
    'gene': ('gene', 'gene'),
    'utr': ('utr', 'utr'),
    'exon_1': ('leader', 'l_part1'),
    'intron': ('leader', 'leader_intron'),
    'l-part2': ('leader', 'l_part2'),
    'region': ('coding', 'region'),
    'heptamer': ('rss', 'heptamer'),
    'spacer': ('rss', 'spacer'),
    'nonamer': ('rss', 'nonamer'),
    'constant': ('constant', 'constant'),
}


def _annotation_dir():
    """Where the BED annotation lives.

    Configured, never hardcoded: these files are a versioned release of an
    external annotation, not repository content. `QTL_ANNOTATION_PATH` in
    `secret.cfg` (or the environment) wins; otherwise the deployment is expected
    to place the release beside the databases it annotates.
    """
    return (app.config.get('QTL_ANNOTATION_PATH')
            or os.environ.get('QTL_ANNOTATION_PATH')
            or os.path.join(app.config['STATIC_PATH'], 'study_data/QTL/annotation'))


_cache = {}


def _bed(name):
    """One BED file as `(contig, start, end, gene)`, in 1-based inclusive coords.

    Static release files, so they are read once and kept. Each is a few hundred
    rows; the whole set is under 100 kB.
    """
    directory = _annotation_dir()
    if not os.path.isdir(directory):
        # returning nothing here would report every variant as intergenic, which
        # is a wrong answer rather than a missing one
        raise FileNotFoundError(
            f'guQTL annotation not found at {directory}. This is the BED release the '
            'analysis was run against; set QTL_ANNOTATION_PATH to point at it.')

    path = os.path.join(directory, name + '.bed')
    stamp = os.path.getmtime(path) if os.path.exists(path) else None

    hit = _cache.get(path)
    if hit is not None and hit[0] == stamp:
        return hit[1]

    rows = []
    if stamp is not None:
        with open(path) as handle:
            for line in handle:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 4 or parts[0].startswith(('#', 'track', 'browser')):
                    continue
                # 0-based half-open -> 1-based inclusive; see the module docstring
                rows.append((parts[0], int(parts[1]) + 1, int(parts[2]), parts[3]))

    _cache[path] = (stamp, rows)
    return rows


def _release(directory):
    """The annotation release, as the release itself states it.

    The directory name is not the release: this one is called `annotation`, while
    the files in it are `immune_receptor_genomics release 251106`. A different
    release moves every feature call, so what travels with the answer has to be
    the version, not where it happens to be mounted.
    """
    try:
        with open(os.path.join(directory, 'RELEASE.txt')) as handle:
            return handle.readline().strip() or None
    except OSError:
        return None


def _overlapping(name, contig, start, end):
    return [row for row in _bed(name)
            if row[0] == contig and row[1] <= end and row[2] >= start]


@ns.route('/region/<string:species>/<string:locus>/<path:variant>')
@api.response(404, 'Species, locus or variant not found')
class QtlRegionApi(Resource):
    @digby_protected()
    def get(self, species, locus, variant):
        """ Returns the genes, functional elements and tested variants around a variant

        The window is centred on the variant, `window` bp wide (default 20 kb).
        With `asc`, the neighbouring variants are plotted against that ASC's scan;
        without, against their strongest signal across any ASC - the same choice
        the Manhattan view makes, so a click through keeps its meaning.

        Coordinates are 1-based inclusive throughout, matching `pos`.
        """

        error = check_species_locus(species, locus)
        if error:
            return error

        session = qtl_session(species, locus)
        record = session.query(Variant).filter(Variant.variant == variant).one_or_none()
        if record is None:
            return {'message': f'No such variant: {variant}'}, 404
        if record.pos is None or record.contig is None:
            return {'message': f'{variant} has no position, so it has no neighbourhood'}, 404

        # Two ways to say which stretch to draw. `window` centres it on the
        # variant, which is what opening the panel wants; an explicit `start`/`end`
        # is what dragging a range on the track wants, and that range does not
        # generally have the variant in the middle of it - or in it at all.
        try:
            window = int(request.args.get('window') or DEFAULT_WINDOW)
            explicit = (int(request.args['start']), int(request.args['end'])) \
                if request.args.get('start') and request.args.get('end') else None
        except ValueError:
            return {'message': 'window, start and end must be whole numbers'}, 400

        if explicit is not None:
            start, end = explicit
            if end <= start:
                return {'message': 'end must be greater than start'}, 400
            start = max(1, start)
            # clamped, not rejected: a drag that runs off the end of the contig is
            # an ordinary gesture, and the same ceiling applies either way
            end = min(end, start + MAX_WINDOW)
            if end - start < MIN_WINDOW:
                end = start + MIN_WINDOW
        else:
            window = max(MIN_WINDOW, min(window, MAX_WINDOW))
            start = max(1, record.pos - window // 2)
            end = start + window

        asc = request.args.get('asc')

        directory = _annotation_dir()
        annotated = os.path.isdir(directory)

        genes, features = [], []
        if annotated:
            for name, (feature, kind) in BED_KINDS.items():
                for _, first, last, gene in _overlapping(name, record.contig, start, end):
                    row = {'name': gene, 'start': first, 'end': last,
                           'feature': feature, 'kind': kind}
                    (genes if kind == 'gene' else features).append(row)

        genes.sort(key=lambda g: g['start'])
        features.sort(key=lambda f: (f['start'], f['end']))

        return {
            'species': species,
            'locus': locus,
            'contig': record.contig,
            'window': {'start': start, 'end': end, 'size': end - start,
                       'centre': record.pos},
            'asc': asc,
            'variant': {'variant': record.variant, 'contig': record.contig,
                        'pos': record.pos, 'maf': record.maf, 'gene': record.gene,
                        'feature': record.feature, 'sub_feature': record.sub_feature,
                        'distance_to_gene': record.distance_to_gene},
            'genes': genes,
            'features': features,
            'variants': _neighbours(session, record, start, end, asc),
            'thresholds': _thresholds(session),
            # said rather than implied: an unconfigured annotation directory draws
            # an empty track, and an empty track must not read as empty sequence
            'annotation': {'available': annotated,
                           'source': _release(directory) if annotated else None},
        }


def _neighbours(session, record, start, end, asc=None):
    """Every tested variant in the window, with how strong its signal is.

    One row per variant, not per variant-and-ASC: this is a track, so a variant is
    one mark. Without an ASC that is the variant's strongest signal against any of
    them, which is what the Manhattan view shows too.
    """
    query = (
        session.query(Variant.variant, Variant.pos, Variant.maf, Variant.gene,
                      Variant.feature, Variant.sub_feature,
                      func.max(UsageAssociation.neglog10_p),
                      cast(UsageAssociation.significant, Integer),
                      Asc.asc)
        .join(UsageAssociation, UsageAssociation.variant_id == Variant.id)
        .join(Asc, Asc.id == UsageAssociation.asc_id)
        .filter(Variant.contig == record.contig,
                Variant.pos >= start, Variant.pos <= end)
        .group_by(Variant.id)
    )
    if asc:
        query = query.filter(Asc.asc == asc)

    rows = [{'variant': name, 'pos': pos, 'maf': maf, 'gene': gene,
             'feature': feature, 'sub_feature': sub_feature,
             'neglog10_p': neglog10_p, 'significant': bool(significant),
             'asc': best_asc, 'selected': name == record.variant}
            for name, pos, maf, gene, feature, sub_feature, neglog10_p,
                significant, best_asc in query.all()]
    rows.sort(key=lambda v: v['pos'])
    return rows


def _self_check():
    """The coordinate frame, checked against facts that fail if it slips by one.

    From the backend directory, with QTL_ANNOTATION_PATH set:

        python -c "import app; from api.qtl.region import _self_check; _self_check()"

    (`app` first: it is the package's entry point, and importing this module
    ahead of it walks into the app/restx circular import.)
    """
    # 1. The elements of a V gene tile its body exactly: UTR, L-PART1, intron,
    #    L-PART2, V-REGION, heptamer, spacer, nonamer, end to end with no gap and
    #    no overlap. Reading the BED as half-open opens a 1 bp gap at every join.
    for gene, contig in (('IGLV3-16', 'chr22'), ('IGKV1-17', 'chr2')):
        parts = sorted((first, last, kind)
                       for name, (_, kind) in BED_KINDS.items() if kind != 'gene'
                       for c, first, last, g in _bed(name)
                       if g == gene and c == contig)
        body = [row for row in _bed('gene') if row[3] == gene and row[0] == contig]
        assert len(body) == 1, gene
        assert parts[0][0] == body[0][1] and parts[-1][1] == body[0][2], (gene, parts)
        for (_, last, kind), (first, _, nxt) in zip(parts, parts[1:]):
            assert first == last + 1, f'{gene}: {kind} -> {nxt} is not contiguous'

    # 2. Two variants whose annotation is independently published. The spacer one
    #    sits on the interval's last base under the correct frame and outside it
    #    under the wrong one, which is what makes it a test rather than a sample.
    def feature_at(contig, pos, gene):
        return {kind for name, (_, kind) in BED_KINDS.items()
                for c, first, last, g in _bed(name)
                if c == contig and g == gene and first <= pos <= last}

    assert 'spacer' in feature_at('chr22', 23170898, 'IGLV3-16')
    assert 'gene' in feature_at('chr22', 23170898, 'IGLV3-16')
    assert not feature_at('chr22', 22756855, 'IGLV9-49')      # 48 bp upstream, outside
    body = [r for r in _bed('gene') if r[3] == 'IGLV9-49'][0]
    assert body[1] - 22756855 == 49, body    # the database says 48: see the docstring

    print('coordinate frame OK')
