"""Build a guQTL database from an igqtl.R run directory.

One database per species and locus. The run tree holds every locus together, so
the loading is filtered as it streams rather than read whole: the largest table
is a million rows and there is no reason to hold it in memory.

Layout expected (see husa_manuscript/docs/analyses/igqtl.md):

    <run>/provenance/manifest.json
    <run>/source_data/usage_associations_<LOCUS>.tsv.gz
    <run>/source_data/{asc_usage,dosage_long}.tsv.gz
    <run>/source_data/{pairing_associations,cell_tests,dj_enrichment}.tsv.gz
    <run>/reports/{usage_thresholds,thresholds,variant_features,usage_leads}.tsv
"""

import csv
import gzip
import json
import math
import os
import sqlite3

from sqlalchemy import create_engine

from db.qtl_model import Base

# rows per executemany; large enough to amortise the call, small enough that a
# million-row table never lands in memory at once
BATCH = 20000

# p-values are reported down to 1e-53 but a hard zero would break the log
MIN_P = 1e-300


def _rows(path, gz=None):
    """Stream a TSV, gzipped or not, as dicts."""
    if gz is None:
        gz = path.endswith('.gz')
    opener = gzip.open if gz else open
    with opener(path, 'rt', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def _num(value, cast=float):
    """A number, or None for the empty strings the run tree uses for absent."""
    if value is None or value == '' or value == 'NA':
        return None
    try:
        return cast(value)
    except ValueError:
        return None


def _bool(value):
    if value is None or value == '':
        return None
    return 1 if str(value).upper() in ('TRUE', 'T', '1', 'YES') else 0


def _neglog10(p):
    return -math.log10(max(p, MIN_P)) if p is not None else None


def _genotype(dosage):
    """The dosage as a called genotype: 0, 1 or 2, which is what the boxplot groups on."""
    return min(2, max(0, int(round(dosage))))


def read_manifest(run_dir):
    path = os.path.join(run_dir, 'provenance', 'manifest.json')
    if not os.path.exists(path):
        return {}
    with open(path) as handle:
        return json.load(handle)


def genotype_matrix(run_dir):
    """The cohort genotype matrix the run was given, or None.

    The manifest records it as `config.genotype`, but relative to the manuscript
    root, which sits above the run tree and is written down nowhere. So the path
    is resolved by walking up from the run directory until it exists - the run is
    a few levels below the root, and the relative path is specific enough that
    nothing else answers to it.
    """
    genotype = (read_manifest(run_dir).get('config') or {}).get('genotype')
    if not genotype:
        return None
    if os.path.isabs(genotype):
        return genotype if os.path.exists(genotype) else None

    base = os.path.realpath(run_dir)
    while True:
        candidate = os.path.join(base, genotype)
        if os.path.exists(candidate):
            return candidate
        parent = os.path.dirname(base)
        if parent == base:
            return None
        base = parent


class QtlBuilder:
    """Loads one locus of a run into a fresh database."""

    def __init__(self, run_dir, locus, path):
        self.run_dir = run_dir
        self.locus = locus
        self.path = path

        if os.path.exists(path):
            os.remove(path)
        os.makedirs(os.path.dirname(path), exist_ok=True)

        # the schema comes from the models, the loading from raw sqlite: the ORM
        # would spend most of the run building objects nobody reads
        Base.metadata.create_all(create_engine('sqlite:///' + path))
        self.con = sqlite3.connect(path)
        self.con.execute('PRAGMA journal_mode = OFF')
        self.con.execute('PRAGMA synchronous = OFF')

        self.variants = {}      # variant name -> id
        self.ascs = {}          # asc name -> id
        self.subjects = {}      # subject -> id

    # ------------------------------------------------------------------ util

    def _insert(self, table, columns, rows):
        if not rows:
            return 0
        sql = (f"INSERT OR IGNORE INTO {table} ({','.join(columns)}) "
               f"VALUES ({','.join('?' * len(columns))})")
        # rows written, not rows offered: OR IGNORE drops duplicates, and a count
        # that ignores that misreports what is actually in the database
        return self.con.executemany(sql, rows).rowcount

    def _source(self, name):
        return os.path.join(self.run_dir, 'source_data', name)

    def _report(self, name):
        return os.path.join(self.run_dir, 'reports', name)

    # --------------------------------------------------------------- loading

    def load_run(self):
        manifest = read_manifest(self.run_dir)

        self.con.execute(
            "INSERT INTO qtl_run (label, generated_at, script, config) VALUES (?,?,?,?)",
            (os.path.basename(os.path.realpath(self.run_dir)),
             manifest.get('generated_at'),
             manifest.get('script'),
             json.dumps(manifest.get('config', {}))))

    def load_thresholds(self):
        rows = []
        for name in ('usage_thresholds.tsv', 'thresholds.tsv'):
            path = self._report(name)
            if not os.path.exists(path):
                continue
            for row in _rows(path):
                if row.get('locus') != self.locus:
                    continue
                # empty string rather than NULL: the usage row appears in both
                # files, and SQLite treats NULLs as distinct, so a UNIQUE over a
                # nullable conditional would not have deduplicated it
                rows.append((
                    row.get('analysis'), row.get('grouped_by') or '',
                    row.get('conditional') or '',
                    _num(row.get('n_subjects'), int), _num(row.get('n_excluded'), int),
                    _num(row.get('n_variants'), int), _num(row.get('n_independent'), int),
                    _num(row.get('n_asc'), int), _num(row.get('threshold')),
                    _num(row.get('n_significant_variants'), int),
                    _num(row.get('n_independent_significant'), int)))

        return self._insert('qtl_threshold',
                            ['analysis', 'grouped_by', 'conditional', 'n_subjects',
                             'n_excluded', 'n_variants', 'n_independent', 'n_asc',
                             'threshold', 'n_significant_variants',
                             'n_independent_significant'], rows)

    def load_usage_associations(self):
        """Variants, ASCs and the association rows, in one pass over the big file.

        The association file is the authority on which variants and ASCs this
        locus tested, so the reference tables are built from it rather than from
        the annotation files, which cover only some of them.

        Note that it holds fewer variants than the run reports as tested: 9,402
        against the 9,792 in qtl_threshold.n_variants for IGH. The remainder were
        genotyped - they are in dosage_long - but produced no usage association,
        which the pipeline docs attribute to the complete-case subject set applied
        inside the scan. Both numbers are kept so a plot showing 9,402 points can
        say what it is 9,402 of.
        """
        path = self._source(f'usage_associations_{self.locus}.tsv.gz')
        if not os.path.exists(path):
            return 0

        associations = []
        written = 0

        for row in _rows(path):
            variant = row['variant']
            if variant not in self.variants:
                self.variants[variant] = len(self.variants) + 1
                self.con.execute(
                    "INSERT INTO qtl_variant (id, variant, contig, pos, maf) VALUES (?,?,?,?,?)",
                    (self.variants[variant], variant,
                     # split on the FIRST underscore, not the last: ids are
                     # contig_pos, but a multi-allelic site adds a third part
                     # (chr22_23161341_2), and rsplit handed that whole
                     # `chr22_23161341` back as the contig. No contig carries an
                     # underscore, so the first one always ends it.
                     variant.split('_', 1)[0] if '_' in variant else None,
                     _num(row.get('pos'), int), _num(row.get('maf'))))

            asc = row['asc']
            if asc not in self.ascs:
                self.ascs[asc] = len(self.ascs) + 1
                self.con.execute(
                    "INSERT INTO qtl_asc (id, asc, segment, asc_position, asc_span, n_member) "
                    "VALUES (?,?,?,?,?,?)",
                    (self.ascs[asc], asc, row.get('segment'),
                     _num(row.get('asc_position')), _num(row.get('asc_span')),
                     _num(row.get('n_member'), int)))

            p = _num(row.get('p_value'))
            associations.append((
                self.variants[variant], self.ascs[asc], _num(row.get('n'), int),
                _num(row.get('beta')), _num(row.get('se')), _num(row.get('t_stat')),
                p, _neglog10(p), _bool(row.get('significant')),
                _num(row.get('distance_to_asc'))))

            if len(associations) >= BATCH:
                written += self._insert('qtl_usage_association', _ASSOC_COLS, associations)
                associations = []

        written += self._insert('qtl_usage_association', _ASSOC_COLS, associations)
        return written

    def annotate_variants(self):
        """Add the genomic feature each variant falls in, where it is known."""
        path = self._report('variant_features.tsv')
        if not os.path.exists(path):
            return 0

        updates = []
        for row in _rows(path):
            if row.get('locus') != self.locus:
                continue
            variant_id = self.variants.get(row['variant'])
            if variant_id:
                updates.append((row.get('gene'), row.get('feature'),
                                row.get('sub_feature'), _num(row.get('distance_to_gene')),
                                variant_id))

        self.con.executemany(
            "UPDATE qtl_variant SET gene=?, feature=?, sub_feature=?, distance_to_gene=? "
            "WHERE id=?", updates)
        return len(updates)

    def mark_leads(self):
        """Flag the lead variants and carry the power columns onto their rows.

        ld_group is skipped: it is the whole standardised dosage vector as a
        comma-string, about 2KB per row, and nothing reads it here.
        """
        path = self._report('usage_leads.tsv')
        if not os.path.exists(path):
            return 0

        updates = []
        for row in _rows(path):
            if row.get('locus') != self.locus:
                continue
            variant_id = self.variants.get(row['variant'])
            asc_id = self.ascs.get(row['asc'])
            if variant_id and asc_id:
                updates.append((
                    _num(row.get('min_genotype_group'), int), _bool(row.get('well_powered')),
                    _bool(row.get('is_cis')), variant_id, asc_id))

        self.con.executemany(
            "UPDATE qtl_usage_association SET min_genotype_group=?, well_powered=?, "
            "is_cis=?, is_lead=1 WHERE variant_id=? AND asc_id=?", updates)
        return len(updates)

    def load_asc_usage(self):
        path = self._source('asc_usage.tsv.gz')
        if not os.path.exists(path):
            return 0

        rows = []
        written = 0
        for row in _rows(path):
            if row.get('locus') != self.locus:
                continue

            asc_id = self.ascs.get(row['asc'])
            if not asc_id:
                continue            # an ASC with a phenotype but no tested variant

            rows.append((self._register_subject(row['subject']), asc_id,
                         _num(row.get('count'), int),
                         _num(row.get('total'), int), _num(row.get('n_asc'), int),
                         _num(row.get('usage')), _num(row.get('logit_usage'))))

            if len(rows) >= BATCH:
                written += self._insert('qtl_asc_usage', _USAGE_COLS, rows)
                rows = []

        return written + self._insert('qtl_asc_usage', _USAGE_COLS, rows)

    def _register_subject(self, subject):
        if subject not in self.subjects:
            self.subjects[subject] = len(self.subjects) + 1
            self.con.execute("INSERT INTO qtl_subject (id, subject) VALUES (?,?)",
                             (self.subjects[subject], subject))
        return self.subjects[subject]

    def load_dosage(self, genotypes=None):
        """Genotypes, from the cohort matrix where there is one.

        Two sources say the same thing in different shapes. dosage_long.tsv.gz
        is written by the run, but only for IGH, so a database built from it
        leaves IGK and IGL with no genotypes at all and no per-variant plot. The
        matrix the run was *given* (`config.genotype`: one row per variant, one
        column per subject) covers every locus, so it is preferred, and
        dosage_long is the fallback for a run whose matrix cannot be found.

        Either way the rows are filtered to this locus by variant membership:
        neither source carries a locus column, so the variants tested here are
        the only way to tell which rows belong.
        """
        if genotypes is None:
            genotypes = genotype_matrix(self.run_dir)
        if genotypes:
            return self._load_dosage_matrix(genotypes)
        return self._load_dosage_long()

    def _load_dosage_matrix(self, path):
        """Read the wide genotype matrix, streaming one variant per row.

        Filtered on both axes as it goes rather than read whole: it holds the
        entire cohort at every locus, and this database wants one locus of it.
        """
        rows = []
        written = 0

        with open(path, newline='') as handle:
            reader = csv.reader(handle, delimiter='\t')
            header = next(reader)

            # load_asc_usage has already registered the subjects this locus
            # phenotyped, and the matrix carries the rest of the cohort besides;
            # a genotype with no usage to plot it against is not worth a row. If
            # no phenotypes were loaded at all, keep the whole cohort rather than
            # silently writing an empty table.
            if not self.subjects:
                for subject in header[1:]:
                    self._register_subject(subject)
            columns = [(index, self.subjects[subject])
                       for index, subject in enumerate(header)
                       if index and subject in self.subjects]

            for record in reader:
                variant_id = self.variants.get(record[0])
                if not variant_id:
                    continue

                for index, subject_id in columns:
                    dosage = _num(record[index])
                    if dosage is None:
                        continue        # NA: no call for this subject
                    rows.append((variant_id, subject_id, dosage, _genotype(dosage)))

                if len(rows) >= BATCH:
                    written += self._insert('qtl_dosage', _DOSAGE_COLS, rows)
                    rows = []

        return written + self._insert('qtl_dosage', _DOSAGE_COLS, rows)

    def _load_dosage_long(self):
        path = self._source('dosage_long.tsv.gz')
        if not os.path.exists(path):
            return 0

        rows = []
        written = 0
        for row in _rows(path):
            variant_id = self.variants.get(row['variant'])
            if not variant_id:
                continue

            rows.append((variant_id, self._register_subject(row['subject']),
                         _num(row.get('dosage')), _num(row.get('genotype'), int)))

            if len(rows) >= BATCH:
                written += self._insert('qtl_dosage', _DOSAGE_COLS, rows)
                rows = []

        return written + self._insert('qtl_dosage', _DOSAGE_COLS, rows)

    def load_pairing(self):
        """The D/J pairing scans. Only IGH has them, and only for its variants."""
        totals = {}

        for name, conditional_default, anchor_col in (
                ('pairing_associations.tsv.gz', 'P(J|D)', 'j_gene'),
                ('pairing_associations_by_d.tsv.gz', 'P(D|J)', 'd_gene')):
            path = self._source(name)
            if not os.path.exists(path):
                continue

            rows = []
            for row in _rows(path):
                variant_id = self.variants.get(row['variant'])
                if not variant_id:
                    continue
                p = _num(row.get('p_value'))
                rows.append((row.get('conditional') or conditional_default, variant_id,
                             row.get(anchor_col), _num(row.get('n'), int),
                             _num(row.get('pillai')), _num(row.get('f_stat')),
                             p, _neglog10(p), _num(row.get('min_genotype_group'), int),
                             _bool(row.get('significant'))))
            totals[name] = self._insert('qtl_pairing_association', _PAIRING_COLS, rows)

        path = self._source('cell_tests.tsv.gz')
        if os.path.exists(path):
            rows = []
            for row in _rows(path):
                variant_id = self.variants.get(row['variant'])
                if not variant_id:
                    continue
                rows.append((row.get('conditional'), variant_id, row.get('d_gene'),
                             row.get('j_gene'), _num(row.get('n'), int),
                             _num(row.get('beta')), _num(row.get('p_value')),
                             _num(row.get('mean_low')), _num(row.get('mean_high')),
                             _num(row.get('delta_mean')), _num(row.get('n_low'), int),
                             _num(row.get('n_high'), int), _num(row.get('omnibus_p_value')),
                             _bool(row.get('omnibus_significant')),
                             _num(row.get('min_genotype_group'), int),
                             _bool(row.get('marked')), _bool(row.get('marked_strict'))))
            totals['cell_tests'] = self._insert('qtl_cell_test', _CELL_COLS, rows)

        # dj_enrichment carries no variant and no locus, only a subject - and the
        # subjects are shared across loci, so matching on subject alone loaded the
        # IGH pairing data into every locus. It belongs to whichever locus actually
        # ran the pairing scan, which is the one that loaded pairing associations.
        ran_pairing = any(totals.get(name) for name in
                          ('pairing_associations.tsv.gz', 'pairing_associations_by_d.tsv.gz'))

        path = self._source('dj_enrichment.tsv.gz')
        if ran_pairing and os.path.exists(path) and self.subjects:
            rows = []
            for row in _rows(path):
                subject_id = self.subjects.get(row['subject'])
                if not subject_id:
                    continue
                rows.append((subject_id, row.get('d_gene'), row.get('j_gene'),
                             _num(row.get('count'), int), _num(row.get('depth'), int),
                             _num(row.get('d_total'), int), _num(row.get('j_total'), int),
                             _num(row.get('expected')), _num(row.get('enrichment')),
                             _num(row.get('p_d')), _num(row.get('p_j')),
                             _num(row.get('p_j_given_d')), _num(row.get('p_d_given_j'))))
            totals['dj_enrichment'] = self._insert('qtl_dj_enrichment', _DJ_COLS, rows)

        return totals

    def finish(self, species):
        self.con.execute(
            "INSERT INTO details (dbtype, species, locus, created_on, created_by) "
            "VALUES (?,?,?,datetime('now'),?)",
            ('guQTL', species, self.locus, 'make_qtl_db'))
        self.con.commit()
        self.con.execute('PRAGMA optimize')
        self.con.close()

        with open(os.path.join(os.path.dirname(self.path), 'db_description.txt'), 'w') as fo:
            fo.write(f'Gene-usage QTL results for {species} {self.locus}')


_ASSOC_COLS = ['variant_id', 'asc_id', 'n', 'beta', 'se', 't_stat', 'p_value',
               'neglog10_p', 'significant', 'distance_to_asc']
_USAGE_COLS = ['subject_id', 'asc_id', 'count', 'total', 'n_asc', 'usage', 'logit_usage']
_DOSAGE_COLS = ['variant_id', 'subject_id', 'dosage', 'genotype']
_PAIRING_COLS = ['conditional', 'variant_id', 'anchor_gene', 'n', 'pillai', 'f_stat',
                 'p_value', 'neglog10_p', 'min_genotype_group', 'significant']
_CELL_COLS = ['conditional', 'variant_id', 'd_gene', 'j_gene', 'n', 'beta', 'p_value',
              'mean_low', 'mean_high', 'delta_mean', 'n_low', 'n_high',
              'omnibus_p_value', 'omnibus_significant', 'min_genotype_group',
              'marked', 'marked_strict']
_DJ_COLS = ['subject_id', 'd_gene', 'j_gene', 'count', 'depth', 'd_total', 'j_total',
            'expected', 'enrichment', 'p_d', 'p_j', 'p_j_given_d', 'p_d_given_j']


def loci_in(run_dir):
    """The loci this run produced usage associations for."""
    source = os.path.join(run_dir, 'source_data')
    if not os.path.isdir(source):
        return []
    return sorted(name[len('usage_associations_'):-len('.tsv.gz')]
                  for name in os.listdir(source)
                  if name.startswith('usage_associations_') and name.endswith('.tsv.gz'))


def build(run_dir, species, locus, static_path, genotypes=None):
    """Build one locus, returning the row counts written.

    `genotypes` overrides the cohort genotype matrix; with none given the run's
    manifest is asked where its own was (see genotype_matrix).
    """
    path = os.path.join(static_path, 'study_data', 'QTL', 'db', species, locus, 'db.sqlite3')
    builder = QtlBuilder(run_dir, locus, path)

    counts = {}
    builder.load_run()
    counts['thresholds'] = builder.load_thresholds()
    counts['usage_associations'] = builder.load_usage_associations()
    counts['variants'] = len(builder.variants)
    counts['ascs'] = len(builder.ascs)
    counts['annotated'] = builder.annotate_variants()
    counts['leads'] = builder.mark_leads()
    counts['asc_usage'] = builder.load_asc_usage()
    counts['dosage'] = builder.load_dosage(genotypes)
    counts['subjects'] = len(builder.subjects)
    counts.update(builder.load_pairing())
    builder.finish(species)

    counts['path'] = path
    return counts
