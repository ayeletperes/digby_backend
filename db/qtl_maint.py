"""Build a guQTL database from an igqtl.R run directory.

One database per species and locus. The run tree holds every locus together, so
the loading is filtered as it streams rather than read whole: the largest table
is a million rows and there is no reason to hold it in memory.

Layout expected: one directory per locus, holding the files the run wrote for it

    manifest.json                       schema, thresholds, provenance
    usage_associations_<LOCUS>.tsv.gz   one row per variant and cluster
    variant_features.tsv.gz             everything true of a variant
    asc_features.tsv.gz                 everything true of a cluster
    genotypes.matrix.gz                 variant x subject, this locus
    asc_usage.tsv.gz                    the usage phenotype
    pairing.tsv.gz, cell_tests.tsv.gz, dj_enrichment.tsv.gz    IGH only
"""

import csv
import gzip
import json
import math
import os
import sqlite3
import sys

from sqlalchemy import create_engine

from db.qtl_model import Base

# rows per executemany; large enough to amortise the call, small enough that a
# million-row table never lands in memory at once
BATCH = 20000

# p-values are reported down to 1e-53 but a hard zero would break the log
MIN_P = 1e-300


# what _rows actually read this build: path -> {'columns': asked, 'n_rows': streamed}.
# Derived rather than declared, so it cannot drift out of step with the loaders.
_READ = {}

# the one column a file may legitimately not carry: only the pairing scans are
# conditional, and the usage tables share their loader
_OPTIONAL = {'conditional'}


class _Row(dict):
    """A row that remembers which columns were asked of it."""

    def __init__(self, row, asked):
        super().__init__(row)
        self.asked = asked

    def __getitem__(self, key):
        self.asked.add(key)
        return super().__getitem__(key)

    def get(self, key, default=None):
        self.asked.add(key)
        return super().get(key, default)


def _rows(path, gz=None):
    """Stream a TSV, gzipped or not, as dicts, recording what was read."""
    if gz is None:
        gz = path.endswith('.gz')
    opener = gzip.open if gz else open
    seen = _READ.setdefault(path, {'columns': set(), 'n_rows': 0})
    count = 0
    with opener(path, 'rt', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            count += 1
            yield _Row(row, seen['columns'])
    seen['n_rows'] = count


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


def read_manifest(run_dir, manifest=None):
    path = manifest or os.path.join(run_dir, 'manifest.json')
    if not os.path.exists(path):
        return {}
    with open(path) as handle:
        return json.load(handle)


def check_schema(run_dir, manifest=None):
    """Stop if the run wrote something other than what the loaders read.

    The manifest declares the columns and row count of every file the run
    produced. Every column asked for has to be one the run says it wrote, and
    every row has to have been read. A column that renames upstream otherwise
    arrives as a table full of NULLs: a build that succeeds and a database that
    is wrong.
    """
    declared = read_manifest(run_dir, manifest).get('schema')
    if not declared:
        print('Manifest declares no schema: this run predates the contract, columns unchecked')
        return

    problems = []
    unchecked = []
    for path, seen in sorted(_READ.items()):
        name = os.path.relpath(path, run_dir).replace(os.sep, '/')
        spec = declared.get(name)
        if spec is None:
            unchecked.append(name)
            continue

        missing = sorted(seen['columns'] - set(spec.get('columns') or ()) - _OPTIONAL)
        if missing:
            problems.append(f'{name}: read {", ".join(missing)}, which the run did not write')
        if spec.get('n_rows') is not None and spec['n_rows'] != seen['n_rows']:
            problems.append(f"{name}: read {seen['n_rows']} rows, the run wrote {spec['n_rows']}")

    if unchecked:
        print('Not declared in the manifest, so unchecked: ' + ', '.join(unchecked))
    if problems:
        sys.exit('The run does not match what this builder reads:\n  ' + '\n  '.join(problems))


def genotype_matrix(run_dir):
    """The cohort genotype matrix, which the run now carries beside its tables."""
    path = os.path.join(run_dir, 'genotypes.matrix.gz')
    return path if os.path.exists(path) else None


class QtlBuilder:
    """Loads one locus of a run into a fresh database."""

    def __init__(self, run_dir, locus, path, manifest=None):
        self.run_dir = run_dir
        self.locus = locus
        self.path = path
        self.manifest = manifest

        # a bare filename has no directory part, and makedirs('') raises
        directory = os.path.dirname(path)
        if directory:
            os.makedirs(directory, exist_ok=True)

        # the schema comes from the models, the loading from raw sqlite: the ORM
        # would spend most of the run building objects nobody reads. create_all
        # leaves an existing database alone, which is what lets a second study
        # be added to the locus without rebuilding the first.
        Base.metadata.create_all(create_engine('sqlite:///' + path))
        self.con = sqlite3.connect(path)
        self.con.execute('PRAGMA journal_mode = OFF')
        self.con.execute('PRAGMA synchronous = OFF')

        self.run_id = None      # set by load_run; scopes everything this build writes
        self.variants = {}      # variant name -> id, within this run
        self.ascs = {}          # asc name -> id, within this run
        self.subjects = {}      # subject -> id, within this run

    # ------------------------------------------------------------------ util

    def _insert(self, table, columns, rows):
        if not rows:
            return 0
        sql = (f"INSERT OR IGNORE INTO {table} ({','.join(columns)}) "
               f"VALUES ({','.join('?' * len(columns))})")
        # rows written, not rows offered: OR IGNORE drops duplicates, and a count
        # that ignores that misreports what is actually in the database
        return self.con.executemany(sql, rows).rowcount

    def _file(self, name):
        return os.path.join(self.run_dir, name)

    # --------------------------------------------------------------- loading

    def load_run(self, project=None):
        """Record this study's run, replacing any earlier build of it.

        Rebuilding one study must not disturb the others, so the old run's rows
        go and the new ones take their place. Everything hangs off the four
        tables that carry run_id, so deleting those and their dependants is the
        whole of it.
        """
        manifest = read_manifest(self.run_dir, self.manifest)

        previous = self.con.execute(
            'SELECT id FROM qtl_run WHERE project IS ?', (project,)).fetchall()
        for (old_id,) in previous:
            self._forget(old_id)

        # everything the manifest states about where the run came from, kept as
        # one blob: provenance travels with the database, not beside it
        provenance = {k: manifest[k] for k in
                      ('source_run', 'annotation_release', 'locus', 'files')
                      if k in manifest}

        cur = self.con.execute(
            "INSERT INTO qtl_run (label, generated_at, script, project, config) "
            "VALUES (?,?,?,?,?)",
            (project or os.path.basename(os.path.realpath(self.run_dir)),
             manifest.get('generated_at'),
             manifest.get('source_run'),
             project,
             json.dumps(provenance)))
        self.run_id = cur.lastrowid

    def _forget(self, run_id):
        """Remove one study's rows, children first."""
        for sql in (
                "DELETE FROM qtl_usage_association WHERE variant_id IN "
                "  (SELECT id FROM qtl_variant WHERE run_id = ?)",
                "DELETE FROM qtl_dosage WHERE variant_id IN "
                "  (SELECT id FROM qtl_variant WHERE run_id = ?)",
                "DELETE FROM qtl_pairing_association WHERE variant_id IN "
                "  (SELECT id FROM qtl_variant WHERE run_id = ?)",
                "DELETE FROM qtl_cell_test WHERE variant_id IN "
                "  (SELECT id FROM qtl_variant WHERE run_id = ?)",
                "DELETE FROM qtl_asc_usage WHERE subject_id IN "
                "  (SELECT id FROM qtl_subject WHERE run_id = ?)",
                "DELETE FROM qtl_dj_enrichment WHERE subject_id IN "
                "  (SELECT id FROM qtl_subject WHERE run_id = ?)",
                "DELETE FROM qtl_variant WHERE run_id = ?",
                "DELETE FROM qtl_asc WHERE run_id = ?",
                "DELETE FROM qtl_subject WHERE run_id = ?",
                "DELETE FROM qtl_threshold WHERE run_id = ?",
                "DELETE FROM qtl_run WHERE id = ?"):
            self.con.execute(sql, (run_id,))

    def load_thresholds(self):
        """The thresholds the run applied, which it states in its manifest."""
        rows = []
        for row in read_manifest(self.run_dir, self.manifest).get('thresholds') or []:
            if row.get('locus') != self.locus:
                continue
            # empty string rather than NULL: SQLite treats NULLs as distinct, so
            # a UNIQUE over a nullable conditional would not deduplicate
            rows.append((
                self.run_id,
                row.get('analysis'), row.get('grouped_by') or '',
                row.get('conditional') or '',
                row.get('n_subjects'), row.get('n_excluded'),
                row.get('n_variants'), row.get('n_independent'),
                row.get('n_asc'), row.get('threshold'),
                row.get('n_significant_variants'),
                row.get('n_independent_significant')))

        return self._insert('qtl_threshold',
                            ['run_id', 'analysis', 'grouped_by', 'conditional', 'n_subjects',
                             'n_excluded', 'n_variants', 'n_independent', 'n_asc',
                             'threshold', 'n_significant_variants',
                             'n_independent_significant'], rows)

    def load_usage_associations(self):
        """Variants, ASCs and the association rows, in one pass over the big file.

        The association file is the authority on which variants and ASCs this
        locus tested, so the reference tables are built from it rather than from
        the annotation files, which cover only some of them.

        It holds fewer variants than the run reports as tested, because some
        were genotyped but produced no usage association: the scan applies a
        complete-case subject set. Both numbers are kept, in qtl_variant and in
        qtl_threshold.n_variants, so a plot can say what its points are a subset of.
        """
        path = self._file(f'usage_associations_{self.locus}.tsv.gz')
        if not os.path.exists(path):
            return 0

        associations = []
        written = 0

        for row in _rows(path):
            variant = row['variant']
            if variant not in self.variants:
                cur = self.con.execute(
                    "INSERT INTO qtl_variant (run_id, variant, contig) VALUES (?,?,?)",
                    (self.run_id, variant,
                     # split on the FIRST underscore, not the last: ids are
                     # contig_pos, but a multi-allelic site adds a third part
                     # (chr22_23161341_2), and rsplit handed that whole
                     # `chr22_23161341` back as the contig. No contig carries an
                     # underscore, so the first one always ends it.
                     variant.split('_', 1)[0] if '_' in variant else None))
                self.variants[variant] = cur.lastrowid

            asc = row['asc']
            if asc not in self.ascs:
                cur = self.con.execute("INSERT INTO qtl_asc (run_id, asc) VALUES (?,?)",
                                       (self.run_id, asc))
                self.ascs[asc] = cur.lastrowid

            p = _num(row.get('p_value'))
            associations.append((
                self.variants[variant], self.ascs[asc],
                _num(row.get('beta')), _num(row.get('se')), _num(row.get('t_stat')),
                p, _neglog10(p), _bool(row.get('significant')),
                _num(row.get('distance_to_asc')),
                _bool(row.get('is_cis')), _bool(row.get('is_lead'))))

            if len(associations) >= BATCH:
                written += self._insert('qtl_usage_association', _ASSOC_COLS, associations)
                associations = []

        written += self._insert('qtl_usage_association', _ASSOC_COLS, associations)
        return written

    def annotate_variants(self):
        """Everything true of a variant rather than of one of its associations.

        The scan reports these per variant, so they arrive once here rather than
        repeated on each of its association rows.
        """
        path = self._file('variant_features.tsv.gz')
        if not os.path.exists(path):
            return 0

        updates = []
        for row in _rows(path):
            if row.get('locus') != self.locus:
                continue
            variant_id = self.variants.get(row['variant'])
            if variant_id:
                updates.append((
                    _num(row.get('pos'), int), _num(row.get('maf')),
                    _num(row.get('missing_rate')), _num(row.get('n'), int),
                    _num(row.get('min_genotype_group'), int),
                    _bool(row.get('well_powered')),
                    row.get('gene'), row.get('feature'), row.get('sub_feature'),
                    _num(row.get('distance_to_gene')), variant_id))

        self.con.executemany(
            "UPDATE qtl_variant SET pos=?, maf=?, missing_rate=?, n=?, "
            "min_genotype_group=?, well_powered=?, gene=?, feature=?, sub_feature=?, "
            "distance_to_gene=? WHERE id=?", updates)
        return len(updates)

    def annotate_ascs(self):
        """Where each cluster sits and how many alleles it holds."""
        path = self._file('asc_features.tsv.gz')
        if not os.path.exists(path):
            return 0

        updates = []
        for row in _rows(path):
            asc_id = self.ascs.get(row['asc'])
            if asc_id:
                updates.append((row.get('segment'), _num(row.get('asc_position')),
                                _num(row.get('asc_span')), _num(row.get('n_member'), int),
                                asc_id))

        self.con.executemany(
            "UPDATE qtl_asc SET segment=?, asc_position=?, asc_span=?, n_member=? "
            "WHERE id=?", updates)
        return len(updates)

    def summarise_variants(self):
        """Store each variant's strongest association, and which ASC gave it.

        The whole-locus Manhattan is one point per variant taken across every
        ASC. Computed on the fly that is 9,402 maxima out of 658,140 rows, and
        the database is immutable once built, so it was several seconds of the
        same arithmetic on every request. Fifty milliseconds here instead.

        Ties are broken arbitrarily, as they were before: two ASCs at the same
        p-value are two equally true answers to "which ASC", and the summary is
        not the place to invent a preference between them.
        """
        self.con.execute("""
            UPDATE qtl_variant SET
              best_neglog10_p = (SELECT MAX(neglog10_p) FROM qtl_usage_association u
                                 WHERE u.variant_id = qtl_variant.id),
              best_asc_id = (SELECT asc_id FROM qtl_usage_association u
                             WHERE u.variant_id = qtl_variant.id
                             ORDER BY neglog10_p DESC LIMIT 1),
              best_significant = (SELECT significant FROM qtl_usage_association u
                                  WHERE u.variant_id = qtl_variant.id
                                  ORDER BY neglog10_p DESC LIMIT 1)
            WHERE run_id = ?""", (self.run_id,))
        return self.con.execute(
            'SELECT COUNT(*) FROM qtl_variant WHERE run_id = ? AND best_neglog10_p IS NOT NULL',
            (self.run_id,)).fetchone()[0]

    def load_asc_usage(self):
        path = self._file('asc_usage.tsv.gz')
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
            cur = self.con.execute("INSERT INTO qtl_subject (run_id, subject) VALUES (?,?)",
                                   (self.run_id, subject))
            self.subjects[subject] = cur.lastrowid
        return self.subjects[subject]

    def load_dosage(self):
        """Genotypes, from the cohort matrix the run carries.

        Filtered to this locus by variant membership: the matrix holds the whole
        cohort at every locus, and this database wants one locus of it.
        """
        matrix = genotype_matrix(self.run_dir)
        return self._load_dosage_matrix(matrix) if matrix else 0

    def _load_dosage_matrix(self, path):
        """Read the wide genotype matrix, streaming one variant per row.

        Filtered on both axes as it goes rather than read whole: it holds the
        entire cohort at every locus, and this database wants one locus of it.
        """
        rows = []
        written = 0

        with gzip.open(path, 'rt', newline='') as handle:
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

    def load_pairing(self):
        """The D/J pairing scans. Only IGH has them, and only for its variants."""
        totals = {}

        path = self._file('pairing.tsv.gz')
        if os.path.exists(path):
            rows = []
            for row in _rows(path):
                variant_id = self.variants.get(row['variant'])
                if not variant_id:
                    continue
                p = _num(row.get('p_value'))
                rows.append((row.get('conditional'), variant_id,
                             row.get('gene'), _num(row.get('n'), int),
                             _num(row.get('pillai')), _num(row.get('f_stat')),
                             p, _neglog10(p), _num(row.get('min_genotype_group'), int),
                             _bool(row.get('significant'))))
            totals['pairing_associations'] = self._insert(
                'qtl_pairing_association', _PAIRING_COLS, rows)

        path = self._file('cell_tests.tsv.gz')
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
        path = self._file('dj_enrichment.tsv.gz')
        if totals.get('pairing_associations') and os.path.exists(path) and self.subjects:
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

    def finish(self, species, project=None):
        # details describes the file, not the study: rewritten so the timestamp
        # tracks the most recent build
        self.con.execute('DELETE FROM details')
        self.con.execute(
            "INSERT INTO details (dbtype, species, locus, created_on, created_by) "
            "VALUES (?,?,?,datetime('now'),?)",
            ('guQTL', species, self.locus, 'make_qtl_db'))
        self.con.commit()
        self.con.execute('PRAGMA optimize')
        self.con.close()

        with open(os.path.join(os.path.dirname(self.path) or '.', 'db_description.txt'), 'w') as fo:
            fo.write(f'Gene-usage QTL results for {species} {self.locus}'
                     + (f', project {project}' if project else ''))


_ASSOC_COLS = ['variant_id', 'asc_id', 'beta', 'se', 't_stat', 'p_value',
               'neglog10_p', 'significant', 'distance_to_asc', 'is_cis', 'is_lead']
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


def build(run_dir, species, locus, path, project=None, manifest=None):
    """Add one study's locus to the database at `path`, returning row counts.

    `project` is the study whose cohort this run scanned. It is stated, not
    inferred: reading an identity out of a filename is how the wrong one gets in.
    """
    builder = QtlBuilder(run_dir, locus, path, manifest)

    _READ.clear()
    counts = {}
    builder.load_run(project)
    counts['thresholds'] = builder.load_thresholds()
    counts['usage_associations'] = builder.load_usage_associations()
    counts['variants'] = len(builder.variants)
    counts['ascs'] = len(builder.ascs)
    counts['annotated'] = builder.annotate_variants()
    counts['clusters_annotated'] = builder.annotate_ascs()
    counts['summarised'] = builder.summarise_variants()
    counts['asc_usage'] = builder.load_asc_usage()
    counts['dosage'] = builder.load_dosage()
    counts['subjects'] = len(builder.subjects)
    counts.update(builder.load_pairing())
    builder.finish(species, project)

    check_schema(run_dir, manifest)
    return counts
