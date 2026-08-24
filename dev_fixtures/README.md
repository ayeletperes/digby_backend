# Local development data

Two ways to get data into a fresh clone: the **real databases**, or a tiny
**synthetic tree** for testing. They occupy the same directory, so use one at a time.

---

## Real study data

Nothing in this repository says where the databases come from. They live in a
separate, private repository:

    https://github.com/ayeletperes/digby_private_data

Its layout does not match what the backend expects, so it cannot simply be copied:

| digby_private_data | backend expects |
|---|---|
| `AIRR-seq/<Species>/<Locus>/db.sqlite3` | `static/study_data/VDJbase/db/<Species>/<Locus>/db.sqlite3` |
| `Genomic/<Species>/<Locus>/db.sqlite3` | `static/study_data/Genomic/db/<Species>/<Locus>/db.sqlite3` |

`install_study_data.py` does the mapping:

```bash
.venv/bin/python dev_fixtures/install_study_data.py /path/to/digby_private_data
```

That installs 14 databases: Human IGH/IGK/IGL/TRB and Rhesus Macaque IGH/IGK/IGL
for AIRR-seq, plus Human IGH/IGHC/IGK/IGL and Rhesus IGH/IGK/IGL for genomic.

**AIRR-seq databases are copied rather than symlinked, deliberately.**
`study_data_db_init` adds an `asc_genotype` column to the `Sample` table when it is
missing, which rewrites the file — and none of the published AIRR-seq databases have
that column. Symlinking them would modify tracked files in `digby_private_data` on the
first start. Genomic databases are symlinked; nothing writes to them and they are ~350MB.

Sample files (FASTA/BAM/GFF for the genome browser) are a separate, much larger set
and are not needed for the API or the dashboard.

---

## Synthetic fixtures

The real study databases are large and are not in the repository, so a fresh clone
has no data and no endpoint can be exercised. `make_test_data.py` writes the
smallest synthetic `static/study_data` tree that still covers the cases the API
branches on:

| Dataset | Held in | Exercises |
|---|---|---|
| Human / IGH | both databases | results merged across databases; a locus listed once |
| Human / IGL | AIRR-seq only | the genomic side absent |
| Rhesus Macaque / IGH | genomic only | a species absent from `vdjbase_dbs` entirely |

It also seeds a pseudogene, an orphon and a `*Del` allele, each of which the
endpoints are expected to filter out.

## Use

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt
cp sample_secret.cfg secret.cfg          # then edit if you need real credentials
.venv/bin/python dev_fixtures/make_test_data.py
.venv/bin/python tests/test_refbook.py
```

Both must be run from the repository root: the app takes `BASE_PATH` from the
working directory.

## Safety

The generator writes a `.synthetic_fixture` marker into each tree it creates and
**refuses to delete any tree without one**, so running it in a checkout that holds
real databases fails loudly instead of destroying them.
