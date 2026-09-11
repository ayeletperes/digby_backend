"""Make a guQTL sqlite database from files in the current directory.

    python make_qtl_db.py <species> <locus>

Run from the dataset directory, as make_vdjbase_db.py and make_genomic_db.py
are. A yml file there names the studies to build and where each one's run
directory is; every study goes into the one db.sqlite3, keyed by its run.
"""

import os
import sys

import yaml

from db.qtl_maint import build


def read_yml_file(dataset_dir):
    yml_files = [entry.name for entry in os.scandir(dataset_dir)
                 if entry.is_file() and os.path.splitext(entry.name)[1] in ('.yml', '.yaml')]
    if not yml_files:
        sys.exit(f'Error: no yml file found in directory {dataset_dir}.')
    if len(yml_files) > 1:
        sys.exit(f'Error: multiple yml files found in directory {dataset_dir}.')
    with open(os.path.join(dataset_dir, yml_files[0])) as fi:
        return yaml.safe_load(fi)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    species, locus = sys.argv[1], sys.argv[2]

    dataset_dir = os.getcwd()
    studies = (read_yml_file(dataset_dir) or {}).get('Studies') or {}
    if not studies:
        sys.exit('The yml file names no studies under `Studies:`.')

    db_file = os.path.join(dataset_dir, 'db.sqlite3')

    for project in sorted(studies):
        run_dir = studies[project].get('run_dir') or os.path.join('studies', project)
        if not os.path.isdir(run_dir):
            sys.exit(f'{project}: no such run directory: {run_dir}')

        # the manifest sits beside the database, named for its study, as the
        # repertoire sets keep their MiAIRR json
        manifest = studies[project].get('manifest_file') or f'{project}_manifest.json'
        if not os.path.isfile(manifest):
            sys.exit(f'{project}: no such manifest: {manifest}')

        print(f'{species} {locus} [{project}]:')
        counts = build(run_dir, species, locus, db_file, project=project,
                       manifest=manifest)
        for name, value in counts.items():
            print(f'    {name:22} {value:>9,}')

    print(f'    -> {db_file} ({os.path.getsize(db_file) / 1024 / 1024:.0f} MB), '
          f'{len(studies)} stud{"y" if len(studies) == 1 else "ies"}')


if __name__ == '__main__':
    main()
