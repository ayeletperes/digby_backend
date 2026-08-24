"""Install real study databases into static/study_data from a digby_private_data checkout.

The backend expects        static/study_data/VDJbase/db/<Species>/<Locus>/db.sqlite3
                           static/study_data/Genomic/db/<Species>/<Locus>/db.sqlite3

digby_private_data holds   AIRR-seq/<Species>/<Locus>/db.sqlite3
                           Genomic/<Species>/<Locus>/db.sqlite3

so the AIRR-seq directory has to be mapped onto VDJbase/db.

AIRR-seq databases are COPIED, not linked. study_data_db_init adds an `asc_genotype`
column to the Sample table when it is missing, which rewrites the file - and the
databases published in digby_private_data do not have that column. Linking them
would modify tracked files in that repository on the first start.

Genomic databases are symlinked: nothing writes to them, and they are ~350MB.

Usage, from the repository root:

    .venv/bin/python dev_fixtures/install_study_data.py ../../digby_private_data
"""

import os
import sys
import shutil

FLAVOURS = [
    # (source directory, destination under study_data, copy?)
    ('AIRR-seq', os.path.join('VDJbase', 'db'), True),
    ('Genomic', os.path.join('Genomic', 'db'), False),
]


def install(source_root, static_root):
    if not os.path.isdir(source_root):
        raise SystemExit(f'No such directory: {source_root}')

    installed = 0

    for source_name, dest_name, copy in FLAVOURS:
        source_dir = os.path.join(source_root, source_name)
        if not os.path.isdir(source_dir):
            print(f'  skipping {source_name}: not present in {source_root}')
            continue

        dest_root = os.path.join(static_root, 'study_data', dest_name)

        for species in sorted(os.listdir(source_dir)):
            species_dir = os.path.join(source_dir, species)
            if not os.path.isdir(species_dir):
                continue

            for locus in sorted(os.listdir(species_dir)):
                database = os.path.join(species_dir, locus, 'db.sqlite3')
                if not os.path.isfile(database):
                    continue

                target_dir = os.path.join(dest_root, species, locus)
                os.makedirs(target_dir, exist_ok=True)
                target = os.path.join(target_dir, 'db.sqlite3')

                if os.path.islink(target) or os.path.exists(target):
                    os.remove(target)

                if copy:
                    shutil.copy2(database, target)
                else:
                    os.symlink(os.path.abspath(database), target)

                description = os.path.join(species_dir, locus, 'db_description.txt')
                if os.path.isfile(description):
                    shutil.copy2(description, target_dir)

                print(f'  {"copied " if copy else "linked "} {species}/{locus}'
                      f'  ->  {dest_name}')
                installed += 1

    print(f'\n{installed} databases installed under {os.path.join(static_root, "study_data")}')
    if installed:
        print('AIRR-seq databases were copied, so the originals are never modified.')


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)

    static_root = os.path.join(os.getcwd(), 'static')
    if not os.path.isdir(static_root):
        raise SystemExit('Run this from the repository root (no ./static here).')

    install(sys.argv[1], static_root)


if __name__ == '__main__':
    main()
