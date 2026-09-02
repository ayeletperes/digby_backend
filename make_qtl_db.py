"""Build the guQTL databases from an igqtl.R run directory.

    python make_qtl_db.py <run_dir> [--species Human] [--locus IGH] [--project P28]

<run_dir> is a dated directory under results/igqtl (or its `current` symlink).
With no --locus, every locus the run produced is built.

--project names the study whose cohort the run scanned. One database holds one
project, and the dashboard offers whichever ones are built, so a second study is
a second build rather than a merge.
"""

import argparse
import os
import sys

from db.qtl_maint import build, loci_in


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('run_dir', help='an igqtl run directory')
    parser.add_argument('--species', default='Human')
    parser.add_argument('--locus', action='append',
                        help='build only this locus; repeatable')
    parser.add_argument('--static', default=None,
                        help='static path to write under (default ./static)')
    parser.add_argument('--genotypes', default=None,
                        help='the cohort genotype matrix to take genotypes from '
                             '(default: the one the run records in its manifest)')
    parser.add_argument('--project', default=None,
                        help='the study whose cohort this run scanned, e.g. P28. '
                             'Recorded in the database and offered as a choice in '
                             'the dashboard. Not guessed from the run: state it')
    args = parser.parse_args()

    if not os.path.isdir(args.run_dir):
        sys.exit(f'No such run directory: {args.run_dir}')

    if args.genotypes and not os.path.exists(args.genotypes):
        sys.exit(f'No such genotype matrix: {args.genotypes}')

    static_path = args.static or os.path.join(os.getcwd(), 'static')
    loci = args.locus or loci_in(args.run_dir)
    if not loci:
        sys.exit(f'No usage_associations_*.tsv.gz found under {args.run_dir}/source_data')

    for locus in loci:
        print(f'{args.species} {locus}' + (f' [{args.project}]' if args.project else '') + ':')
        counts = build(args.run_dir, args.species, locus, static_path,
                       args.genotypes, args.project)
        path = counts.pop('path')
        for name, value in counts.items():
            print(f'    {name:22} {value:>9,}')
        print(f'    -> {path} ({os.path.getsize(path) / 1024 / 1024:.0f} MB)')


if __name__ == '__main__':
    main()
