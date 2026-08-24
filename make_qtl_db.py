"""Build the guQTL databases from an igqtl.R run directory.

    python make_qtl_db.py <run_dir> [--species Human] [--locus IGH]

<run_dir> is a dated directory under results/igqtl (or its `current` symlink).
With no --locus, every locus the run produced is built.
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
    args = parser.parse_args()

    if not os.path.isdir(args.run_dir):
        sys.exit(f'No such run directory: {args.run_dir}')

    static_path = args.static or os.path.join(os.getcwd(), 'static')
    loci = args.locus or loci_in(args.run_dir)
    if not loci:
        sys.exit(f'No usage_associations_*.tsv.gz found under {args.run_dir}/source_data')

    for locus in loci:
        print(f'{args.species} {locus}:')
        counts = build(args.run_dir, args.species, locus, static_path)
        path = counts.pop('path')
        for name, value in counts.items():
            print(f'    {name:22} {value:>9,}')
        print(f'    -> {path} ({os.path.getsize(path) / 1024 / 1024:.0f} MB)')


if __name__ == '__main__':
    main()
