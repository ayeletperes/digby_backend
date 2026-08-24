"""Regression tests for the refbook (Reference Dashboard) endpoints.

Runs against the synthetic tree written by dev_fixtures/make_test_data.py, which
is shaped to cover the cases these endpoints branch on:

    Human / IGH            both databases     - results must be merged, locus listed once
    Human / IGL            AIRR-seq only
    Rhesus Macaque / IGH   genomic only       - species absent from vdjbase_dbs entirely

Run from the repo root:

    .venv/bin/python dev_fixtures/make_test_data.py
    .venv/bin/python tests/test_refbook.py
"""

import os
import sys
import logging
import warnings

warnings.filterwarnings('ignore')
logging.disable(logging.CRITICAL)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import app as appmod  # noqa: E402

# These assertions describe the synthetic tree, not the real databases, so refuse
# to report a misleading failure when a developer has real study data installed.
MARKER = os.path.join(appmod.app.config['STATIC_PATH'],
                      'study_data', 'VDJbase', 'db', '.synthetic_fixture')
if not os.path.exists(MARKER):
    print('SKIP: these tests need the synthetic tree.\n'
          '      Run: .venv/bin/python dev_fixtures/make_test_data.py\n'
          '      (real study data appears to be installed; move it aside first)')
    sys.exit(0)

client = appmod.app.test_client()
failures = []


def check(label, url, predicate, describe):
    try:
        response = client.get(url)
    except Exception as e:
        # an endpoint that raises must be one reported failure, not the end of the run
        failures.append(f'{label}: {type(e).__name__}: {e}')
        print(f'  FAIL  {label}  ({type(e).__name__}: {e})')
        return

    if response.status_code != 200:
        failures.append(f'{label}: HTTP {response.status_code} from {url}')
        print(f'  FAIL  {label}  (HTTP {response.status_code})')
        return

    body = response.get_json()
    if predicate(body):
        print(f'  ok    {label}')
    else:
        failures.append(f'{label}: expected {describe}, got {body}')
        print(f'  FAIL  {label}  (expected {describe})')


def main():
    print('refbook endpoints:')

    # A locus present in both databases was appended once per database.
    check('species_and_loci lists each locus once',
          '/api/refbook/species_and_loci',
          lambda b: b['loci']['Human'] == ['IGH', 'IGL'],
          "Human loci == ['IGH', 'IGL']")

    # The ASC list was re-derived from the last query alone, dropping the other
    # database's genes; IGHV1-2 is AIRR-seq only and IGHV4-34 genomic only.
    check('ascs_in_locus merges both databases',
          '/api/refbook/ascs_in_locus/Human/IGH',
          lambda b: b['ascs'] == ['IGHV1-2', 'IGHV3-23', 'IGHV4-34'],
          'all three ASCs')

    check('ascs_in_locus reports both source flags',
          '/api/refbook/ascs_in_locus/Human/IGH',
          lambda b: b['genomic'] is True and b['airr_seq'] is True,
          'genomic and airr_seq both true')

    check('ascs_in_locus excludes pseudogenes and orphons',
          '/api/refbook/ascs_in_locus/Human/IGH',
          lambda b: not any('OR' in a or a == 'IGHV3-30' for a in b['ascs']),
          'no IGHV3-30 and no /OR gene')

    check('ascs_in_locus handles an AIRR-seq-only locus',
          '/api/refbook/ascs_in_locus/Human/IGL',
          lambda b: b['ascs'] == ['IGLV1-40', 'IGLV2-14'] and b['genomic'] is False,
          'the two IGL genes, genomic false')

    # These four raised KeyError for a species held in only one of the two
    # databases, because dbs[species] was indexed before it was known to exist.
    genomic_only = [
        ('ascs_in_locus', '/api/refbook/ascs_in_locus/Rhesus Macaque/IGH',
         lambda b: b['ascs'] == ['IGHV1-2', 'IGHV4-34'], 'the two genomic ASCs'),
        ('ascs_overview', '/api/refbook/ascs_overview/Rhesus Macaque/IGH/IGHV1-2',
         lambda b: b['total'] == 2, 'two alleles'),
        ('asc_seqs', '/api/refbook/asc_seqs/Rhesus Macaque/IGH/IGHV1-2',
         lambda b: len(b['alleles']) == 2, 'two sequences'),
        ('asc_usage', '/api/refbook/asc_usage/Rhesus Macaque/IGH/IGHV1-2',
         lambda b: b['alleles'] == [], 'an empty list, not an error'),
        ('asc_zygousity', '/api/refbook/asc_zygousity/Rhesus Macaque/IGH/IGHV1-2',
         lambda b: b['samples'] == [], 'an empty list, not an error'),
    ]
    for name, url, predicate, describe in genomic_only:
        check(f'{name} survives a genomic-only species', url, predicate, describe)

    check('asc_seqs excludes deleted alleles',
          '/api/refbook/asc_seqs/Human/IGH/IGHV1-2',
          lambda b: not any('Del' in a['name'] for a in b['alleles']),
          'no *Del allele')

    check('asc_seqs returns gapped and ungapped forms',
          '/api/refbook/asc_seqs/Human/IGH/IGHV1-2',
          lambda b: all('.' in a['seq_gapped'] and '.' not in a['seq'] for a in b['alleles']),
          'seq_gapped dotted, seq undotted')

    # An unknown species or locus is a 404, not a crash.
    for url in ('/api/refbook/ascs_in_locus/Nonesuch/IGH',
                '/api/refbook/ascs_in_locus/Human/XYZ'):
        try:
            code = client.get(url).status_code
        except Exception as e:
            failures.append(f'{url}: {type(e).__name__}: {e}')
            print(f'  FAIL  {url} raised {type(e).__name__}')
            continue
        if code == 404:
            print(f'  ok    404 for {url}')
        else:
            failures.append(f'{url}: expected 404, got {code}')
            print(f'  FAIL  expected 404 for {url}, got {code}')

    print()
    if failures:
        print(f'{len(failures)} FAILED')
        for failure in failures:
            print('  -', failure)
        return 1

    print('all passed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
