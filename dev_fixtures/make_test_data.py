"""Build a small synthetic study_data tree so the API can be exercised offline.

The real databases are large and are not in the repo, which makes even a smoke
test of an endpoint impossible on a fresh clone. This writes the smallest tree
that still covers the cases the refbook endpoints actually branch on:

  Human / IGH            in BOTH databases  - the merge path, and locus de-duplication
  Human / IGL            AIRR-seq only      - the genomic side missing
  Rhesus Macaque / IGH   genomic only       - a species absent from vdjbase_dbs entirely

Run from the repo root:  .venv/bin/python dev_fixtures/make_test_data.py
"""

import os
import sys
import shutil
import datetime

# run from anywhere; the models are imported as top-level packages
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from db.vdjbase_model import (Base as VdjBase, Details as VdjDetails,
                              Gene as VdjGene, Allele as VdjAllele,
                              AllelesSample as VdjAllelesSample)
from db.vdjbase_airr_model import Sample as VdjSample
from db.genomic_db import (Base as GenBase, Details as GenDetails,
                           Gene as GenGene, Sequence as GenSequence)
# imported for its side effect: RefSeq relates to Sample, so the mapper cannot
# be configured until both AIRR metadata modules are in the class registry
import db.genomic_airr_model  # noqa: F401

STATIC = os.path.join(os.getcwd(), 'static', 'study_data')

# written beside every tree this script creates. Nothing without it is ever deleted,
# so pointing this at a checkout that holds real databases refuses rather than wipes.
MARKER = '.synthetic_fixture'

# a V-REGION long enough to look real; gapped form carries IMGT dots
SEQ_GAPPED = ('CAGGTGCAGCTGGTGCAGTCTGGGGCT...GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTC'
              'TCCTGCAAGGCTTCTGGATACACCTTC............ACCGGCTACTATATGCAC')
SEQ = SEQ_GAPPED.replace('.', '')


def _session(path, base):
    if os.path.exists(path):
        os.remove(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    engine = create_engine('sqlite:///' + path)
    base.metadata.create_all(engine)
    return sessionmaker(bind=engine)()


def _describe(directory, text):
    with open(os.path.join(directory, 'db_description.txt'), 'w') as fo:
        fo.write(text)


def _clear(db_root):
    """Remove a previously generated tree, refusing to touch anything else."""
    if not os.path.exists(db_root):
        return
    if not os.path.exists(os.path.join(db_root, MARKER)):
        raise SystemExit(
            f'Refusing to delete {db_root}: it holds databases this script did not create.\n'
            f'Move them aside first, or run against a checkout with no real study data.')
    shutil.rmtree(db_root)


def build_vdjbase(species, locus, genes, pseudo_genes=(), orphons=()):
    directory = os.path.join(STATIC, 'VDJbase', 'db', species, locus)
    session = _session(os.path.join(directory, 'db.sqlite3'), VdjBase)

    session.add(VdjDetails(dbtype='AIRR-seq', species=species, locus=locus,
                           created_on=datetime.datetime.now(), created_by='fixture'))

    # one sample, with dummy ids for the metadata tables we do not populate
    session.add(VdjSample(id=1, sample_name=f'{species[:3]}_{locus}_S1', patient_id=1,
                          seq_protocol_id=1, study_id=1, tissue_pro_id=1, data_pro_id=1))

    allele_id = 0
    for gene_id, (gene_name, pseudo) in enumerate(
            [(g, 0) for g in genes] + [(g, 1) for g in pseudo_genes] + [(g, 0) for g in orphons], start=1):
        session.add(VdjGene(id=gene_id, name=gene_name, type=locus + 'V',
                            family=gene_name.split('-')[0], species=species, pseudo_gene=pseudo))

        for suffix in ('01', '02'):
            allele_id += 1
            session.add(VdjAllele(id=allele_id, name=f'{gene_name}*{suffix}',
                                  pipeline_name=f'{gene_name}*{suffix}', seq=SEQ_GAPPED,
                                  seq_len=str(len(SEQ)), similar='', appears=5,
                                  gene_id=gene_id, is_single_allele=True,
                                  low_confidence=False, novel=(suffix == '02'), max_kdiff=0.0))
            session.add(VdjAllelesSample(allele_id=allele_id, patient_id=1, sample_id=1,
                                         hap='geno', kdiff=0.0, count=10 * allele_id,
                                         total_count=1000))

    # a deleted allele: asc_seqs must exclude it
    allele_id += 1
    session.add(VdjAllele(id=allele_id, name=f'{genes[0]}*Del', pipeline_name='Del',
                          seq=SEQ_GAPPED, seq_len=str(len(SEQ)), similar='', appears=1,
                          gene_id=1, is_single_allele=True, low_confidence=False,
                          novel=False, max_kdiff=0.0))

    session.commit()
    _describe(directory, f'Synthetic AIRR-seq fixture: {species} {locus}')
    print(f'  VDJbase  {species}/{locus}: {len(genes)} genes '
          f'(+{len(pseudo_genes)} pseudo, +{len(orphons)} orphon)')


def build_genomic(species, locus, genes, pseudo_genes=()):
    directory = os.path.join(STATIC, 'Genomic', 'db', species, locus)
    session = _session(os.path.join(directory, 'db.sqlite3'), GenBase)

    session.add(GenDetails(dbtype='Genomic', species=species, locus=locus,
                           created_on=datetime.datetime.now(), created_by='fixture'))

    seq_id = 0
    for gene_id, (gene_name, pseudo) in enumerate(
            [(g, 0) for g in genes] + [(g, 1) for g in pseudo_genes], start=1):
        # the genomic Gene has no species column - species is implied by the dataset path
        session.add(GenGene(id=gene_id, name=gene_name, type=locus + 'V',
                            family=gene_name.split('-')[0], pseudo_gene=pseudo))

        for suffix in ('01', '02'):
            seq_id += 1
            session.add(GenSequence(id=seq_id, name=f'{gene_name}*{suffix}',
                                    imgt_name=f'{gene_name}*{suffix}', type=locus + 'V',
                                    novel=(suffix == '02'), appearances=3, deleted=False,
                                    functional='Functional', sequence=SEQ,
                                    gapped_sequence=SEQ_GAPPED, gene_id=gene_id))

    session.commit()
    _describe(directory, f'Synthetic genomic fixture: {species} {locus}')
    print(f'  Genomic  {species}/{locus}: {len(genes)} genes (+{len(pseudo_genes)} pseudo)')


def main():
    for flavour in ('VDJbase', 'Genomic'):
        _clear(os.path.join(STATIC, flavour, 'db'))

    for flavour in ('VDJbase', 'Genomic'):
        db_root = os.path.join(STATIC, flavour, 'db')
        os.makedirs(db_root, exist_ok=True)
        open(os.path.join(db_root, MARKER), 'w').close()

    print('Building synthetic study data...')
    # Human IGH: in both databases, with only IGHV3-23 shared
    build_vdjbase('Human', 'IGH', ['IGHV1-2', 'IGHV3-23'],
                  pseudo_genes=['IGHV3-30'], orphons=['IGHV1-2/OR15-1'])
    build_genomic('Human', 'IGH', ['IGHV3-23', 'IGHV4-34'])
    # Human IGL: AIRR-seq only
    build_vdjbase('Human', 'IGL', ['IGLV1-40', 'IGLV2-14'])
    # Rhesus IGH: genomic only, and this species is absent from the AIRR-seq tree
    build_genomic('Rhesus Macaque', 'IGH', ['IGHV1-2', 'IGHV4-34'])
    print('Done. Root:', STATIC)


if __name__ == '__main__':
    main()
