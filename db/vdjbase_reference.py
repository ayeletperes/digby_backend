import os.path
import os
from Bio import SeqIO
import importlib.util

from db.vdjbase_model import Gene, Allele, AlleleAliasSets
from db.vdjbase_exceptions import DbCreationError


def read_reference(filename):
    records = {}

    for rec in SeqIO.parse(filename, 'fasta'):
        if '|' in rec.description:
            rd = rec.description.split('|')[1]      # assume IMGT convention
        else:
            rd = rec.description

        records[rd] = rec.seq.lower()

    return records


def import_reference_alleles(reference_dir, session, species):
    result = []
    if os.path.isfile(os.path.join(reference_dir, 'gene_order.py')):
        spec = importlib.util.spec_from_file_location("gene_order", os.path.join(reference_dir, 'gene_order.py'))
        gene_order = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(gene_order)
    else:
        raise DbCreationError('gene_order.py not found - skipped')

    added_genes = []
    extra_locus = len(gene_order.LOCUS_ORDER)
    extra_alpha = len(gene_order.ALPHA_ORDER)
    for file in os.listdir(reference_dir):
        if os.path.splitext(file)[1] == '.fasta':
            recs = read_reference(os.path.join(reference_dir, file))

            for allele, sequence in recs.items():
                gene_name = allele.split('*')[0] if '*' in allele else allele

                if gene_name not in added_genes:
                    (extra_alpha, extra_locus) = add_gene(extra_alpha, extra_locus, gene_name, gene_order, session, species)
                    added_genes.append(gene_name)

                save_allele(allele, gene_name, sequence, session)

        session.commit()

    if len(added_genes) == 0:
        raise DbCreationError('No genes added from reference set - skipped')

    result.append('Reference alleles added')
    return result


def save_allele(allele_name, gene_name, sequence, session):
    similar = session.query(Allele).filter(Allele.seq == str(sequence)).one_or_none()

    if similar is not None:
        if similar.similar is None or len(similar.similar) == 0:
            similar.similar = '|%s|' % allele_name
        else:
            similar.similar += ', ' + '|%s|' % allele_name
    else:
        g = session.query(Gene).filter(Gene.name == gene_name).one_or_none()
        a = Allele(
            name=allele_name,
            seq=str(sequence),
            seq_len=str(len(sequence)),
            appears=0,
            gene_id=g.id,
            is_single_allele=True,
            low_confidence=False,
            novel=False,
            max_kdiff=0,
            similar='',
            pipeline_name='',
        )
        session.add(a)
    session.flush()


def add_gene(extra_alpha, extra_locus, gene, gene_order, session, species):
    if gene in gene_order.LOCUS_ORDER:
        locus_order = gene_order.LOCUS_ORDER.index(gene)
    else:
        locus_order = extra_locus
        extra_locus += 1
    if gene in gene_order.ALPHA_ORDER:
        alpha_order = gene_order.ALPHA_ORDER.index(gene)
    else:
        alpha_order = extra_alpha
        extra_alpha += 1
    g = Gene(
        name=gene,
        type=gene[0:4],
        family=gene.split('-')[0] if '-' in gene else gene,
        species=species,
        locus_order=locus_order,
        alpha_order=alpha_order,
        pseudo_gene=1 if gene in gene_order.PSEUDO_GENES else 0
    )
    session.add(g)
    return (extra_alpha, extra_locus)


# scan subdirs of reference_dir for alias sets
# add an alias set if we find one or more .fasta files in a subdir
def add_alias_sets(reference_dir, session):
    results = []
    alias_num = 1

    alleles = session.query(Allele).all()
    ungapped_seqs = {a.seq.replace('.', ''): a for a in alleles}

    for subdir in os.listdir(reference_dir):
        subdir_path = os.path.join(reference_dir, subdir)
        if os.path.isdir(subdir_path):
            fasta_files = [f for f in os.listdir(subdir_path) if f.endswith('.fasta')]
            if fasta_files:
                alias_set = AlleleAliasSets(
                    set_name=subdir,
                    alias_number=alias_num
                )
                session.add(alias_set)
                alias_name = f'alias_{alias_num}'
                alias_num += 1
                results.append(f'Processing alias set {alias_num}')

                for fasta_file in fasta_files:
                    recs = read_reference(os.path.join(reference_dir, subdir_path, fasta_file))

                    for allele, sequence in recs.items():
                        if sequence in ungapped_seqs:
                            similar = ungapped_seqs[sequence]

                            # add the allele name to the alias_name column of the similar allele
                            sn = getattr(similar, alias_name)
                            if sn is None or len(sn) == 0:
                                setattr(similar, alias_name, allele)
                            else:
                                setattr(similar, alias_name, f'{sn}|{allele}')

                session.commit()

    results.append('Alias sets processed')
    return results
 