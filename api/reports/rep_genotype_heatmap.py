# Haplotype heatmap for VDJbase samples

from werkzeug.exceptions import BadRequest

from api.reports.report_utils import trans_df, collate_samples, chunk_list
from api.reports.reports import SYSDATA, run_rscript, send_report
from api.reports.report_utils import make_output_file, find_primer_translations, translate_primer_alleles, translate_primer_genes
from app import app, vdjbase_dbs
from db.vdjbase_model import Sample, Gene
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params, get_multiple_order_file
import pandas as pd


HEATMAP_GENOTYPE_SCRIPT = "genotype_heatmap.R"
SAMPLE_CHUNKS = 400


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    html = (format == 'html')
    chain, samples_by_dataset = collate_samples(rep_samples)
    genotypes = pd.DataFrame()

    for dataset in samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        primer_trans, gene_subs = find_primer_translations(session)

        sample_list = []
        for sample_chunk in chunk_list(samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list.extend(session.query(Sample.name, Sample.genotype, Sample.patient_id).filter(Sample.name.in_(sample_chunk)).all())

        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)

        if len(wanted_genes) > 0:
            for (name, genotype, patient_id) in sample_list:
                sample_path = os.path.join(VDJBASE_SAMPLE_PATH, species, dataset, genotype.replace('samples/', ''))

                if not os.path.isfile(sample_path):
                    raise BadRequest('Genotype file for sample %s/%s is missing.' % (dataset, name))

                genotype = pd.read_csv(sample_path, sep='\t', dtype=str)

                genotype = trans_df(genotype)

                # translate pipeline allele names to VDJbase allele names
                for col in ['alleles', 'GENOTYPED_ALLELES']:
                    genotype[col] = [translate_primer_alleles(x, y, primer_trans) for x, y in zip(genotype['gene'], genotype[col])]

                genotype['gene'] = [translate_primer_genes(x, gene_subs) for x in genotype['gene']]
                genotype = genotype[genotype.gene.isin(wanted_genes)]

                subject_name = name if len(samples_by_dataset) == 1 else dataset + '_' + name

                if 'subject' not in genotype.columns.values:
                    genotype.insert(0, 'subject', subject_name)
                else:
                    genotype.subject = subject_name

                genotypes = genotypes.append(genotype)[genotype.columns.tolist()]

    if len(genotypes) == 0:
        raise BadRequest('No records matching the filter criteria were found.')

    geno_path = make_output_file('csv')
    genotypes.to_csv(geno_path, sep='\t')

    if format == 'pdf':
        attachment_filename = '%s_genotype.pdf' % species
    else:
        attachment_filename = None

    locus_order = ('sort_order' in params and params['sort_order'] == 'Locus')
    gene_order_file = get_multiple_order_file(species, samples_by_dataset.keys(), locus_order=locus_order)

    output_path = make_output_file('html' if html else 'pdf')
    file_type = 'T' if html else 'F'
    cmd_line = ["-i", geno_path,
                "-o", output_path,
                "-t", file_type,
                "-k", str(params['f_kdiff']),
                "-c", chain,
                "-g", gene_order_file
                ]

    if run_rscript(HEATMAP_GENOTYPE_SCRIPT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, attachment_filename)
    else:
        raise BadRequest('No output from report')



