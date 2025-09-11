import simple_bio_seq as simple
import csv
import os

refs = simple.read_imgt_fasta('imgt_gapped_21_03_2021.fa', ['Homo sapiens'], chains=['IGLV', 'IGLJ'])
refs = refs['Homo sapiens']

def gene_number(name):
    if 'OR' in name:
        return 900

    if '-' in name:
        name = name.split('-')[1]
    else:
        name = name[4:]

    if '*' in name:
        name = name.split('*')[0]
        
    try:
        ret = int(name)
    except:
        ret = 900

    return ret


def allele_number(name):
    try:
        ret = int(name.split('*')[1])
    except:
        ret = 50

    return ret
    
    
def receptor_chain(name):
    return name[:4]


def sortkey(name):
    cname = receptor_chain(name)
    if 'V' in cname:
        c = 200000
    elif 'D' in cname:
        c = 100000
    else:
        c = 0
        
    return c + 100*gene_number(name) + allele_number(name)
    

if os.path.isfile('changed_refrence.csv'):
    with open('changed_refrence.csv', 'r') as fi:
        reader = csv.DictReader(fi)
        for row in reader:
            if row['OLD_GENE'] in refs['IGKV']:
                if 'DELETED' in row['NEW_GENE']:
                    del refs['IGKV'][row['OLD_GENE']]
                else:
                    refs['IGKV'][row['NEW_GENE']] = refs['IGKV'][row['OLD_GENE']]
                    del refs['IGKV'][row['OLD_GENE']]
                    

seqnames = []

for chain in refs.keys():
    simple.write_fasta(refs[chain], chain + '.fasta')
    seqnames.extend(refs[chain].keys())

seqnames.sort(key=sortkey)

with open('gene_order.py', 'w') as fo:
    fo.write('LOCUS_ORDER = [\n')
    for seq in seqnames:
        fo.write('"%s",\n' % seq)
    fo.write(']\n')
    fo.write('ALPHA_ORDER = [\n')
    for seq in seqnames:
        fo.write('"%s",\n' % seq)
    fo.write(']\n')
    fo.write('PSEUDO_GENES = [\n')
    for seq in seqnames:
        if 'OR' in seq or 'P' in seq:
            fo.write('"%s",\n' % seq)
    fo.write(']\n')
    
    

    