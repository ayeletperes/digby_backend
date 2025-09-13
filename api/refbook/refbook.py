from app import vdjbase_dbs, app, genomic_dbs
from flask import request
from flask_restx import Resource, reqparse
from api.restx import api
from api.system.system import digby_protected
from sqlalchemy import or_, func, cast, Float

from api.vdjbase.vdjbase import get_vdjbase_species, find_datasets as get_vdjbase_datasets
from api.genomic.genomic import get_genomic_species, get_genomic_datasets
from db.vdjbase_model import Gene as VDJbaseGene, Allele as VDJbaseAllele, AllelesSample as VDJbaseAllelesSample
from db.genomic_db import Gene as GenomicGene, Sequence as GenomicSequence
from db.vdjbase_airr_model import Sample as VDJbaseSample

ns = api.namespace('refbook', description='Refbook related operations')


def chains_in_loc(loc):
    ret = [loc + 'V', loc + 'J']
    if loc in ['IGH', 'TRB', 'TRA']:
        ret.append(loc + 'D')
    return ret


@ns.route('/species_and_chains')
@api.response(404, 'No species available!')
class SpeciesApi(Resource):
    @digby_protected()
    def get(self):
        """ Returns the list of species and chains for which information is held """

        ret = {'species': [], 'chains': {}}

        genomic_sp = get_genomic_species()

        for sp in genomic_sp:
            ret['species'].append(sp)
            for rec in get_genomic_datasets(sp):
                loc = rec['dataset']
                if 'IGHC' not in loc:
                    if sp not in ret['chains']:
                        ret['chains'][sp] = []
                    ret['chains'][sp].extend(chains_in_loc(loc))

        vdjbase_sp = get_vdjbase_species()

        for sp in vdjbase_sp:
            if sp not in ret['species']:
                ret['species'].append(sp)
            for rec in get_vdjbase_datasets(sp):
                loc = rec['dataset']
                if 'IGHC' not in loc:
                    if sp not in ret['chains']:
                        ret['chains'][sp] = []
                    for ch in chains_in_loc(loc):
                        if ch not in ret['chains'][sp]:
                            ret['chains'][sp].append(ch)

        return ret


@ns.route('/ascs_in_chain/<string:species>/<string:chain>')
@api.response(404, 'Species or chain not found')
class AscsInChainApi(Resource):
    @digby_protected()
    def get(self, species, chain):
        """ Returns the list of ASCs in a given chain for a given species """

        species_chains = SpeciesApi.get(self)
        
        if species not in species_chains['species']:
            return {'message': 'Species not found'}, 404
        if chain not in species_chains['chains'][species]:
            return {'message': 'Chain not found'}, 404

        ret = []
        ds = chain[:3]

        if ds in vdjbase_dbs[species]:
            session = vdjbase_dbs[species][ds].session
            genes = session.query(VDJbaseGene.name).filter(VDJbaseGene.type == chain).all()
            ret.extend([g[0] for g in genes])

        if ds in genomic_dbs[species]:
            session = genomic_dbs[species][ds].session
            genes = session.query(GenomicGene.name).filter(GenomicGene.type == chain).all()
            ret.extend(genes)

        ret = sorted(set([g[0] for g in genes]))
        return ret


@ns.route('/ascs_overview/<string:species>/<string:chain>/<string:asc>')
@api.response(404, 'Species or chain not found')
class AscsOverview(Resource):
    @digby_protected()
    def get(self, species, chain, asc):
        """ Returns data for the overview refbook component """

        species_chains = SpeciesApi.get(self)
        
        if species not in species_chains['species']:
            return {'message': 'Species not found'}, 404
        if chain not in species_chains['chains'][species]:
            return {'message': 'Chain not found'}, 404
        
        ret = {}
        alleles = {}

        ds = chain[:3]
        if ds in vdjbase_dbs[species]:
            session = vdjbase_dbs[species][ds].session
            vdjbase_alleles = session.query(VDJbaseAllele.name, VDJbaseAllele.novel, VDJbaseAllele.appears) \
                .join(VDJbaseGene) \
                .filter(VDJbaseGene.type == chain, VDJbaseGene.name == asc) \
                .all()
            
            for a, novel, appearances in vdjbase_alleles:
                alleles[a] = {'VDJbase': appearances, 'Genomic': 0, 'novel': novel}
        
        if ds in genomic_dbs[species]:
            session = genomic_dbs[species][ds].session
            genomic_alleles = session.query(GenomicSequence.name, GenomicSequence.novel, GenomicSequence.appearances) \
                .join(GenomicGene) \
                .filter(GenomicGene.type == chain, GenomicGene.name == asc, or_(GenomicSequence.functional == 'Functional', GenomicSequence.functional == 'ORF')) \
                .all()
            
            for a, novel, appearances in genomic_alleles:
                if a not in alleles:
                    alleles[a] = {'VDJbase': 0, 'Genomic': appearances, 'novel': novel}
                else:
                    alleles[a]['Genomic'] = appearances

        ret['total'] = len(alleles)
        ret['novel'] = sum(rec['novel'] for rec in alleles.values())
        ret['baseline'] = ret['total'] - ret['novel']
        alleles = dict(sorted(alleles.items()))
        ret['alleles'] = list(alleles.keys())
        # these counts may not be exactly what we want, I am not sure what to do if there are samples
        # for which we don't have both genomic and airr-seq results
        ret['genomic_only_counts'] = [rec['Genomic'] if rec['VDJbase'] == 0 else 0 for rec in alleles.values()]
        ret['vdjbase_only_counts'] = [rec['VDJbase'] if rec['Genomic'] == 0 else 0 for rec in alleles.values()]
        ret['both_counts'] = [min(rec['Genomic'], rec['VDJbase']) if rec['Genomic'] > 0 and rec['VDJbase'] > 0 else 0 for rec in alleles.values()]

        return ret


@ns.route('/asc_seqs/<string:species>/<string:chain>/<string:asc>')
@api.response(404, 'Species or chain not found')
class AscSeqs(Resource):
    @digby_protected()
    def get(self, species, chain, asc):
        """ Returns sequences of all alleles in an ASC """

        species_chains = SpeciesApi.get(self)
        
        if species not in species_chains['species']:
            return {'message': 'Species not found'}, 404
        if chain not in species_chains['chains'][species]:
            return {'message': 'Chain not found'}, 404
        
        ret = {}
        alleles = []
        recs = []

        ds = chain[:3]
        if ds in vdjbase_dbs[species]:
            session = vdjbase_dbs[species][ds].session
            vdjbase_alleles = session.query(VDJbaseAllele.name, VDJbaseAllele.seq) \
                .join(VDJbaseGene) \
                .filter(VDJbaseGene.type == chain, VDJbaseGene.name == asc, 'Del' not in VDJbaseAllele.name) \
                .all()
            
            for a, seq_gapped in vdjbase_alleles:
                recs.append({'name': a, 'seq_gapped': seq_gapped.upper(), 'seq': seq_gapped.upper().replace('.', '')})
                alleles.append(a)
        
        if ds in genomic_dbs[species]:
            session = genomic_dbs[species][ds].session
            genomic_alleles = session.query(GenomicSequence.name, GenomicSequence.gapped_sequence, GenomicSequence.sequence) \
                .join(GenomicGene) \
                .filter(GenomicGene.type == chain, GenomicGene.name == asc, or_(GenomicSequence.functional == 'Functional', GenomicSequence.functional == 'ORF')) \
                .all()
            
            for a, gapped, ungapped in genomic_alleles:
                if a not in alleles:
                    recs.append({'name': a, 'seq_gapped': gapped.upper(), 'seq': ungapped.upper()})

        return {'alleles': recs}

@ns.route('/asc_usage/<string:species>/<string:chain>/<string:asc>')
@api.response(404, 'Species or chain not found')
class AscUsage(Resource):
    @digby_protected()
    def get(self, species, chain, asc):
        """ Returns usage statistics for all alleles in an ASC """

        species_chains = SpeciesApi.get(self)
        
        if species not in species_chains['species']:
            return {'message': 'Species not found'}, 404
        if chain not in species_chains['chains'][species]:
            return {'message': 'Chain not found'}, 404
        
        
        ds = chain[:3]

        if ds in vdjbase_dbs[species]:
            session = vdjbase_dbs[species][ds].session
            
            totals = (
                session.query(
                    VDJbaseAllelesSample.patient_id,
                    VDJbaseAllelesSample.sample_id,
                    func.sum(VDJbaseAllelesSample.total_count).label("total")
                )
                .group_by(VDJbaseAllelesSample.patient_id, VDJbaseAllelesSample.sample_id)
                .subquery()
            )
            
            fraction_expr = (
                cast(VDJbaseAllelesSample.count, Float) /
                func.nullif(cast(totals.c.total, Float), 0.0)
            )

            vdjbase_alleles = (
                session.query(
                    VDJbaseAllele.name,
                    func.group_concat(VDJbaseSample.sample_name).label("samples"),
                    func.group_concat(func.coalesce(fraction_expr, 0.0)).label("fraction"),
                )
                .join(VDJbaseAllele, VDJbaseAllele.id == VDJbaseAllelesSample.allele_id)
                .join(VDJbaseSample, VDJbaseSample.id == VDJbaseAllelesSample.sample_id)
                .join(
                    totals,
                    (VDJbaseAllelesSample.patient_id == totals.c.patient_id)
                    & (VDJbaseAllelesSample.sample_id == totals.c.sample_id),
                )
                .join(VDJbaseGene, VDJbaseGene.id == VDJbaseAllele.gene_id)
                .group_by(VDJbaseAllele.name)
                .filter(
                    VDJbaseGene.type == chain,
                    VDJbaseGene.name == asc,
                    VDJbaseAllelesSample.count.isnot(None)
                )
                .all()
            )

        recs = [{'name': allele, 'usage': list(usages.split(',') if usages else []), 'samples': list(samples.split(',') if samples else [])} for allele, samples, usages in vdjbase_alleles]

        return {'alleles': recs}

@ns.route('/asc_zygousity/<string:species>/<string:chain>/<string:asc>')
@api.response(404, 'Species or chain not found')
class AscZygosity(Resource):
    @digby_protected()
    def get(self, species, chain, asc):
        """ Returns zygosity statistics for all subjects in a given ASC """

        species_chains = SpeciesApi.get(self)
        
        if species not in species_chains['species']:
            return {'message': 'Species not found'}, 404
        if chain not in species_chains['chains'][species]:
            return {'message': 'Chain not found'}, 404
        
        recs = []

        ds = chain[:3]

        if ds in vdjbase_dbs[species]:
            session = vdjbase_dbs[species][ds].session

            alleles_per_sample = (
                session.query(
                    VDJbaseAllelesSample.sample_id,
                    #VDJbaseSample.sample_name,
                    func.group_concat(func.distinct(VDJbaseAllele.name)).label("alleles"),
                )
                .join(VDJbaseAllele, VDJbaseAllele.id == VDJbaseAllelesSample.allele_id)
                .join(VDJbaseGene, VDJbaseGene.id == VDJbaseAllele.gene_id)
                #.join(VDJbaseSample, VDJbaseSample.sample_id == VDJbaseAllelesSample.sample_id)
                .filter(
                    VDJbaseGene.type == chain,
                    VDJbaseGene.name == asc,
                )
                .group_by(VDJbaseAllelesSample.sample_id)
                .all()
            )

            
            for sample_id, alleles in alleles_per_sample:
                allele_list = alleles.split(',') if alleles else []
                recs.append({
                    "name": sample_id,
                    "sets": list(allele_list)
                })

        return {'samples': recs}
