#
# backend maintenance functions for VDJbase-style databases
#
import os.path
import os
import sys
import traceback
import zipfile
import datetime

from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
from sqlalchemy.orm import sessionmaker


from db.vdjbase_confidence import check_novel_confidence
from db.vdjbase_cross_sample import update_alleles_appearance, calculate_gene_frequencies, calculate_patterns
from db.vdjbase_model import *
from db.vdjbase_reference import import_reference_alleles
from db.vdjbase_projects import import_studies
from db.vdjbase_genotypes import process_genotypes, add_deleted_alleles, process_haplotypes_and_stats
from db.vdjbase_exceptions import *
from db.source_details import db_source_details


Session = sessionmaker()

def create_single_database(job, species, dataset, upload_path, in_situ=False):
    result = []

    # if called from UI, create/unpack in uploads directory
    # if called from command line, create in situ from files without unpacking

    if not in_situ:
        ds_dir = os.path.join(upload_path, species, dataset)
        try:
            job.update_state(state='PENDING', meta={'value': 'Unzipping data files'})
            result.extend(extract_files(job, ds_dir, species, dataset))
        except Exception as e:
            result.append(e.args[0])
            traceback.print_exc(limit=2, file=sys.stdout)
            return False, result
    else:
        ds_dir = upload_path

    db_file = os.path.join(ds_dir, 'db.sqlite3')

    if os.path.isfile(db_file):
        os.remove(db_file)

    engine = create_engine('sqlite:///' + db_file, echo=False, poolclass=NullPool)
    Base.metadata.create_all(engine)
    db_connection = engine.connect()
    engine.session = Session(bind=db_connection)
    session = engine.session

    commit_id, branch = db_source_details()
    details = session.query(Details).one_or_none()
    if not details:
        details = Details(
            dbtype='airrseq',
            species=species, 
            locus=dataset, 
            created_on=datetime.datetime.now(),
            created_by='digby_backend',
            software_commit_id=commit_id,
            software_branch=branch
            )
        session.add(details)

    success = True

    try:
        job.update_state(state='PENDING', meta={'value': 'Importing reference alleles'})
        result.extend(import_reference_alleles(os.path.join(ds_dir, 'reference'), session, species))
        job.update_state(state='PENDING', meta={'value': 'Importing studies'})
        result.extend(import_studies(ds_dir, species, dataset, session))
        job.update_state(state='PENDING', meta={'value': 'Processing deleted alleles'})
        result.extend(add_deleted_alleles(session))
        job.update_state(state='PENDING', meta={'value': 'Processing genotypes'})
        result.extend(process_genotypes(ds_dir, species, dataset, session))
        job.update_state(state='PENDING', meta={'value': 'Processing haplotypes'})
        result.extend(process_haplotypes_and_stats(ds_dir, species, dataset, session))
        job.update_state(state='PENDING', meta={'value': 'Analyzing allele appearances'})
        result.extend(update_alleles_appearance(session))
        job.update_state(state='PENDING', meta={'value': 'Calculating gene frequencies'})
        result.extend(calculate_gene_frequencies(ds_dir, session))
        job.update_state(state='PENDING', meta={'value': 'Calculating patterns'})
        result.extend(calculate_patterns(session))

        job.update_state(state='PENDING', meta={'value': 'Creating confidence reports'})
        result.extend(check_novel_confidence(ds_dir, session))
    except Exception as e:
        result.append(e.args[0])
        traceback.print_exc(limit=2, file=sys.stdout)
        success = False
    finally:
        db_connection.close()
        engine.dispose()

    return success, result


def extract_files(job, ds_dir, species, dataset):
    result = []

    with zipfile.ZipFile(os.path.join(ds_dir, 'data.zip')) as zip_ref:
        zip_ref.extractall(ds_dir)

    for d in ['reference', 'samples']:
        if not os.path.isdir(os.path.join(ds_dir, d)):
            result.append('datafile directory %s not found' % d)
            raise DbCreationError('datafile directory %s not found')

    result.append('Processing %s/%s' %(species, dataset))

    return result

