# Manage a list of available vdjbase-style databases

from os.path import join, isdir, isfile
from os import listdir

from sqlalchemy import create_engine, inspect, text
from sqlalchemy.orm import sessionmaker

from db.vdjbase_maint import create_single_database
from db.vdjbase_model import Details



Session = sessionmaker()


class ContentProvider():
    db = None
    connection = None
    session = None
    description = None
    binomial = None
    taxid = None
    created = None

    def __init__(self, path, description):
        self.db = create_engine('sqlite:///' + path + '?check_same_thread=false', echo=False)
        self.connection = self.db.connect()
        self.session = Session(bind=self.connection)
        self.description = description

    def get_session(self):
        if self.session is None or not self.session.is_active:
            print('renewing session')
            self.session = Session(bind=self.connection)
        return self.session
    
    def close_session(self):
        if self.session is not None and self.session.is_active:
            self.session.close()
            self.session = None

    def close(self):
        self.connection.close()
        self.db.dispose()


# TODO: when we change to biomial species names, update this lookup and corresponding code
species_lookup = {
    'Human': ('Homo sapiens', 'NCBITAXON:9606'),
    'Rhesus Macaque': ('Macaca mulatta', 'NCBITAXON:9544'),
    'Crab-eating Macaque': ('Macaca fascicularis', 'NCBITAXON:9541'),
}


def study_data_db_init(vdjbase_db_path):
    sqlite_dbs = {}

    for species in listdir(vdjbase_db_path):
        p = join(vdjbase_db_path, species)
        if isdir(p) and species[0] != '.':
            for name in listdir(p):
                if isdir(join(p, name)) and name[0] != '.' and '.txt' not in name:
                    description = ''
                    if isfile(join(p, name, 'db_description.txt')):
                        with open(join(p, name, 'db_description.txt'), 'r') as fi:
                            description = ' '.join(fi.readlines())
                    if species not in sqlite_dbs:
                        sqlite_dbs[species] = {}
                    sqlite_dbs[species][name] = ContentProvider(join(p, name, 'db.sqlite3'), description)

                    session = sqlite_dbs[species][name].session
                    try:
                        sqlite_dbs[species][name].created = session.query(Details.created_on).one_or_none()[0]
                    except Exception as e:
                        print('Error querying Details table for %s/%s: %s' % (species, name, e))

                    if species in species_lookup:
                        sqlite_dbs[species][name].binomial = species_lookup[species][0]
                        sqlite_dbs[species][name].taxid = species_lookup[species][1]
                    else:
                        sqlite_dbs[species][name].binomial = species
                        sqlite_dbs[species][name].taxid = ''
                        print('Species %s not found in species lookup' % species)

    # temp fix: add asc_genotype column to sample table if not there already

    for species in sqlite_dbs:
        for locus in sqlite_dbs[species]:
            inspector = inspect(sqlite_dbs[species][locus].db)

            # Was gated on the path not containing "genomic", which assumed every
            # other kind of database has a Sample table. Ask the database instead,
            # so a flavour without one is simply skipped.
            if 'sample' not in {name.lower() for name in inspector.get_table_names()}:
                continue

            cols = inspector.get_columns('Sample')
            if 'asc_genotype' not in [col['name'] for col in cols]:
                # Not `with ... as con`: exiting that block closes the connection,
                # and this dataset's long-lived session is bound to it, so every
                # later query raised ResourceClosedError. The branch only runs when
                # the column is missing, so it broke the first start after a
                # database was rebuilt and looked fine on the next one.
                sqlite_dbs[species][locus].connection.execute(
                    'ALTER TABLE Sample ADD COLUMN asc_genotype text')
                sqlite_dbs[species][locus].session.commit()


    # sort datasets of each species

    for species in sqlite_dbs:
        sqlite_dbs[species] = dict(sorted(sqlite_dbs[species].items(), key=lambda kv: kv[0]))

    # put Human at the front
    sqlite_dbs = dict(sorted(sqlite_dbs.items(), key=lambda kv: 'aaaaa' if kv[0] == 'Human' else kv[0]))

    return sqlite_dbs


