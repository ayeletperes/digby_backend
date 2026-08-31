import os
import custom_logging
import yaml

from werkzeug.middleware.proxy_fix import ProxyFix
from flask import Flask, Blueprint, send_from_directory
from flask_migrate import Migrate
from flask_mail import Mail
from flask_bootstrap import Bootstrap
from flask_cors import CORS

from flask_swagger_ui import get_swaggerui_blueprint
from extensions import celery

from db.vdjbase_db import study_data_db_init
from db.madc import madc_init

sql_db = None

app = Flask(__name__, static_folder=None)
app.static_url_path = '/static'
bootstrap = Bootstrap(app)
app.static_url_path = None
app.config.from_pyfile('config.cfg')
app.config.from_pyfile('secret.cfg')

# configure/initialize all your extensions

sql_db = celery.init_app(app)

app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1, x_prefix=1)

CORS(app)

app.config['BASE_PATH'] = os.getcwd()

if 'STATIC_PATH' not in app.config:
    app.config['STATIC_PATH'] = os.path.join(app.config['BASE_PATH'], 'static')

if 'OUTPUT_PATH' not in app.config:
    app.config['OUTPUT_PATH'] = os.path.join(app.config['STATIC_PATH'], 'output')

if 'UPLOAD_PATH' not in app.config:
    app.config['UPLOAD_PATH'] = os.path.join(app.config['BASE_PATH'], 'uploads')

app.config['R_SCRIPT_PATH'] = os.path.join(app.config['BASE_PATH'], 'api/reports/R_scripts')

if 'R_LIBS' not in os.environ or os.environ['R_LIBS'] is None or len(os.environ['R_LIBS']) < 1:
    os.environ['R_LIBS'] = app.config['R_SCRIPT_PATH']

# TODO - make this work on Windows as well as Unix, the D:\ screws it up
"""
else:
    r_libs = os.environ['R_LIBS'].split(':')
    r_libs = r_libs.append(app.config['R_SCRIPT_PATH'])
    os.environ['R_LIBS'] = ':'.join(r_libs)
"""


mail = Mail(app)
custom_logging.init_logging(app, mail)

vdjbase_dbs = study_data_db_init(os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/db'))
genomic_dbs = study_data_db_init(os.path.join(app.config['STATIC_PATH'], 'study_data/Genomic/db'))
qtl_dbs = study_data_db_init(os.path.join(app.config['STATIC_PATH'], 'study_data/QTL/db')) \
    if os.path.isdir(os.path.join(app.config['STATIC_PATH'], 'study_data/QTL/db')) else {}
madc_index = madc_init(app)

from api.restx import api
from api.genomic.genomic import ns as genomic
from api.vdjbase.vdjbase import ns as vdjbase
from api.reports.reports import ns as reports
from api.refbook.refbook import ns as refbook
from api.refbook.sunburst import ns as refbook_sunburst
from api.qtl.qtl import ns as qtl
# tree and region build their namespace with api.namespace(), which registers it
# on import, so importing them is all the wiring they need. Bind them under a
# different name: `import api.refbook.tree` would rebind `api` here to the
# package, shadowing the Api object imported from api.restx.
from api.refbook import tree as _refbook_tree      # noqa: F401
from api.qtl import region as _qtl_region          # noqa: F401
from api.system.system import ns as system, digby_protected

from db.genomic_db import *
import db.vdjbase_maint
import db.vdjbase_export
from db.vdjbase_igsnper import do_igsnper

migrate = Migrate(app, sql_db)

blueprint = Blueprint('api', __name__, url_prefix='/api')
api.init_app(blueprint)
api.add_namespace(genomic)
api.add_namespace(vdjbase)
api.add_namespace(reports)
api.add_namespace(system)
api.add_namespace(refbook)
api.add_namespace(refbook_sunburst)
api.add_namespace(qtl)
app.register_blueprint(blueprint)

from api_v1.open_api import api_bp

# PROXY_PREFIX is needed to get around a bug in swagger_ui: https://github.com/PWZER/swagger-ui-py/issues/34
# When this is fixed, ProxyFix will take care of things and the special location for /admin/api_v1 in nginix confs can be removed

app.register_blueprint(api_bp, url_prefix="/api/v1")

SWAGGER_URL = "/api/v1"
API_URL = '/static/vdjbase-api-openapi3.yaml'
swaggerui_blueprint = get_swaggerui_blueprint(
    SWAGGER_URL,
    API_URL,
    config={
        'app_name': "VDJbase Open API v1"
    }
)


@app.route('/schema/<path:filename>')
def serve_schema(filename):
    return send_from_directory('schema', filename)


app.register_blueprint(swaggerui_blueprint, url_prefix=SWAGGER_URL)

with open('schema/vdjbase-api-openapi3.yaml', 'r') as f:
    openapi_schema = yaml.safe_load(f)

from api.reports.reports import load_report_defs
load_report_defs()

from flask_jwt_extended import JWTManager

app.config["JWT_TOKEN_LOCATION"] = ["headers", "cookies"]
jwt = JWTManager(app)

@app.route('/static/<path:path>', methods=['GET', 'POST'])
def static(path):
    if '/gff' in path:
        return send_from_gff(path)
    else:
        return send_from_directory(app.config['STATIC_PATH'], path)

@digby_protected()
def send_from_gff(path):
    return send_from_directory(app.config['STATIC_PATH'], path)

