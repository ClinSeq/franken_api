import logging.config
import click
from flask import Flask, Blueprint

from franken_api import settings
from franken_api.api.franken.endpoints.franken_api import ns as franken_namespace
from franken_api.api.franken.endpoints.referral_db_api import ns2 as referral_namespace
from franken_api.api.franken.endpoints.curation_db_api import ns3 as curation_namespace
from franken_api.api.restplus import api
from franken_api.database import db

def configure_app(flask_app):
    flask_app.config['SQLALCHEMY_DATABASE_URI'] = settings.SQLALCHEMY_DATABASE_URI
    flask_app.config['SQLALCHEMY_BINDS'] = {'curation': settings.SQLALCHEMY_BINDS}
    flask_app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = settings.SQLALCHEMY_TRACK_MODIFICATIONS
    flask_app.config['SWAGGER_UI_DOC_EXPANSION'] = settings.RESTPLUS_SWAGGER_UI_DOC_EXPANSION
    flask_app.config['RESTPLUS_VALIDATE'] = settings.RESTPLUS_VALIDATE
    flask_app.config['RESTPLUS_MASK_SWAGGER'] = settings.RESTPLUS_MASK_SWAGGER
    flask_app.config['ERROR_404_HELP'] = settings.RESTPLUS_ERROR_404_HELP
    flask_app.config['PROBIO'] = settings.MOUNT_POINT_PROBIO
    flask_app.config['PSFF'] = settings.MOUNT_POINT_PSFF
    flask_app.config['LARS_DUCTAL'] = settings.MOUNT_POINT_LARS_DUCTAL
    flask_app.config['AZ_RINGTRIAL'] = settings.MOUNT_POINT_AZ_RINGTRIAL
    flask_app.config['msHSPC'] = settings.MOUNT_POINT_MS_HSPC
    flask_app.config['RESBIO'] = settings.MOUNT_POINT_RESBIO
    flask_app.config['HD_C3'] = settings.MOUNT_POINT_HD_C3
    flask_app.config['LPC'] = settings.MOUNT_POINT_LPC
    flask_app.config['CHEERS'] = settings.MOUNT_POINT_CHEERS



def initialize_app(flask_app):
    configure_app(flask_app)
    blueprint = Blueprint('api', __name__, url_prefix='/api')
    api.init_app(blueprint)
    api.add_namespace(franken_namespace)
    api.add_namespace(referral_namespace)
    api.add_namespace(curation_namespace)
    flask_app.register_blueprint(blueprint)
    db.init_app(flask_app)


app = Flask(__name__)

@click.command()
@click.option('-p', '--port', type=int, default=5000)
def main(port):
    initialize_app(app)
    logging.info('>>>>> Starting development server at http://{}/api/ <<<<<'.format('localhost'))
    app.run(debug=settings.FLASK_DEBUG, port=port, host='0.0.0.0')
