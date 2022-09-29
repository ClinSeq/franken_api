import logging.config
import click
from flask import Flask, Blueprint, request, redirect

from franken_api import settings
from franken_api.api.franken.endpoints.authentication_db_api import au as auth_namespace
from franken_api.api.franken.endpoints.franken_api import ns as franken_namespace
from franken_api.api.franken.endpoints.referral_db_api import ns2 as referral_namespace
from franken_api.api.franken.endpoints.curation_db_api import ns3 as curation_namespace
from franken_api.api.restplus import api
from franken_api.database import db

def configure_app(flask_app):
    flask_app.config['SQLALCHEMY_DATABASE_URI'] = settings.SQLALCHEMY_DATABASE_URI
    flask_app.config['SQLALCHEMY_BINDS'] = settings.SQLALCHEMY_BINDS
    flask_app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = settings.SQLALCHEMY_TRACK_MODIFICATIONS
    # flask_app.config['SWAGGER_UI_DOC_EXPANSION'] = settings.RESTX_SWAGGER_UI_DOC_EXPANSION
    flask_app.config['RESTX_VALIDATE'] = settings.RESTX_VALIDATE
    flask_app.config['RESTX_MASK_SWAGGER'] = settings.RESTX_MASK_SWAGGER
    flask_app.config['BUNDLE_ERRORS'] = settings.BUNDLE_ERRORS
    flask_app.config['MTBP_SCRIPT'] = settings.MTBP_SCRIPT_PATH
    flask_app.config['PDF_SCRIPT'] = settings.PDF_SCRIPT_PATH

def initialize_app(flask_app):
    configure_app(flask_app)
    blueprint = Blueprint('api', __name__, url_prefix='/api')
    api.init_app(blueprint)
    api.add_namespace(auth_namespace)
    api.add_namespace(franken_namespace)
    api.add_namespace(referral_namespace)
    api.add_namespace(curation_namespace)
    flask_app.register_blueprint(blueprint)
    db.init_app(flask_app)


app = Flask(__name__)

'''
@app.before_request
def before_request():

    if request.scheme == 'https' and request.is_secure:
        url = request.url.replace("http://", "https://", 1)
        code = 301
        return redirect(url, code=code)
'''

@click.command()
@click.option('-p', '--port', type=int, default=5000)
def main(port):
    initialize_app(app)
    logging.info('>>>>> Starting development server at http://{}/api/ <<<<<'.format('localhost'))
    app.run(debug=settings.FLASK_DEBUG, port=port, host='0.0.0.0')

if __name__ == '__main__':
    app.run(debug=True)
