import logging.config
import click
from flask import Flask, Blueprint

from franken_api import settings
from franken_api.api.franken.endpoints.authentication_db_api import au as auth_namespace
from franken_api.api.franken.endpoints.franken_api import ns as franken_namespace
from franken_api.api.franken.endpoints.referral_db_api import ns2 as referral_namespace
from franken_api.api.franken.endpoints.curation_db_api import ns3 as curation_namespace
from franken_api.api.restplus import api
from franken_api.database import db
from flask_migrate import Migrate, MigrateCommand
from flask_script import Manager

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
initialize_app(app)
migrate = Migrate(app, db)
manager = Manager(app)

manager.add_command('db', MigrateCommand)

if __name__ == '__main__':
    manager.run()
