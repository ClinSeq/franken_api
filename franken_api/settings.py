import os

# Flask settings
FLASK_DEBUG = True  # Do not use debug mode in production

# Flask-Restplus settings
RESTPLUS_SWAGGER_UI_DOC_EXPANSION = 'list'
RESTPLUS_VALIDATE = True
RESTPLUS_MASK_SWAGGER = False
RESTPLUS_ERROR_404_HELP = False

# SQLAlchemy settings
SQLALCHEMY_DATABASE_URI = os.environ['DATABASE_URL']  #'mysql+pymysql://username@host/database_name'
SQLALCHEMY_BINDS = os.environ['CURATION_DB_URL']

SQLALCHEMY_TRACK_MODIFICATIONS = False


#path to franken json files
MOUNT_POINT_PROBIO = '/nfs/PROBIO/autoseq-output'
# MOUNT_POINT_PROBIO = '/home/karman/probio/test/PROBIO/autoseq-output'
MOUNT_POINT_PSFF = '/nfs/PSFF/autoseq-output'
# MOUNT_POINT_PSFF = '/nfs/CLINSEQ/PSFF/autoseq-output'
MOUNT_POINT_LARS_DUCTAL= '/nfs/CLINSEQ/LARS_DUCTAL/autoseq-output'