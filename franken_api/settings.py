import os

# Flask settings
FLASK_DEBUG = True  # Do not use debug mode in production

# Flask-Restplus settings
# RESTX_SWAGGER_UI_DOC_EXPANSION = 'list'
RESTX_VALIDATE = True
RESTX_MASK_SWAGGER = False
BUNDLE_ERRORS = False

# SQLAlchemy settings
SQLALCHEMY_DATABASE_URI = os.environ['DATABASE_URL']  #'mysql+pymysql://username@host/database_name'
SQLALCHEMY_BINDS = { "curation": os.environ['CURATION_DB_URL'], "leaderboard" : os.environ['LEADERBOARD_DB_URL'] , "ipcmLeaderboard": os.environ['IPCM_DB_URL']}
SQLALCHEMY_TRACK_MODIFICATIONS = False