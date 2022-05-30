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
SQLALCHEMY_BINDS = { "curation": os.environ['CURATION_DB_URL'], "leaderboard" : os.environ['LEADERBOARD_DB_URL'] , "ipcmLeaderboard": os.environ['IPCM_DB_URL']}

SQLALCHEMY_TRACK_MODIFICATIONS = False

# MTBP Python script path 
MTBP_SCRIPT_PATH = os.environ['MTBP_SCRIPT_PATH']
PDF_SCRIPT_PATH = os.environ['PDF_SCRIPT_PATH']


#path to franken json files
MOUNT_POINT_PROBIO = '/nfs/PROBIO/autoseq-output'
MOUNT_POINT_PSFF = '/nfs/PSFF/autoseq-output'
MOUNT_POINT_LARS_DUCTAL = '/nfs/CLINSEQ/LARS_DUCTAL/autoseq-output'
MOUNT_POINT_AZ_RINGTRIAL = '/nfs/CLINSEQ/AZ_ringtrial/autoseq-output'
MOUNT_POINT_MS_HSPC = '/nfs/PROBIO/autoseq-output/msHSPC'
MOUNT_POINT_RESBIO = '/nfs/RESBIO/autoseq-output'
MOUNT_POINT_HD_C3 = '/nfs/PROBIO2/healthy_donors_comprehensive/comprehensive_v3/autoseq-output'
MOUNT_POINT_LPC = '/nfs/PROBIO3/LPC/autoseq-output'
MOUNT_POINT_CHEERS = '/nfs/CLINSEQ/JanPieter/CHEERS/autoseq-output'
MOUNT_POINT_ULLEN = '/nfs/ULLEN/autoseq-output'
MOUNT_POINT_CRC_REFLEX = '/nfs/CLINSEQ/CRC_REFLEX/autoseq-output'
MOUNT_POINT_IPCM = '/nfs/IPCM/autoseq-output'
MOUNT_POINT_EVAL_SAMPLES = '/nfs/PROBIO3/EVAL_SAMPLES/autoseq-output'
