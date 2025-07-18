# Franken Plot  API

### Setup and Prerequisite installation

#### Back-End
    cd franken_api
```
    python setup.py  install
    pip install -r requirements.txt
    pip install .
``` 
### Run application
```
    franken_api -p 5000
```
### Run DB Migration script:

```
sudo docker exec -it probio_ui_docker_postgress_1 bash
su postgres
psql
CREATE USER referral_writer WITH ENCRYPTED PASSWORD '<password: this should match with password in flask app settings.py>';
CREATE DATABASE referrals;
CREATE DATABASE curation;
GRANT ALL PRIVILEGES ON DATABASE referrals TO referral_writer;
GRANT ALL PRIVILEGES ON DATABASE curation TO referral_writer;
exit

python franken_api/migrate.py db init --multidb
python franken_api/migrate.py db migrate
python franken_api/migrate.py db upgrade
```
### Generate Profile Summary for Project 
```
    cd franken_api/genomic_profile_lookup/
    python read_curation_info_store_postgresql_final.py '<nfs-autoseq-path>' '<project-name>
```

### Export MTBP json format 
    --
```
    export MTBP_SCRIPT_PATH="<MTBP_Script_path>"

```