#!/usr/bin/env python
# coding: utf-8

# ### File name format
#    **MTBP_iPCM_KI_3405_CFDNA03884517**
#
#     *MTBP - molecular tumor board portal
#     iPCM - clinical trial/study/project
#     KI - Karolinska Hospital (site) (other sites ST - St Goran, SOS - Sodersjuksjuset)
#     3405 - Patient Identifier
#     CFDNA03884517 - sample type with ID number*

## Usage : python3 MTBP_json_format.py --path <nfs-output-path> --project_name <project-name>  --capture_format <project-prefix-code> --sample_ids <sample-ids>
## Eg: 
# export user_name='<MTBP-UserName>'
# export user_pwd='<MTBP-Password>' 
# python3 MTBP_json_format.py --path /data/demo/PROBIO/autoseq-output --project_name PROBIO  --capture_format PB --sample_ids P-00477166 P-586956 P-00466944


import os
import pandas as pd
import numpy as np
import json
import re
import logging
from logging.handlers import QueueHandler
import psycopg2
import psycopg2.extras
from sqlalchemy import create_engine
from configparser import ConfigParser
import yaml
import argparse
import sys
from datetime import date, datetime
import math
from decimal import Decimal
import subprocess
from subprocess import Popen
import csv

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')

logger1 = logging.getLogger('MTBP Overall')
logger1.setLevel(logging.INFO)

logger2 = logging.getLogger('MTBP Samplewise')
logger2.setLevel(logging.INFO)


# ### Get the current directory path
def path():
	return os.path.dirname(__file__)

# ### Read a DB information from yml file
def readConfig(section,filename=path()+"/config.yml"):
	with open(filename, "r") as ymlfile:
		cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
		section = cfg[section]
	return section


# ### Fetch SQL Query
def fetch_sql_query(section, sql):
	conn = None
	try:
		params = readConfig(section)
		conn = psycopg2.connect(**params)
		cur = conn.cursor(cursor_factory = psycopg2.extras.RealDictCursor) # RealDictCursor, NamedTupleCursor
		cur.execute(sql)
		rows = cur.fetchall()
		data = [dict(r) for r in rows]
		cur.close()
		return data
	except (Exception, psycopg2.DatabaseError) as error:
		print("SQL Fetch Exception", str(error))
	finally:
		if conn is not None:
			conn.close()

def json_serial(obj):
	if isinstance(obj, (datetime, date)):
		return obj.isoformat()
	elif isinstance(obj, Decimal):
		return float(obj)
	elif isinstance(obj, str) and obj.isdigit():
		return int(obj)
	raise TypeError("Type %s not serializable" % type(obj))

def convert_to_numeric(string):
	try:
		return int(string)
	except ValueError:
		try:
			return float(string)
		except ValueError:
			return 'NA'

def validate_tissue_type(ecrf_tissue_type, sample_type):
	tissue_type = 'NA'
	if ecrf_tissue_type == 'Cytology':
		tissue_type = ecrf_tissue_type.lower()
	elif  ecrf_tissue_type == 'FFPE|cfDNA':
		tissue_type = 'FFPE' if sample_type == 'T' else 'cfdna'
	elif  ecrf_tissue_type == 'FFPE|Cytology':
		tissue_type = 'FFPE' if sample_type == 'T' else 'cfdna'
	elif  ecrf_tissue_type == 'Cytology|cfDNA':
		tissue_type = 'cytology' if sample_type == 'T' else 'cfdna'
	else:
		tissue_type = ecrf_tissue_type
	return tissue_type

# ### Fetch the sample information from ipcm referral table
def build_ipcm_sample_details(normal_cfdna, cfdna, capture_format, sample_type, seq_date, germline_dna):

	identifier_status = False
	res_json_data = []
	try:

		## tissue cfdna
		cfdna = re.sub(r'[a-zA-Z]', '', cfdna)
		sql = "select rf.pnr, rf.cdk from ipcm_referral_t as rf WHERE rf.rid like'%{}%' OR rf.blood like'%{}%' OR rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' OR rf.dna3 like'%{}%'".format(cfdna, cfdna, cfdna, cfdna, cfdna)
		res_data = fetch_sql_query('ipcmLeaderboard', sql)

		## normal cfdna
		sql2 = "select rf.pnr, rf.cdk from ipcm_referral_t as rf WHERE rf.rid like'%{}%' OR rf.blood like'%{}%' OR rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' OR rf.dna3 like'%{}%'".format(normal_cfdna, normal_cfdna, normal_cfdna, normal_cfdna, normal_cfdna)
		res_normal_data = fetch_sql_query('ipcmLeaderboard', sql2)

		t_pnr = res_data[0]['pnr'] if len(res_data) > 0 else ''
		t_cdk = res_data[0]['cdk'] if len(res_data) > 0 else ''

		n_pnr = res_normal_data[0]['pnr'] if len(res_normal_data) > 0 else''
		n_cdk = res_normal_data[0]['cdk'] if len(res_normal_data) > 0 else ''
		
		## Compare two pnr number
		if (t_pnr == n_pnr or t_pnr == ''):
			study_id = n_cdk
			pnr = n_pnr

			extra_cond = "" if study_id in ['None', ''] else " and ec.study_id='{}'".format(study_id)

			query_ecrf = "SELECT CONCAT('MTBP_{}_', ec.study_id,'_{}{}') as identifier, TO_DATE(rf.datum::text, 'YYYYMMDD') as referral_date, '{}' as seq_date, TO_DATE(rf.date_birth::text, 'YYYYMMDD') as birthdate, get_hospital_code(ec.site_id) as hospital, 'oncotree' as cancer_taxonomy, CASE WHEN ec.cancer_type_id != 0 THEN get_tissue_name(ec.cancer_type_id,ec.cancer_type_code) ELSE 'NA' END as tissue, CASE WHEN ec.cancer_type_code !='' THEN ec.cancer_type_code ELSE 'NA' END as cancer_code, 'primary' as tissue_source, CASE WHEN ec.tissue_type ='' THEN 'NA' ELSE ec.tissue_type END as tissue_type, CASE WHEN cell_fraction <> '' THEN CAST(cell_fraction AS FLOAT) ELSE NULL END AS pathology_ccf, CASE WHEN ec.germline_dna = '0' THEN {} ELSE to_number(ec.germline_dna::text, '9'::text)::integer END as germline_dna  from ipcm_referral_t as rf INNER JOIN ipcm_ecrf_t as ec ON regexp_replace(CAST(ec.birth_date AS VARCHAR), '-', '', 'g') =  LEFT(rf.pnr, 8)  WHERE rf.pnr='{}' {} limit 1 ".format(capture_format, sample_type, cfdna, seq_date, germline_dna, pnr, extra_cond)
			res_ecrd_data = fetch_sql_query('ipcmLeaderboard', query_ecrf)
			if len(res_ecrd_data)>0:
				res_json = json.dumps(res_ecrd_data, default = json_serial)
				json_data = json.loads(res_json, parse_int=lambda x: int(x) if x.isdigit() else x)
				identifier_status = True
				res_json_data = json_data[0]
				res_json_data['tissue_type'] = validate_tissue_type(res_json_data['tissue_type'], sample_type)
				res_json_data['pathology_ccf'] = 'NA' if res_json_data['pathology_ccf'] == None else res_json_data['pathology_ccf']
			else:
				identifier_status = False
		else:
			pnr = t_pnr
			query_ecrf_2 = "SELECT CONCAT('MTBP_{}_', ec.study_id,'_{}{}') as identifier, TO_DATE(rf.datum::text, 'YYYYMMDD') as referral_date, '{}' as seq_date, TO_DATE(rf.date_birth::text, 'YYYYMMDD') as birthdate, get_hospital_code(ec.site_id) as hospital, 'oncotree' as cancer_taxonomy, CASE WHEN ec.cancer_type_id != 0 THEN get_tissue_name(ec.cancer_type_id,ec.cancer_type_code) ELSE 'NA' END as tissue, CASE WHEN ec.cancer_type_code !='' THEN ec.cancer_type_code ELSE 'NA' END as cancer_code, 'primary' as tissue_source, CASE WHEN ec.tissue_type ='' THEN 'NA' ELSE ec.tissue_type END as tissue_type, CASE WHEN cell_fraction <> '' THEN CAST(cell_fraction AS FLOAT) ELSE NULL END AS pathology_ccf, CASE WHEN ec.germline_dna = '0' THEN {} ELSE to_number(ec.germline_dna::text, '9'::text)::integer END as germline_dna from ipcm_referral_t as rf INNER JOIN ipcm_ecrf_t as ec ON regexp_replace(CAST(ec.birth_date AS VARCHAR), '-', '', 'g') =  LEFT(rf.pnr, 8)  WHERE rf.pnr='{}'".format(capture_format, sample_type, cfdna, seq_date, germline_dna, pnr)
			res_ecrd_data_2 = fetch_sql_query('ipcmLeaderboard', query_ecrf_2)
			if len(res_ecrd_data_2) > 0:
				res_json_2 = json.dumps(res_ecrd_data_2, default = json_serial)
				json_data_2 = json.loads(res_json_2)
				identifier_status = True
				res_json_data = json_data_2[0]
				res_json_data['tissue_type'] = validate_tissue_type(res_json_data['tissue_type'], sample_type)
				res_json_data['pathology_ccf'] = 'NA' if res_json_data['pathology_ccf'] == None else res_json_data['pathology_ccf']
			else:
				identifier_status = False

		return res_json_data, identifier_status
	
	except Exception as e:
		print("Build iPCM Exception", str(e))
		return res_json_data, identifier_status

def fetch_cancer_type_code(disease):
	cancer_code = 'NA'
	tissue = 'NA'
	sql = "SELECT ts.sub_type_code as cancer_code, tt.tissue_name  FROM ipcm_tissue_subtype as ts INNER JOIN ipcm_tissue_type as tt ON tt.t_id=ts.t_id where ts.sub_type_name ~* '^({})$' order by ts.sub_type_level asc limit 1".format(disease)
	res_data = fetch_sql_query('ipcmLeaderboard', sql)
	if(res_data):
		cancer_code = res_data[0]["cancer_code"]
		tissue = res_data[0]["tissue_name"]
	else:
		tissue = disease
	return tissue, cancer_code

# ### Fetch the genomic profile information
def build_genomic_profile_sample_details(project_name, cfdna, sample_id, capture_id, capture_format, sample_type, seq_date, germline_dna):

	hospital_lookup = { "Karolinska": "KS", "Karolinska Sjukhuset": "KS", "Södersjukhuset": "SO", "St Göran": "ST" }

	sql = "SELECT study_code, study_site, dob, disease FROM genomic_profile_summary where project_name='{}' and sample_id='{}' and capture_id='{}'".format(project_name, sample_id, capture_id)
	res_data = fetch_sql_query('curation', sql)
	res_json = json.dumps(res_data, default = json_serial)

	## Fetch tissue type and germline_dna from capture-id if eCRF information not available from leaderboard
	capture_arr = capture_id.split("_")
	tissue_type =  'cfDNA' if 'CFDNA' in capture_arr[0] else 'NA'

	sample_data = {}

	sample_data["identifier"] = "NA"
	sample_data["referral_date"] = "NA"
	sample_data["seq_date"] = seq_date
	sample_data["birthdate"] = "NA"
	sample_data["hospital"] = "NA"
	sample_data["cancer_taxonomy"] = "oncotree"
	sample_data["tissue"] = "NA"
	sample_data["cancer_code"] = "NA"
	sample_data["tissue_source"] = "primary"


	if(res_data):
		for key, val in enumerate(res_data):
			disease = res_data[key]["disease"]
			study_code = res_data[key]["study_code"]
			dob = res_data[key]["dob"]
			study_site = res_data[key]["study_site"]
		
			if(study_code):
				identifier_name =  "MTBP_"+capture_format+"_"+study_code+"_"+sample_type+cfdna
				sample_data["identifier"] = identifier_name
			if(dob and dob !='NA'):
				pattern = re.compile(r'^\d{4}\d{2}\d{2}$')  # '%Y%m%d'
				if pattern.match(dob):
					dob_format = "%Y%m%d"
				else:
					dob_format = "%Y-%m-%d"

				sample_data["birthdate"] = datetime.strptime(dob.strip(), dob_format).date().strftime("%Y-%m-%d")
			if(study_site and study_site != 'NA'):
				study_site = study_site.strip()
				sample_data["hospital"] = hospital_lookup[study_site] if (study_site in hospital_lookup) else study_site

			if(disease != ""):
				tissue, cancer_code = fetch_cancer_type_code(disease)
				sample_data["tissue"] = tissue
				sample_data["cancer_code"] = cancer_code

	sample_data["tissue_type"] = tissue_type
	sample_data["pathology_ccf"] = "NA"
	sample_data["germline_dna"] = germline_dna

	return sample_data

# ### Fetch cancer code in the genomic profile 
def fetch_cancer_code(project_name,  sample_id, capture_id):
	cancer_code = ''
	sql = "SELECT study_code, disease FROM genomic_profile_summary where project_name='{}' and sample_id='{}' and capture_id='{}' limit 1".format(project_name, sample_id, capture_id)
	res_data = fetch_sql_query('curation', sql)
	if len(res_data) > 0:
		disease = res_data[0]["disease"]
		cancer_code = disease
	return cancer_code


# ### Fetch the sample information from biobank referral table
def build_sample_details(project_name, capture_id, cfdna, capture_format, sample_type, seq_date, germline_dna):

	## Fetch tissue type and germline_dna from capture-id if eCRF information not available from leaderboard
	capture_arr = capture_id.split("_")
	tissue_type =  'cfDNA' if 'CFDNA' in capture_arr[0] else 'NA'

	sample_data = {}

	identifier_status = False

	sample_data["identifier"] = "NA"
	sample_data["referral_date"] = "NA"
	sample_data["seq_date"] = seq_date
	sample_data["birthdate"] = "NA"
	sample_data["hospital"] = "NA"

	sql = "SELECT pnr, CAST(datum AS VARCHAR) as datum, rid, tid from probio_bloodreferrals WHERE cf_dna1 like '%{}%' OR cf_dna2 like '%{}%' or cf_dna3 like '%{}%' or kommentar like '%{}%'".format(cfdna, cfdna, cfdna, cfdna)
	res_data = fetch_sql_query('referral', sql)

	if(res_data):
		sample_data["identifier"] = "MTBP_"+capture_format+"_"+res_data[0]['tid']+"_"+sample_type+cfdna
		sample_data["referral_date"] = res_data[0]["datum"]
		if res_data[0]["pnr"] != None :
			pnr = res_data[0]["pnr"][0:8]
			sample_data["birthdate"] =  datetime.strptime(pnr, "%Y%m%d").date().strftime("%Y-%m-%d")
		else:
			pnr='NA'
			sample_data["birthdate"] ='NA'
		identifier_status = True

		sql = "SELECT subject_id, CAST(dob as VARCHAR), site_name from sample_status_t WHERE pnr like '%{}%'".format(pnr)
		glb_data_1 = fetch_sql_query('leaderboard', sql)
		if(glb_data_1):
			sample_data["birthdate"] = glb_data_1[0]["dob"]
			sample_data["hospital"] = glb_data_1[0]["site_name"]
	else:
		cfdna_rid = re.sub("^0+(?!$)", "", cfdna)
		sql = "SELECT DISTINCT rid, subjectid from biobank_t WHERE regexp_replace(referenceID, '[^a-zA-Z0-9]+', '','g') like '{}' or regexp_replace(referenceID, '[^a-zA-Z0-9]+', '','g') like '{}'  or regexp_replace(referenceID, '[^a-zA-Z0-9]+', '','g') IN ('{}') or regexp_replace(rid, '[^a-zA-Z0-9]+', '','g') IN('{}')".format(cfdna, cfdna, cfdna, cfdna_rid)
		bio_data = fetch_sql_query('leaderboard', sql)
		if(bio_data):
			subject_id = bio_data[0]['subjectid']
			sql = "SELECT subject_id, CAST(dob as VARCHAR), site_name from sample_status_t WHERE regexp_replace(subject_id, '[^a-zA-Z0-9]+', '','g') like regexp_replace('{}', 'P-', '','g') or regexp_replace(subject_id, '[^a-zA-Z0-9]+', '','g') like regexp_replace('{}', '-', '','g')".format(subject_id, subject_id)
			glb_data_2 = fetch_sql_query('leaderboard', sql)

			if(glb_data_2):
				sample_data["birthdate"] = glb_data_2[0]["dob"]
				sample_data["hospital"] = glb_data_2[0]["site_name"]

	sample_data["cancer_taxonomy"] = "oncotree"
	sample_data["tissue"] = "NA"
	sample_data["cancer_code"] = "NA"
	sample_data["tissue_source"] = "primary"
	sample_data["tissue_type"] = tissue_type
	sample_data["pathology_ccf"] = "NA"
	sample_data["germline_dna"] = germline_dna

	return sample_data, identifier_status


# ### Fetch MSI & TMB from genomic profile information
def build_penotypes(project_name, sample_id, capture_id):
	msi_status = "NA"
	tumörmutationsbörda = 'NA'
	pathogenic_gDNA_variant = 'NA'
	td='NA'
	purity_val ='NA'
	ploidy_val = 'NA'

	sql = "SELECT ctdna_param, ctdna_method, genome_wide FROM genomic_profile_summary where project_name='{}' and sample_id='{}' and capture_id='{}'".format(project_name, sample_id, capture_id)
	res_data = fetch_sql_query('curation', sql)
	if(len(res_data) > 0):
		genome_wide = res_data[0]['genome_wide'].replace("\'", "\"")
		genome_wide_json = json.loads(genome_wide)

		for j in genome_wide_json:
			title = genome_wide_json[j]['title']
			if (title == "FRACTION OF CANCER DNA"):
				purity_val = genome_wide_json[j]['result'] if ('result' in genome_wide_json[j] and genome_wide_json[j]['result'] != [] and genome_wide_json[j]['result'] != "") else 'NA'
			elif (title == "PLOIDY"):
				ploidy_val = genome_wide_json[j]['result'] if ('result' in genome_wide_json[j] and genome_wide_json[j]['result'] != [] and genome_wide_json[j]['result'] != "") else 'NA'
			elif (title == "TUMOR MUTATIONAL BURDEN"):
				tumörmutationsbörda = genome_wide_json[j]['result'] if ('result' in genome_wide_json[j] and genome_wide_json[j]['result'] != [] and genome_wide_json[j]['result'] != "") else 'NA'
			elif (title == "MSI STATUS"):
				msi_status = genome_wide_json[j]['result'] if ('result' in genome_wide_json[j] and genome_wide_json[j]['result'] != [] and genome_wide_json[j]['result'] != "") else 'NA'
			elif (title == "PATHOGENIC GERMLINE VARIANTS"):
				pathogenic_gDNA_variant = genome_wide_json[j]['result'] if ('result' in genome_wide_json[j] and genome_wide_json[j]['result'] != [] and genome_wide_json[j]['result'] != "") else 'NA'
			elif (title == "OTHER GENOMIC PHENOTYPE"):
				td = genome_wide_json[j]['result'] if ('result' in genome_wide_json[j] and genome_wide_json[j]['result'] != [] and genome_wide_json[j]['result'] != "") else 'NA'

	return purity_val, ploidy_val, msi_status, tumörmutationsbörda, td, pathogenic_gDNA_variant


# ### Fetch the pipeline version from autoseq-snakemake folder
def fetch_pipeline_version():
	params = readConfig('pipeline')
	dir_path = params['path']
	if os.path.exists(dir_path):
		file_path = dir_path + '.bumpversion.cfg'
		config = ConfigParser()
		config.read_file(open(r'{}'.format(file_path)))
		curr_version = config.get('bumpversion', 'current_version')
		return curr_version
	else:
		return ''

# ### Build a QC Json 
def build_qc(root_path, ecrf_tissue_type, sample_type):
	file_path = root_path + '/qc/'

	msi_list = ''
	qc_df_data = ''

	try:
		qc_file_list = list(filter(lambda x: x.endswith('.qc_overview.txt') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))

		if(len(qc_file_list) > 0):
			qc_filename = file_path + qc_file_list[0]

			if len(qc_filename) > 0:
				qc_df_data = pd.read_csv(qc_filename, delimiter = "\t")

				column_list = ['SAMP', 'MEAN_TARGET_COVERAGE', 'contamination_%']
				column_dict = {'SAMP': 'sample', 'MEAN_TARGET_COVERAGE': 'coverage', 'contamination_%': 'contamination'}

				if 'Overall_QC' in qc_df_data.columns:
					column_list.append('Overall_QC')
					column_dict['Overall_QC'] = 'overall'
				else:
					column_list.append('overall')
					qc_df_data['overall'] = "Ok"

				# if 'MEAN_TARGET_COVERAGE' in qc_df_data.columns:
				# 	column_list.append('MEAN_TARGET_COVERAGE')
				# 	column_dict['MEAN_TARGET_COVERAGE'] = 'coverage'

				msi_list = "NA"
				if 'msing_score' in qc_df_data.columns:
					msi_arr = qc_df_data['msing_score'].dropna().tolist()
					msi_list = msi_arr[0] if len(msi_arr) > 0 else 'NA'

				qc_df_data = qc_df_data[column_list]
				qc_df_data = qc_df_data.rename(columns=column_dict)

				qc_df_data["contamination"] = qc_df_data["contamination"] / 100

				if 'coverage' in qc_df_data.columns:
					qc_df_data["coverage"] = qc_df_data["coverage"].round(0).astype(int)
				else:
					qc_df_data["coverage"] = 'NA'
				
				for idx, row in qc_df_data.iterrows():
					sample_name = row['sample']

					if '-N-' in sample_name:
						qc_df_data.loc[idx,'tissue_type'] = "gDNA"
					elif '-CFDNA-' in sample_name:
						if ecrf_tissue_type != 'NA' and ecrf_tissue_type !='':
							qc_df_data.loc[idx,'tissue_type'] = validate_tissue_type(ecrf_tissue_type, sample_type)
						else:
							qc_df_data.loc[idx,'tissue_type'] =  "cfDNA"
					elif '-T-' in sample_name:
						if ecrf_tissue_type != 'NA' and ecrf_tissue_type !='':
							qc_df_data.loc[idx,'tissue_type'] = validate_tissue_type(ecrf_tissue_type, sample_type)
						else:
							qc_df_data.loc[idx,'tissue_type'] = 'NA'

					if 'overall' in row:
						qc_status = row['overall']
						if qc_status == 'Borderline':
							qc_df_data["overall"] = "Ok"
				
				qc_df_data.fillna('NA', inplace=True)
				
				qc_json = qc_df_data.to_json(orient = 'index')

				return msi_list, qc_json

	except Exception as e:
		print("Build QC Exception", str(e))
		#raise

	return msi_list, qc_df_data


# ### Fetch a Ploidy information from purecn
def build_ploidy(root_path):
	file_path = root_path + '/purecn/'

	regex = '[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+.csv'

	ploidy_list = "NA"
	purity_list = "NA"

	try:
		if os.path.exists(file_path):
			purecn_file_list = list(filter(lambda x: (re.match(regex, x) ),os.listdir(file_path)))
			if len(purecn_file_list) > 0:
				purecn_filename = file_path + purecn_file_list[0]
				if(os.stat(purecn_filename).st_size != 0):
					purecn_df_data = pd.read_csv(purecn_filename, delimiter = ",")
					if(purecn_df_data.empty):
						purecn_df_data = purecn_df_data[['Sampleid', 'Purity', 'Ploidy', 'Sex', 'Contamination']]
						ploidy_list = round(purecn_df_data['Ploidy'].tolist()[0],4)
						purity_list = round(purecn_df_data['Purity'].tolist()[0],4)

	except Exception as e:
		print("Build Polidy Exception", str(e))
		#raise

	return purity_list, ploidy_list


# ### Build a Small variant Json (Somatic & Germline)
def build_small_variants(root_path):

	file_path = root_path + '/'
	smv_df = pd.DataFrame()

	try:
		smv_file_list = list(filter(lambda x: x.endswith('-igvnav-input.txt') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))
		regex = '^(?:(?!-(CFDNA|T)).)*igvnav-input.txt$'
		regex2 = '(.*)-(CFDNA|T)-(\w.*)(germline-igvnav-input).*txt$'

		for i in smv_file_list:
			smv_filename = file_path + i
			smv_df_data = pd.read_csv(smv_filename, delimiter = "\t")

			if 'CALL' in smv_df_data.columns:

				column_list = ['CHROM', 'END', 'REF', 'ALT']
				if(re.match(regex, i) or re.match(regex2, i)):
					reads_column_list = ['N_DP', 'N_ALT']
					column_dict = {'CHROM': 'chr', 'END': 'pos', 'REF' : 'ref', 'ALT' : 'alt', 'N_DP' : 'ref_reads',  'N_ALT' : 'alt_reads'}
				else:
					reads_column_list = ['T_DP', 'T_ALT']
					column_dict = {'CHROM': 'chr', 'END': 'pos', 'REF' : 'ref', 'ALT' : 'alt', 'T_DP' : 'ref_reads',  'T_ALT' : 'alt_reads'}

				column_list.extend(reads_column_list)

				smv_df_data = smv_df_data.loc[(smv_df_data['CALL'] == "S") | (smv_df_data['CALL'] == "G")]
				#smv_df_data["CHROM"] = smv_df_data["CHROM"].map(int)
				#smv_df_data["CHROM"] =  smv_df_data['CHROM'].astype(str).str.isdigit().map(int)
				#smv_df_data['CHROM'] = pd.to_numeric(smv_df_data['CHROM'], errors='coerce')


				if 'CLONALITY' in smv_df_data.columns:
					column_list.append('CLONALITY')
					column_dict['CLONALITY'] = 'clonality'

				if 'SECONDHIT' in smv_df_data.columns:
					column_list.append('SECONDHIT')
					column_dict['SECONDHIT'] = 'second_hit'

				smv_df_data = smv_df_data[column_list]
				smv_df_data = smv_df_data.rename(columns=column_dict)

				if(re.match(regex, i)):
					smv_df_data['origin'] = "Germline"
					smv_df_data['ref_reads'] = smv_df_data['ref_reads'] - smv_df_data['alt_reads']
				else:
					smv_df_data['origin'] = "Somatic"
					smv_df_data['ref_reads'] = smv_df_data['ref_reads'] - smv_df_data['alt_reads']
					if 'clonality' in smv_df_data.columns:
						smv_df_data["clonality"] = smv_df_data["clonality"].fillna('')
					smv_df_data["clonality"] = smv_df_data["clonality"].apply(lambda x: x.capitalize()) if 'clonality' in smv_df_data.columns else ''

				smv_df_data['alt'] = smv_df_data['alt'].map(lambda x: x.lstrip('[').rstrip(']'))

				smv_df_data['strand'] = '+'
				if 'second_hit' in smv_df_data.columns:
					smv_df_data['second_hit'] = smv_df_data['second_hit'].replace('-', 'NA')

			smv_df = pd.concat([smv_df_data, smv_df])

		smv_df.fillna('NA', inplace=True)
		smv_df.reset_index(drop=True, inplace=True)
		return smv_df.to_json(orient = 'index')

	except Exception as e:
		print("Build Small Variants Exception", str(e))
		raise

	return smv_df.to_json(orient = 'index')


# ### Build a CNV Json (Somatic & Germline)
def build_cnv(root_path):

	cnv_df = pd.DataFrame()
	file_path = root_path + '/cnv/'

	try:
		cnv_file_list = list(filter(lambda x: x.endswith('_curated.cns') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))

		if len(cnv_file_list) > 0:
			regex = '^(?:(?!-(CFDNA|T)).)*_curated.cns$'
			for i in cnv_file_list:
				cnv_filename = file_path + i
				cnv_df_data = pd.read_csv(cnv_filename, delimiter = "\t")

				column_list = ['chromosome', 'start', 'end', 'gene', 'ASSESSMENT', 'origin']
				column_dict = {'chromosome': 'chr', 'gene': 'genes', 'ASSESSMENT': 'type'}

				if 'ASSESSMENT' in cnv_df_data.columns:
					cnv_df_data = cnv_df_data.loc[(cnv_df_data['ASSESSMENT'].notnull())]
					if cnv_df_data.empty:
						print("CNV Data frame empty : {}".format(i))
					else:
						# cnv_df_data['chromosome'] = pd.to_numeric(cnv_df_data['chromosome'], errors='coerce')

						if(re.match(regex, i)):
							cnv_df_data['origin'] = "Germline"
						else:
							cnv_df_data['origin'] = "Somatic"

						if 'COPY_NUMBER' in cnv_df_data.columns:
							column_list.append('COPY_NUMBER')
							column_dict['COPY_NUMBER'] = 'copy_number'
						else:
							column_list.append('copy_number')
							cnv_df_data['copy_number'] = 0

						cnv_df_data = cnv_df_data[column_list]
						cnv_df_data = cnv_df_data.rename(columns=column_dict)
						# "NAN" to 0 in the copy number column
						cnv_df_data['copy_number'] = cnv_df_data['copy_number'].fillna(0)

						cnv_df_data['copy_number'] = cnv_df_data['copy_number'].round(0).astype(int,  errors='ignore')
						# cnv_df_data['genes'] = cnv_df_data.genes.apply(lambda x: 'NA' if x.isnull() else x.split(', '))
						cnv_df_data['genes'] = cnv_df_data['genes'].apply(lambda x: 'NA' if pd.isna(x) else x.split(', '))

						cnv_df = pd.concat([cnv_df_data, cnv_df])

		cnv_df.fillna('NA', inplace=True)
		cnv_df.reset_index(drop=True, inplace=True)
		return cnv_df.to_json(orient = 'index')

	except Exception as e:
		print("Build CNV Exception", str(e))
		#raise

	return cnv_df.to_json(orient = 'index')


# ### Build a SVS Json
def build_svs(root_path):
	
	file_path = root_path + '/svs/igv/'

	svs_filter = pd.DataFrame()

	try:
		svs_file_list = list(filter(lambda x: (re.match('[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+-sv-annotated.txt', x) or x.endswith('_annotate_combined_SV.txt')) and not x.startswith('.') and not x.endswith('.out'),os.listdir(file_path)))

		if len(svs_file_list) > 0:
			svs_filename = file_path + svs_file_list[0]

			sample_list = ['germline', 'somatic', 'tumor', 'cfdna']
			svs_filter = pd.read_csv(svs_filename, delimiter = "\t")

			column_list = ['SAMPLE', 'SVTYPE', 'strand', 'vaf', 'GENE_A', 'IGV_COORD', 'GENE_B', 'CLONALITY', 'SECONDHIT']
			column_dict = {'SAMPLE': 'origin', 'SVTYPE': 'sv_type', 'GENE_A': 'gene_a' , 'IGV_COORD' : 'variant', 'GENE_B' : 'gene_b', 'CLONALITY': 'clonality', 'SECONDHIT' : 'second_hit'}

			if 'CALL' in svs_filter.columns:
				svs_filter = svs_filter.loc[(svs_filter['CALL'] == True ) | (svs_filter['CALL'] == 'true')]

				if not svs_filter.empty:
					svs_filter = svs_filter[['SAMPLE', 'SVTYPE', 'IGV_COORD', 'GENE_A', 'GENE_B','SECONDHIT', 'CLONALITY', 'CONSEQUENCE', 'VARIANT_STRING']]

					svs_filter = svs_filter[svs_filter['SAMPLE'].isin(sample_list)]

					svs_filter["SAMPLE"] = np.where(svs_filter["SAMPLE"] == "germline", 'germline', 'somatic')

					svs_filter['IGV_COORD'] = svs_filter['IGV_COORD'].str.strip()

					svs_filter['IGV_COORD'] = svs_filter['IGV_COORD'].str.split(' ')

					svs_filter["strand"] = '+'
					svs_filter["vaf"] = ''

					if 'CONSEQUENCE' in svs_filter.columns:
						column_list.append('CONSEQUENCE')
						column_dict['CONSEQUENCE'] = 'consequence'
					else:
						svs_filter["consequence"] = 'NA'

					if 'VARIANT_STRING' in svs_filter.columns:
						column_list.append('VARIANT_STRING')
						column_dict['VARIANT_STRING'] = 'variant_string'
					else:
						svs_filter["variant_string"] = "NA"

					svs_filter = svs_filter[column_list]
					svs_filter = svs_filter.rename(columns=column_dict)

					svs_filter.fillna('NA', inplace=True)
					svs_filter["origin"] = svs_filter["origin"].apply(lambda x: x.capitalize())
					svs_filter["clonality"] = svs_filter["clonality"].apply(lambda x: x.capitalize())
					svs_filter["consequence"] = svs_filter["consequence"].apply(lambda x: x.capitalize())

			else:
				svs_filter = pd.DataFrame()

	except Exception as e:
		print("Build SVS Exception", str(e))
		# svs_filter = pd.DataFrame()
		raise

	return svs_filter.to_json(orient = 'index')

def check_cancer_code(json_file_path, key, expected_value):
	# Load JSON data from file
	with open(json_file_path, 'r') as file:
		data = json.load(file)
	
	cancer_code_status = False
	if key in data['sample']:
		# Check if the value for the key matches the expected value
		if data['sample'][key] == expected_value:
			cancer_code_status = False
		else:
			cancer_code_status = True

	return cancer_code_status

def upload_json_MTBP(nfs_path, json_file_path, project_name, sample_id, cfdna, identifier_study_id):

	mtbp_file_name = nfs_path + "/MTBP_report_status.csv"

	if os.path.isfile(mtbp_file_name) == False:
		with open(mtbp_file_name, 'w', newline='', encoding='utf-8') as file:
			writer = csv.writer(file)
			writer.writerow(['PROJECT', 'SAMPLE_ID', 'BARCODE', 'STUDY_ID', 'MTBP_JSON_PATH', 'CURL_STATUS'])
	
	curl_res = ''

	if os.path.isfile(json_file_path):
		user_name = os.environ['user_name']
		user_pwd = os.environ['user_pwd']
		key = 'cancer_code'
		expected_value ='NA'
		cancer_code_status = check_cancer_code(json_file_path, key, expected_value)
		if cancer_code_status:
			curl_cmd = "curl -sk --form 'fileToUpload=@{}' --form 'username={}' --form 'password={}' --form 'Proj=1' --form 'seqdata=1' https://cloud-mtb.scilifelab.se/UploadClinicalDataAPI.php".format(json_file_path, user_name, user_pwd)
			logger2.info("message : {}".format(curl_cmd))
			logger1.info('Curl Command {}'.format(curl_cmd))
			proc = subprocess.run(curl_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, capture_output=False, text=True)
			output= proc.stdout
			res_error = proc.stderr
			res_status_code = proc.returncode
			if res_status_code:
				logger2.error("message : {}".format(str(res_error)))
				curl_res = str(res_error)
			else:
				logger2.info("message : {}".format(str(output)))
				curl_res = str(output)
		else:
			curl_res = 'Cancer Code NA, Please check in the iPCM Leaderboard'
			logger2.info("message : {}".format(str(curl_res)))

		data= [project_name,  sample_id, cfdna, identifier_study_id, json_file_path, curl_res]
		with open(mtbp_file_name, 'a+') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(data)
	else:
		logger2.error("{} : Json File not generated".format(json_file_path))

	return 0
	
# ## Build Json from output files
def build_json(nfs_path, root_path, output_path, project_name, normal_cfdna, cfdna, sample_id, capture_id, capture_format, sample_type, seq_date, file_handler2):

	print("--- MTBP Json Format Started ---\n")

	## Check the sample was germline_dna or not from capture-id with sample-id
	capture_arr = capture_id.split("_")
	sdid_0 = sample_id in capture_arr[0]
	sdid_1 = sample_id in capture_arr[1]
	germline_dna = 1 if sdid_0 == sdid_1 else 0

	project_json = {}

	# Sample Information
	logger2.info('--- Sample fetching started ---')
	itendifiter_status = True
	if(project_name == "IPCM" or capture_format == "iPCM"):
		sample_details_json, itendifiter_status = build_ipcm_sample_details(normal_cfdna, cfdna, capture_format, sample_type, seq_date, germline_dna)
	else:
		sample_details_json, itendifiter_status = build_sample_details(project_name, capture_id, cfdna, capture_format, sample_type, seq_date, germline_dna)

	if(sample_details_json and itendifiter_status):
		project_json["sample"] = sample_details_json
	else:
		sample_details_json = build_genomic_profile_sample_details(project_name, cfdna, sample_id, capture_id, capture_format, sample_type, seq_date, germline_dna)
		project_json["sample"] = sample_details_json

	if sample_details_json["identifier"] == "NA" :
		identifier_study_id = sample_id
	else:
		identifier_study_id = sample_details_json["identifier"].split("_")[2]

	if sample_details_json["cancer_code"] == "NA" :
		cancer_code = fetch_cancer_code(project_name, sample_id, capture_id)
		if cancer_code != '':
			sample_details_json["cancer_code"] = cancer_code
	
	ecrf_tissue_type = sample_details_json['tissue_type']
	logger1.info("Study id : {}".format(identifier_study_id))
	logger1.info("Cancer code : {}".format(sample_details_json['cancer_code']))

	logger2.info('--- Sample fetching completed ---')

	# Pipeline
	curr_version = fetch_pipeline_version()
	piln_version = curr_version if curr_version != '' else '3.2.0'
	project_json["pipeline"] = { "genome_reference": "hg19", "version": piln_version}

	# QC
	logger2.info('--- QC started ---')
	msi_val, qc_json = build_qc(root_path, ecrf_tissue_type, sample_type)
	project_json["qc"] =  "NA" if qc_json == "" else json.loads(qc_json)
	logger2.info('--- QC completed ---')

	# Get the Ploidy value
	#purity_val, ploidy_val = build_ploidy(root_path)

	# Phenotype
	purity_val, ploidy_val, msi, tmb, td, pathogenic_gDNA_variant = build_penotypes(project_name, sample_id, capture_id)
	project_json["phenotypes"] = {"purity" : convert_to_numeric(purity_val) ,"ploidy": convert_to_numeric(ploidy_val), "msi": msi, "tmb": tmb, "td": td, "pathogenic_gDNA_variant": pathogenic_gDNA_variant}

	# Small Variant (Somatic & Germline)
	logger2.info('--- Small Variant started ---')
	small_variant_json = build_small_variants(root_path)
	project_json["small_variants"] = json.loads(small_variant_json)
	logger2.info('--- Small Variant completed ---')

	# CNVs
	logger2.info('--- CNVs started ---')
	cnv_json = build_cnv(root_path)
	project_json["cnas"] = json.loads(cnv_json)
	logger2.info('--- CNVs completed ---')

	# SVS
	logger2.info('--- SVS started ---')
	svs_json = build_svs(root_path)
	project_json["gsr"] = json.loads(svs_json)
	logger2.info('--- SVS completed ---')

	#final_json = json.dumps(project_json)

	file_name = output_path+"/MTBP_"+capture_format+"_"+identifier_study_id+"_"+sample_type+cfdna+".json"

	print("Path : ", root_path, file_name)

	with open(file_name, 'w') as f:
		json.dump(project_json, f, indent=4)

	logger2.info('--- Generated Json format successfuly ---\n')

	print("\n----  MTBP Json Format Completed -----\n")
	upload_json_MTBP(nfs_path, file_name, project_name, sample_id, cfdna, identifier_study_id)
	logger1.info('Sample {} Completed'.format(sample_id))
	logger1.info('================\n')
	logger2.removeHandler(file_handler2)

def main(nfs_path, project_name, capture_format, sample_ids):
	
	regx_pat = '^P-.*[0-9]$'
	regx_capture_pat = '^{}-P-.*[0-9]$'.format(capture_format)

	overall_log_name = nfs_path+"/MTBP_"+project_name+"_overall.log"

	file_handler1 = logging.FileHandler(overall_log_name)
	file_handler1.setLevel(logging.INFO)
	file_handler1.setFormatter(formatter)
	logger1.addHandler(file_handler1)

	# formatter = logging.basicConfig(format = '%(asctime)s  %(levelname)-10s %(name)s %(message)s', level=logging.INFO , datefmt="%Y-%m-%d %H:%M:%S")

	logger1.info('Overall MTBP JSON')

	for sid in sample_ids:
		sample_id = sid
		logger1.info('Sample ID : {}'.format(sample_id))
		sample_dir_path = os.path.join(nfs_path, sid)
		if os.path.isdir(sample_dir_path):
			inner_dir_list = os.listdir(sample_dir_path)
			if(len(inner_dir_list) > 0):
				for indir in inner_dir_list:
					capture_id = indir
					logger1.info('Capture ID : {}'.format(capture_id))
					capture_check = re.match(regx_capture_pat, capture_id)
					if(capture_check):
						cfdna = capture_id.split("-")[4]
						root_path = os.path.join(sample_dir_path,capture_id)
						output_path = root_path+"/MTBP";
						
						capture_arr = capture_id.split("-")

						normal_idx = capture_arr.index("N")
						normal_cfdna = capture_arr[normal_idx+1]

						cfdna_idx = capture_arr.index("T") if 'T' in capture_arr else capture_arr.index("CFDNA")
						sample_type = capture_arr[cfdna_idx]
						cfdna_id = capture_arr[cfdna_idx+1]

						cfdna = re.sub(r'[a-zA-Z]', '', cfdna_id)

						capture_format = capture_arr[0]

						capture_arr = capture_id.split("-")
						seq_date_str = re.findall(r'\d+', capture_arr[5])[0]

						if len(seq_date_str) != 8:
							seq_date_str = '20' + str(seq_date_str)
						
						seq_date = datetime.strptime(seq_date_str, "%Y%m%d").date().strftime("%Y-%m-%d")
						
						log_name = output_path+"/MTBP_"+project_name+"_"+sample_id+"_"+sample_type+cfdna+".log"

						if(not os.path.exists(output_path)):
							os.mkdir(output_path)
						else:
							for f in os.listdir(output_path):
								os.remove(os.path.join(output_path, f))

						file_handler2 = logging.FileHandler(log_name)
						file_handler2.setLevel(logging.INFO)
						file_handler2.setFormatter(formatter)
						logger2.addHandler(file_handler2)

						# formatter2 = logging.basicConfig(format = '%(asctime)s  %(levelname)-10s %(name)s %(message)s', level=logging.INFO ,  datefmt="%Y-%m-%d %H:%M:%S")

						logger2.info('--- Generated Json format Started---')

						logger2.info("Sample Id : {} || Capture Id : {} ".format(sample_id,capture_id))
						try:
							print("Sample Id : {} || Capture Id : {} ".format(sample_id,capture_id))
							build_json(nfs_path, root_path, output_path, project_name, normal_cfdna, cfdna, sample_id, capture_id, capture_format, sample_type, seq_date, file_handler2)
						except Exception as e:
							print("Exception", str(e))
							logger2.error("Failed : {}".format(str(e)))
							logger2.error('--- Generated Json format Failed ---')
							raise
			else:
				logger1.info("No Files & Folder found for {}".format(sid))
				logger1.info('-------------------')
		else:
			logger1.info("Sample {} not found".format(sid))
			logger1.info('-------------------')
		
	
	
	
if __name__ == "__main__":
		
	# Create the parser
	profile_parser = argparse.ArgumentParser(description='Generate MTBP Json')

	# Add the arguments
	profile_parser.add_argument('--path', metavar='nfs root path', type=str, help='define the sample autoseq-output path')
	profile_parser.add_argument('--project_name', metavar='project name', type=str, help='define the project name (eg: PROBIO, PSFF, IPCM )')
	profile_parser.add_argument('--capture_format', metavar='project name', type=str, help='define the capture name (eg: PB, PSFF, iPCM )')
	profile_parser.add_argument('--sample_ids', metavar='Sample Id List', type=str, nargs='+', action='store', help='define the capture name (eg: ["P-002546", "P-00898"] )')
	
	args = profile_parser.parse_args()

	nfs_path = args.path
	project_name = args.project_name
	capture_format = args.capture_format
	sample_ids = args.sample_ids

	print("sample_ids", sample_ids)

	if not os.path.isdir(nfs_path):
		print('The path specified does not exist')
		sys.exit()
	else:
		main(nfs_path,project_name, capture_format, sample_ids)
	
