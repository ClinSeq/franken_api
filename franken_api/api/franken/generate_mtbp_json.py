#!/usr/bin/env python
# coding: utf-8

### File name format for IPCM and Default 
#    **MTBP_iPCM_KI_3405_CFDNA03884517**
#
#     *MTBP - molecular tumor board portal
#     iPCM - clinical trial/study/project
#     KI - Karolinska Hospital (site) (other sites ST - St Goran, SOS - Sodersjuksjuset)
#     3405 - Patient Identifier
#     CFDNA03884517 - sample type with ID number*

### IDENTIFIER Format for PREDDLUNG
#    ** MOL2347-25_T-DNA **
#     MOL432-25_
#
# Starts with: MOL, followed by: one or more digits
# Then: a hyphen (-)
# Then: two digits (year)
# Then: an optional suffix starting with _.  The suffix can contain:
#  * Uppercase and lowercase letters and digits
#  * URL-safe special characters (it must not contain characters that require percent-encoding in a URL)

import os
import re
from datetime import date, datetime
import logging
from decimal import Decimal
import json
from configparser import ConfigParser
import pandas as pd # type: ignore
import numpy as np
from pathlib import Path

from franken_api.database import db
from flask import current_app

from sqlalchemy.exc import SQLAlchemyError # type: ignore
from sqlalchemy import and_, or_, text
from sqlalchemy.orm import scoped_session, sessionmaker

SessionFactory = sessionmaker()
Session = scoped_session(SessionFactory)


def generate_list_to_dict(result):
	d, row = {}, []
	for rowproxy in result:
		#row_as_dict = rowproxy
		row_as_dict = rowproxy._mapping
		for column, value in row_as_dict.items():
			if isinstance(value, datetime):
				value = datetime.strftime(value, "%Y-%m-%d")
			else:
				value = str(value)
			d = {**d, **{column: value}}
		row.append(d)
	
	return row
	

def create_db_session(db_name, query):
	session = None
	try:
		# engine1 = db.get_engine(db_name)
		engine = db.get_engine(current_app, db_name) if db_name !='' else db.get_engine(current_app)
		SessionFactory.configure(bind=engine)
		session = Session() 
		sql_query = text(query)
		res = session.execute(sql_query)
		session.commit()
		return res
	except SQLAlchemyError as e:
		if session:
			session.rollback()
		raise
	finally:
		if session:
			Session.remove()


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

def fetch_cancer_type_code(disease):
	cancer_code = 'NA'
	tissue = 'NA'
	sql = "SELECT ts.sub_type_code as cancer_code, tt.tissue_name  FROM ipcm_tissue_subtype as ts INNER JOIN ipcm_tissue_type as tt ON tt.t_id=ts.t_id where ts.sub_type_name ~* '^({})$' order by ts.sub_type_level asc limit 1".format(disease)
	res = create_db_session('ipcmLeaderboard', sql)
	res_data = generate_list_to_dict(res)

	if(res_data):
		cancer_code = res_data[0]["cancer_code"]
		tissue = res_data[0]["tissue_name"]
	else:
		tissue = disease
	return tissue, cancer_code

### Fetch the sample information from ipcm referral table
def build_ipcm_sample_details(normal_cfdna, cfdna, capture_format, sample_type, seq_date, germline_dna):

	identifier_status = False
	res_json_data = []
	try:

		## tissue cfdna
		cfdna = re.sub(r'[a-zA-Z]', '', cfdna)
		sql = "select rf.pnr, rf.cdk from ipcm_referral_t as rf WHERE rf.rid like'%{}%' OR rf.blood like'%{}%' OR rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' OR rf.dna3 like'%{}%'".format(cfdna, cfdna, cfdna, cfdna, cfdna)
		res = create_db_session('ipcmLeaderboard', sql)
		res_data = generate_list_to_dict(res)

		## normal cfdna
		sql2 = "select rf.pnr, rf.cdk from ipcm_referral_t as rf WHERE rf.rid like'%{}%' OR rf.blood like'%{}%' OR rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' OR rf.dna3 like'%{}%'".format(normal_cfdna, normal_cfdna, normal_cfdna, normal_cfdna, normal_cfdna)
		res_2 = create_db_session('ipcmLeaderboard', sql2)
		res_normal_data = generate_list_to_dict(res_2)
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
			res_ecrf = create_db_session('ipcmLeaderboard', query_ecrf)
			res_ecrd_data = generate_list_to_dict(res_ecrf)

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
			res_ecrf_2 = create_db_session('ipcmLeaderboard', query_ecrf_2)
			res_ecrd_data_2 = generate_list_to_dict(res_ecrf_2)

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

### Fetch the sample information from biobank referral table
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
	res = create_db_session('', sql)
	res_data = generate_list_to_dict(res)

	if(res_data):
		sample_data["identifier"] = "MTBP_"+capture_format+"_"+res_data[0]['tid']+"_"+sample_type+cfdna
		sample_data["referral_date"] = res_data[0]["datum"]
		if res_data[0]["pnr"] != None and res_data[0]["pnr"] != 'None':
			pnr = res_data[0]["pnr"][0:8]
			sample_data["birthdate"] =  datetime.strptime(pnr, "%Y%m%d").date().strftime("%Y-%m-%d")
		else:
			pnr='NA'
			sample_data["birthdate"] ='NA'
		identifier_status = True

		sql2 = "SELECT subject_id, CAST(dob as VARCHAR), site_name from sample_status_t WHERE pnr like '%{}%'".format(pnr)
		res_glb = create_db_session('leaderboard', sql2)
		glb_data_1 = generate_list_to_dict(res_glb)

		if(glb_data_1):
			sample_data["birthdate"] = glb_data_1[0]["dob"]
			sample_data["hospital"] = glb_data_1[0]["site_name"]
	else:
		cfdna_rid = re.sub("^0+(?!$)", "", cfdna)
		sql_bio = "SELECT DISTINCT rid, subjectid from biobank_t WHERE regexp_replace(referenceID, '[^a-zA-Z0-9]+', '','g') like '{}' or regexp_replace(referenceID, '[^a-zA-Z0-9]+', '','g') like '{}'  or regexp_replace(referenceID, '[^a-zA-Z0-9]+', '','g') IN ('{}') or regexp_replace(rid, '[^a-zA-Z0-9]+', '','g') IN('{}')".format(cfdna, cfdna, cfdna, cfdna_rid)
		res_bio = create_db_session('leaderboard', sql_bio)
		bio_data = generate_list_to_dict(res_bio)

		if(bio_data):
			subject_id = bio_data[0]['subjectid']
			sql_glb_2 = "SELECT subject_id, CAST(dob as VARCHAR), site_name from sample_status_t WHERE regexp_replace(subject_id, '[^a-zA-Z0-9]+', '','g') like regexp_replace('{}', 'P-', '','g') or regexp_replace(subject_id, '[^a-zA-Z0-9]+', '','g') like regexp_replace('{}', '-', '','g')".format(subject_id, subject_id)
			res_glb_2 = create_db_session('leaderboard', sql_glb_2)
			glb_data_2 = generate_list_to_dict(res_glb_2)

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

### Fetch the genomic profile information
def build_genomic_profile_sample_details(project_name, cfdna, sample_id, capture_id, capture_format, sample_type, seq_date, germline_dna):

	hospital_lookup = { "Karolinska": "KS", "Karolinska Sjukhuset": "KS", "Södersjukhuset": "SO", "St Göran": "ST" }

	sql = "SELECT study_code, study_site, dob, disease FROM genomic_profile_summary where project_name='{}' and sample_id='{}' and capture_id='{}'".format(project_name, sample_id, capture_id)
	res = create_db_session('curation', sql)
	res_data = generate_list_to_dict(res)

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
				identifier_name =  "MTBP_"+capture_format+"_"+study_code+"_"+sample_type+cfdna if project_name !='PREDDLUNG' else "NA"
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

### Fetch cancer code in the genomic profile 
def fetch_cancer_code(project_name,  sample_id, capture_id):
	cancer_code = ''
	sql = "SELECT study_code, disease FROM genomic_profile_summary where project_name='{}' and sample_id='{}' and capture_id='{}' limit 1".format(project_name, sample_id, capture_id)
	res = create_db_session('curation', sql)
	res_data = generate_list_to_dict(res)

	if len(res_data) > 0:
		disease = res_data[0]["disease"]
		cancer_code = disease
	return cancer_code

### Fetch the pipeline version from autoseq-snakemake folder
def fetch_pipeline_version():
	dir_path = '/nfs/PIPELINE/autoseq-snakemake/'
	if os.path.exists(dir_path):
		file_path = dir_path + '.bumpversion.cfg'
		config = ConfigParser()
		config.read_file(open(r'{}'.format(file_path)))
		curr_version = config.get('bumpversion', 'current_version')
		return curr_version
	else:
		return ''
### Build a QC Json 
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

	return msi_list, qc_df_data

### Fetch a Ploidy information from purecn
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

	return purity_list, ploidy_list

def parse_tumor_val(tumor_val):
	if tumor_val == 'NA':
		return ''
	
	tmb_snv = 'NA'
	tmb_indel = 'NA'

	parts = [p.strip() for p in tumor_val.split(',')]
	tmb_dict = {}
	for part in parts:
		if 'SNV' in part:
			match_snv = re.findall(r'\d+\.\d+|\d+', part)
			if match_snv:
				tmb_snv = match_snv[0]
				tmb_dict['tmb_snv'] = int(tmb_snv)
		elif 'INDEL' in part:
			match_indel = re.findall(r'\d+\.\d+|\d+', part)
			if match_indel:
				tmb_indel = match_indel[0]
				tmb_dict['tmb_indel'] = int(tmb_indel)

	if len(tmb_dict) > 0:
		return tmb_dict
	else:
		return ''

### Fetch MSI & TMB from genomic profile information
def build_penotypes(project_name, sample_id, capture_id):
	msi_status = "NA"
	tumörmutationsbörda = 'NA'
	pathogenic_gDNA_variant = 'NA'
	td='NA'
	purity_val ='NA'
	ploidy_val = 'NA'
	tumor_val = ''

	sql = "SELECT ctdna_param, ctdna_method, genome_wide FROM genomic_profile_summary where project_name='{}' and sample_id='{}' and capture_id='{}'".format(project_name, sample_id, capture_id)
	res = create_db_session('curation', sql)
	res_data = generate_list_to_dict(res)

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
				tumor_arr = genome_wide_json[j]['comment'] if ('comment' in genome_wide_json[j] and genome_wide_json[j]['comment'] != [] and genome_wide_json[j]['comment'] != "") else 'NA'
				tumor_val = parse_tumor_val(tumor_arr)
				
			elif (title == "MSI STATUS"):
				msi_status = genome_wide_json[j]['result'] if ('result' in genome_wide_json[j] and genome_wide_json[j]['result'] != [] and genome_wide_json[j]['result'] != "") else 'NA'
			elif (title == "PATHOGENIC GERMLINE VARIANTS"):
				pathogenic_gDNA_variant = genome_wide_json[j]['result'] if ('result' in genome_wide_json[j] and genome_wide_json[j]['result'] != [] and genome_wide_json[j]['result'] != "") else 'NA'
			elif (title == "OTHER GENOMIC PHENOTYPE"):
				td = genome_wide_json[j]['result'] if ('result' in genome_wide_json[j] and genome_wide_json[j]['result'] != [] and genome_wide_json[j]['result'] != "") else 'NA'

	return purity_val, ploidy_val, msi_status, tumörmutationsbörda, td, pathogenic_gDNA_variant, tumor_val

### Build a Small variant Json (Somatic & Germline)
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
				
				smv_df_data["clonality"] = smv_df_data['clonality'].replace('', 'NA')
				
				smv_df_data['strand'] = '+'
				if 'second_hit' in smv_df_data.columns:
					smv_df_data['second_hit'] = smv_df_data['second_hit'].replace('-', 'NA')

			smv_df = pd.concat([smv_df_data, smv_df])

		smv_df.fillna('NA', inplace=True)
		smv_df.reset_index(drop=True, inplace=True)
		return smv_df.to_json(orient = 'index')

	except Exception as e:
		print("Build Small Variants Exception", str(e))

	return smv_df.to_json(orient = 'index')

### Build a CNV Json (Somatic & Germline)
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

	return cnv_df.to_json(orient = 'index')

### Build a SVS Json
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
		svs_filter = pd.DataFrame()

	return svs_filter.to_json(orient = 'index')

### Build Json from output files
def build_json(root_path, output_path, project_name, normal_cfdna, cfdna, sample_id, capture_id, capture_format, sample_type, seq_date):

	print("--- MTBP Json Format Started ---\n")

	## Check the sample was germline_dna or not from capture-id with sample-id
	capture_arr = capture_id.split("_")
	sdid_0 = sample_id in capture_arr[0]
	sdid_1 = sample_id in capture_arr[1]
	germline_dna = 1 if sdid_0 == sdid_1 else 0

	project_json = {}
	sample_details_json = {}
	itendifiter_status = True

	# Sample Information
	logging.info('--- Sample fetching started ---')

	if(project_name == "IPCM" or capture_format == "iPCM"):
		sample_details_json, itendifiter_status = build_ipcm_sample_details(normal_cfdna, cfdna, capture_format, sample_type, seq_date, germline_dna)
	elif project_name != 'PREDDLUNG':
		sample_details_json, itendifiter_status = build_sample_details(project_name, capture_id, cfdna, capture_format, sample_type, seq_date, germline_dna)

	if(sample_details_json and itendifiter_status):
		project_json["sample"] = sample_details_json
	else:
		sample_details_json = build_genomic_profile_sample_details(project_name, cfdna, sample_id, capture_id, capture_format, sample_type, seq_date, germline_dna)
		project_json["sample"] = sample_details_json

	if sample_details_json["identifier"] == "NA" :
		if project_name == 'PREDDLUNG':
			identifier_study_id = sample_id.replace("P-","P")
			sample_details_json["identifier"] = "MTBP_PreDDLung_{}_{}{}".format(identifier_study_id, sample_type, cfdna)
		else:
			identifier_study_id = sample_id
	else:
		identifier_study_id = sample_details_json["identifier"].split("_")[2]

	if sample_details_json["cancer_code"] == "NA" :
		disease = fetch_cancer_code(project_name, sample_id, capture_id)
		if disease != '':
			tissue, cancer_code_new = fetch_cancer_type_code(disease)
			sample_details_json["cancer_code"] = cancer_code_new if cancer_code_new !='NA' else sample_details_json['tissue']

	
	ecrf_tissue_type = sample_details_json['tissue_type']

	logging.info('--- Sample fetching completed ---')

	# Pipeline
	curr_version = fetch_pipeline_version()
	piln_version = curr_version if curr_version != '' else '3.2.0'
	project_json["pipeline"] = { "genome_reference": "hg19", "version": piln_version}

	# QC
	logging.info('--- QC started ---')
	msi_val, qc_json = build_qc(root_path, ecrf_tissue_type, sample_type)
	project_json["qc"] =  "NA" if qc_json == "" else json.loads(qc_json)
	logging.info('--- QC completed ---')

	# Phenotype
	purity_val, ploidy_val, msi, tmb, td, pathogenic_gDNA_variant, tumor_val = build_penotypes(project_name, sample_id, capture_id)
	project_json["phenotypes"] = {"purity" : convert_to_numeric(purity_val) ,"ploidy": convert_to_numeric(ploidy_val), "msi": msi, "tmb": tmb, "td": td, "pathogenic_gDNA_variant": pathogenic_gDNA_variant, **(tumor_val if isinstance(tumor_val, dict) else {})}

	# Small Variant (Somatic & Germline)
	logging.info('--- Small Variant started ---')
	small_variant_json = build_small_variants(root_path)
	project_json["small_variants"] = json.loads(small_variant_json)
	logging.info('--- Small Variant completed ---')

	# CNVs
	logging.info('--- CNVs started ---')
	cnv_json = build_cnv(root_path)
	project_json["cnas"] = json.loads(cnv_json)
	logging.info('--- CNVs completed ---')

	# SVS
	logging.info('--- SVS started ---')
	svs_json = build_svs(root_path)
	project_json["gsr"] = json.loads(svs_json)
	logging.info('--- SVS completed ---')

	# final_json = json.dumps(project_json)
	file_name = f"MTBP_{capture_format}_{identifier_study_id}_{sample_type}{cfdna}.json"
	file_path = os.path.join(output_path,file_name)

	with open(file_path, 'w') as f:
		json.dump(project_json, f, indent=4)

	logging.info('--- Generated Json format successfully ---\n')

	print("\n----  MTBP Json Format Completed -----\n")


### Main Function
def generate_json(nfs_path, project_name, sample_id, capture_id):

	nfs_folder_path = os.path.join(nfs_path, sample_id, capture_id)

	output_dir = os.path.join(nfs_folder_path, "MTBP")

	capture_arr = capture_id.split("-")
	
	normal_idx = capture_arr.index("N")
	normal_cfdna = capture_arr[normal_idx+1]

	cfdna_idx = capture_arr.index("T") if 'T' in capture_arr else capture_arr.index("CFDNA")
	sample_type = capture_arr[cfdna_idx]
	cfdna_id = capture_arr[cfdna_idx+1]

	cfdna = re.sub(r'[a-zA-Z]', '', cfdna_id)

	capture_format = capture_arr[0]
	
	# capture_arr = capture_id.split("-")
	seq_date_str = re.findall(r'\d+', capture_arr[5])[0]
	seq_date = datetime.strptime(seq_date_str, "%Y%m%d").date().strftime("%Y-%m-%d")

	log_name = f"{output_dir}/MTBP_{project_name}_{sample_id}_{sample_type}{cfdna}.log"

	output_path = Path(output_dir)

	if output_path.exists():
		# Remove only files (not subdirectories) inside the folder
		for item in output_path.iterdir():
			if item.is_file():
				item.unlink()
	else:
		output_path.mkdir(parents=True)

	logging.basicConfig(format = '%(asctime)s  %(levelname)-10s %(name)s %(message)s', level=logging.INFO , filename=log_name, datefmt="%Y-%m-%d %H:%M:%S")
	logging.info('--- Generated Json format Started---')

	logging.info("Sample Id : {} || Capture Id : {} ".format(sample_id,capture_id))

	try:
		build_json(nfs_folder_path, output_dir, project_name, normal_cfdna, cfdna, sample_id, capture_id, capture_format, sample_type, seq_date)
	except Exception as e:
		print("Main Exception", str(e))
		logging.error("Failed : {}".format(str(e)))
		logging.error('--- Generated Json format Failed ---\n')
		raise
