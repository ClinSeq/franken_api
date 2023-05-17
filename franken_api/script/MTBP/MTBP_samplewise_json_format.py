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

import os
import pandas as pd
import numpy as np
import json
import re
import logging
import psycopg2
import psycopg2.extras
from sqlalchemy import create_engine
from configparser import ConfigParser
import yaml
import argparse
import sys
from datetime import date, datetime
import math

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
	raise TypeError("Type %s not serializable" % type(obj))


# ### Fetch the sample information from ipcm referral table
def build_icpm_sample_details(normal_cfdna, cfdna):

	identifier_status = False
	res_json_data = ''
	try:

		## tissue cfdna
		cfdna = re.sub(r'[a-zA-Z]', '', cfdna)
		sql = "select rf.pnr from ipcm_referral_t as rf WHERE rf.rid like'%{}%' OR rf.blood like'%{}%' OR rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' OR rf.dna3 like'%{}%'".format(cfdna, cfdna, cfdna, cfdna, cfdna)
		res_data = fetch_sql_query('ipcmLeaderboard', sql)

		## normal cfdna
		sql2 = "select rf.pnr from ipcm_referral_t as rf WHERE rf.rid like'%{}%' OR rf.blood like'%{}%' OR rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' OR rf.dna3 like'%{}%'".format(normal_cfdna, normal_cfdna, normal_cfdna, normal_cfdna, normal_cfdna)
		res_normal_data = fetch_sql_query('ipcmLeaderboard', sql2)

		t_pnr = res_data[0]['pnr'] if len(res_data) else ''
		n_pnr = res_normal_data[0]['pnr'] if len(res_normal_data) else''

		## Compare two pnr number
		pnr = n_pnr if (t_pnr == n_pnr or t_pnr == '') else t_pnr

		if(pnr):
			# dob = pnr[0:8]
			sql3 = "SELECT CONCAT('MTBP_iPCM_', ec.study_id,'_CFDNA{}') as identifier, TO_DATE(rf.datum::text, 'YYYYMMDD') as sample_date, TO_DATE(rf.date_birth::text, 'YYYYMMDD') as birthdate, get_hospital_code(ec.site_id) as hospital, 'oncotree' as cancer_taxonomy, CASE WHEN ec.cancer_type_code !='' THEN ec.cancer_type_code ELSE get_tissue_name(ec.cancer_type_id, ec.cancer_type_code) END as cancer_code,  CASE WHEN ec.cancer_type_code !='' THEN ec.cancer_type_code ELSE 'N/A' END as cancer_sub_code, 'primary' as tissue_source, get_tissue_name(ec.cancer_type_id, ec.cancer_type_code) as tissue_type, ec.cell_fraction as pathology_ccf, ec.germline_dna  from ipcm_referral_t as rf INNER JOIN ipcm_ecrf_t as ec ON regexp_replace(CAST(ec.birth_date AS VARCHAR), '-', '', 'g') =  LEFT(rf.pnr, 8)  WHERE rf.pnr='{}'".format(cfdna, pnr)
			res_data2 = fetch_sql_query('ipcmLeaderboard', sql3)
			if len(res_data2)>0:
				res_json = json.dumps(res_data2, default = json_serial)
				identifier_status = True
				json_data = json.loads(res_json)
				res_json_data = json_data[0]
			else:
				identifier_status = False
		return res_json_data, identifier_status
	except Exception as e:
		print("Build iPCM Exception", str(e))
		return [], identifier_status


# ### Fetch the genomic profile information
def build_genomic_profile_sample_details(project_name, cfdna, sample_id, capture_format):

	hospital_lookup = { "Karolinska": "KS", "Karolinska Sjukhuset": "KS", "Södersjukhuset": "SO", "St Göran": "ST" }

	sql = "SELECT study_code, study_site, dob, disease FROM genomic_profile_summary where project_name='{}' and sample_id='{}' and capture_id='{}'".format(project_name, sample_id, capture_id)
	res_data = fetch_sql_query('curation', sql)
	res_json = json.dumps(res_data, default = json_serial)

	sample_data = {}

	sample_data["identifier"] = "N/A"
	sample_data["sample_date"] = "N/A"
	sample_data["birthdate"] = "N/A"
	sample_data["hospital"] = "N/A"

	if(res_data):
		for key, val in enumerate(res_data):
			if(res_data[key]["study_code"]):
				identifier_name =  "MTBP_"+capture_format+"_"+res_data[key]["study_code"]+"_CFDNA"+cfdna
				sample_data["identifier"] = identifier_name
			if(res_data[key]["dob"] and res_data[key]["dob"] !='NA'):
				sample_data["birthdate"] = datetime.strptime(res_data[key]["dob"], "%Y-%m-%d").date().strftime("%Y-%m-%d")
			if(res_data[key]["study_site"]):
				sample_data["hospital"] = hospital_lookup[res_data[key]["study_site"]]


	sample_data["sample_date"] = "N/A"
	sample_data["seq_date"] = "N/A"
	sample_data["cancer_taxonomy"] = "N/A"
	sample_data["cancer_code"] = "N/A"
	sample_data["cancer_sub_code"] = "N/A"
	sample_data["tissue_source"] = "N/A"
	sample_data["tissue_type"] = "N/A"
	sample_data["pathology_ccf"] = "N/A"
	sample_data["bioinf_ccf"] = "N/A"
	sample_data["germline_dna"] = "N/A"

	return sample_data


# ### Fetch the sample information from biobank referral table
def build_sample_details(project_name, cfdna):

	sample_data = {}

	identifier_status = False

	sample_data["identifier"] = "N/A"
	sample_data["sample_date"] = "N/A"
	sample_data["birthdate"] = "N/A"
	sample_data["hospital"] = "N/A"

	sql = "SELECT pnr, CAST(datum AS VARCHAR) as datum, rid, tid from probio_bloodreferrals WHERE cf_dna1 like '%{}%' OR cf_dna2 like '%{}%' or cf_dna3 like '%{}%' or kommentar like '%{}%'".format(cfdna, cfdna, cfdna, cfdna)
	res_data = fetch_sql_query('referral', sql)

	if(res_data):
		sample_data["identifier"] = "MTBP_"+project_name+"_"+res_data[0]['tid']+"_CFDNA"+cfdna
		sample_data["sample_date"] = res_data[0]["datum"]
		pnr = res_data[0]["pnr"][0:8]
		sample_data["birthdate"] =  datetime.strptime(pnr, "%Y%m%d").date().strftime("%Y-%m-%d")
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

	sample_data["seq_date"] = "NA"
	sample_data["cancer_taxonomy"] = "NA"
	sample_data["cancer_code"] = "NA"
	sample_data["tissue_source"] = "NA"
	sample_data["tissue_type"] = "NA"
	sample_data["pathology_ccf"] = "NA"
	sample_data["bioinf_ccf"] = "NA"
	sample_data["germline_dna"] = "NA"

	return sample_data, identifier_status


# ### Build a QC Json 
def build_qc(root_path):
	file_path = root_path + '/qc/'

	msi_list = ''
	qc_df_data = ''

	try:
		qc_filename = file_path + list(filter(lambda x: x.endswith('.qc_overview.txt') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))[0]

		if len(qc_filename) > 0:
			qc_df_data = pd.read_csv(qc_filename, delimiter = "\t")

			column_list = ['SAMP', 'MEAN_TARGET_COVERAGE', 'contamination_%']
			column_dict = {'SAMP': 'sample', 'MEAN_TARGET_COVERAGE': 'coverage', 'contamination_%': 'contamination'}

			if 'Overall_QC' in qc_df_data.columns:
				column_list.append('Overall_QC')
				column_dict['Overall_QC'] = 'overall'

			msi_list = "N/A"
			if 'msing_score' in qc_df_data.columns:
				msi_arr = qc_df_data['msing_score'].dropna().tolist()
				msi_list = msi_arr[0] if len(msi_arr) > 0 else 'NA'

			qc_df_data = qc_df_data[column_list]
			qc_df_data = qc_df_data.rename(columns=column_dict)

			qc_df_data["contamination"] = qc_df_data["contamination"] / 100

			qc_df_data["coverage"] = qc_df_data["coverage"].round(0).astype(int)

			qc_df_data.fillna('NA', inplace=True)

			qc_json = qc_df_data.to_json(orient = 'index')

			return msi_list, qc_json

	except Exception as e:
		print("Build QC Exception", str(e))

	return msi_list, qc_df_data


# ### Fetch a Ploidy information from purecn
def build_ploidy(root_path):
	file_path = root_path + '/purecn/'

	regex = '[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+.csv'

	ploidy_list = "N/A"
	purity_list = "N/A"

	try:
		purecn_filename = file_path + list(filter(lambda x: (re.match(regex, x) ),os.listdir(file_path)))[0]

		if len(purecn_filename) > 0:
			purecn_df_data = pd.read_csv(purecn_filename, delimiter = ",")

			purecn_df_data = purecn_df_data[['Sampleid', 'Purity', 'Ploidy', 'Sex', 'Contamination']]
			ploidy_list = round(purecn_df_data['Ploidy'].tolist()[0],4)
			purity_list = round(purecn_df_data['Purity'].tolist()[0],4)

	except Exception as e:
		print("Build Polidy Exception", str(e))

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

				column_list = ['CHROM', 'START', 'END', 'REF', 'ALT']
				if(re.match(regex, i) or re.match(regex2, i)):
					reads_column_list = ['N_DP', 'N_ALT', 'N_VAF']
					column_dict = {'CHROM': 'chr', 'START': 'start', 'END': 'end', 'REF' : 'ref', 'ALT' : 'alt', 'N_DP' : 'ref_reads',  'N_ALT' : 'alt_reads', 'N_VAF' : 'vaf'}
				else:
					reads_column_list = ['T_DP', 'T_ALT', 'T_VAF']
					column_dict = {'CHROM': 'chr', 'START': 'start', 'END': 'end', 'REF' : 'ref', 'ALT' : 'alt', 'T_DP' : 'ref_reads',  'T_ALT' : 'alt_reads', 'T_VAF' : 'vaf'}

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
					smv_df_data["clonality"] = smv_df_data["clonality"].apply(lambda x: x.capitalize()) if 'clonality' in smv_df_data.columns else '' 

				smv_df_data['alt'] = smv_df_data['alt'].map(lambda x: x.lstrip('[').rstrip(']'))

				smv_df_data['strand'] = '+'

			smv_df = pd.concat([smv_df_data, smv_df])

		smv_df.fillna('NA', inplace=True)
		smv_df.reset_index(drop=True, inplace=True)
		return smv_df.to_json(orient = 'index')

	except Exception as e:
		print("Build Small Variants Exception", str(e))

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

					if(re.match(regex, i)):
						cnv_df_data['origin'] = "Germline"
					else:
						cnv_df_data['origin'] = "Somatic"

					if 'COPY_NUMBER' in cnv_df_data.columns:
						column_list.append('COPY_NUMBER')
						column_dict['COPY_NUMBER'] = 'copy_number'
					else:
						column_list.append('copy_number')
						cnv_df_data['copy_number'] = 'NA'

					cnv_df_data = cnv_df_data[column_list]
					cnv_df_data = cnv_df_data.rename(columns=column_dict)

					cnv_df_data['copy_number'] = cnv_df_data['copy_number'].round(0).astype(int,  errors='ignore')

					cnv_df_data['genes'] = cnv_df_data.genes.apply(lambda x: x.split(', '))
					cnv_df = pd.concat([cnv_df_data, cnv_df])

		cnv_df.fillna('NA', inplace=True)
		cnv_df.reset_index(drop=True, inplace=True)
		return cnv_df.to_json(orient = 'index')

	except Exception as e:
		print("Build CNV Exception", str(e))

	return cnv_df.to_json(orient = 'index')


# ### Build a SVS Json
def build_svs(root_path):
	
	file_path = root_path + '/svs/igv/'

	svs_filter = pd.DataFrame()

	try:
		svs_filename = file_path + list(filter(lambda x: (re.match('[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+-sv-annotated.txt', x) or x.endswith('_annotate_combined_SV.txt')) and not x.startswith('.') and not x.endswith('.out'),os.listdir(file_path)))[0]

		sample_list = ['germline', 'somatic', 'tumor', 'cfdna']
		svs_filter = pd.read_csv(svs_filename, delimiter = "\t")

		column_list = ['SAMPLE', 'SVTYPE', 'strand', 'vaf', 'GENE_A', 'IGV_COORD', 'GENE_B', 'CLONALITY', 'SECONDHIT']
		column_dict = {'SAMPLE': 'origin', 'SVTYPE': 'sv_type', 'GENE_A': 'gene_a' , 'IGV_COORD' : 'variant', 'GENE_B' : 'gene_b', 'CLONALITY': 'clonality', 'SECONDHIT' : 'secondhit'}

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
					svs_filter["variant_string"] = "N/A"

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


# ## Build Json from output files
def build_json(root_path, output_path, project_name, normal_cfdna, cfdna, sample_id, capture_format):

	print("--- MTBP Json Format Started ---\n")

	project_json = {}

	# Sample Information
	logging.info('--- Sample fetching started ---')
	itendifiter_status = True
	if(project_name == "IPCM" or capture_format == "iPCM"):
		sample_details_json, itendifiter_status = build_icpm_sample_details(normal_cfdna, cfdna)
	else:
		sample_details_json, itendifiter_status = build_sample_details(project_name, cfdna)

	if(sample_details_json and itendifiter_status):
		project_json["sample"] = sample_details_json
	else:
		sample_details_json = build_genomic_profile_sample_details(project_name, cfdna, sample_id, capture_format)
		project_json["sample"] = sample_details_json

	if sample_details_json["identifier"] == "N/A" :
		identifier_study_id = sample_id
	else:
		identifier_study_id = sample_details_json["identifier"].split("_")[2]
		
	logging.info('--- Sample fetching completed ---')

	# Pipeline
	project_json["pipeline"] = { "genome_reference": "hg19", "version": "1.0"}

	# QC
	logging.info('--- QC started ---')
	msi_val, qc_json = build_qc(root_path)
	project_json["qc"] =  "N/A" if qc_json == "" else json.loads(qc_json)
	logging.info('--- QC completed ---')

	# Get the Ploidy value
	purity_val, ploidy_val = build_ploidy(root_path)

	# Phenotype
	project_json["phenotypes"] = {"purity" : purity_val ,"ploidy": ploidy_val, "msi": "Low",  "hrrd":"N/A",  "tmb": "N/A"}


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

	final_json = json.dumps(project_json)

	file_name = output_path+"/MTBP_"+capture_format+"_"+identifier_study_id+"_CFDNA"+cfdna+".json"

	print("Path : ", root_path, "/", file_name)

	with open(file_name, 'w') as f:
		json.dump(project_json, f, indent=4)

	logging.info('--- Generated Json format successfuly ---\n')

	print("\n----  MTBP Json Format Completed -----\n")


# ## Main Function
def main(nfs_path, project_name, sample_id, capture_id):

	root_path = os.path.join(nfs_path,sample_id,capture_id)

	output_path = root_path+"/MTBP"

	capture_arr = capture_id.split("-")
	normal_idx = capture_arr.index("N")
	normal_cfdna = capture_arr[normal_idx+1]

	cfdna_idx = capture_arr.index("T") if 'T' in capture_arr else capture_arr.index("CFDNA")
	cfdna_id = capture_arr[cfdna_idx+1]

	#cfdna = capture_id.split("-")[4]
	cfdna = re.sub(r'[a-zA-Z]', '', cfdna_id)

	capture_format = capture_arr[0]

	log_name = output_path+"/MTBP_"+project_name+"_"+sample_id+"_CFDNA"+cfdna+".log"

	if(not os.path.exists(output_path)):
		os.mkdir(output_path)
	else:
		for f in os.listdir(output_path):
			os.remove(os.path.join(output_path, f))


	logging.basicConfig(format = '%(asctime)s  %(levelname)-10s %(name)s %(message)s', level=logging.INFO , filename=log_name, datefmt="%Y-%m-%d %H:%M:%S")
	logging.info('--- Generated Json format Started---')

	logging.info("Sample Id : {} || Capture Id : {} ".format(sample_id,capture_id))

	try:
		build_json(root_path, output_path, project_name, normal_cfdna, cfdna, sample_id, capture_format)
	except Exception as e:
		print("Main Exception", str(e))
		logging.error("Failed : {}".format(str(e)))
		logging.error('--- Generated Json format Failed ---\n')
		raise

if __name__ == "__main__":

	# Create the parser
	profile_parser = argparse.ArgumentParser(description='Generate MTBP Json')

	# Add the arguments
	profile_parser.add_argument('--path', type=str, help='Define the nfs path name',required=True)
	profile_parser.add_argument('--project', type=str, help='Define the project name (eg: PROBIO, PSFF, IPCM )',required=True)
	profile_parser.add_argument('--sample', type=str, help='Define the sample id (eg: P-00UZG3006 )',required=True)
	profile_parser.add_argument('--capture', type=str, help='Define the capture id (eg: PB-P-00UZG3006-CFDNA-0328277-KH20201014-C320201014_PB-P-00UZG3006-N-0328277-KH20201014-C320201014 )',required=True)

	args = profile_parser.parse_args()

	project_name = args.project
	sample_id = args.sample
	capture_id = args.capture
	nfs_path = args.path

	#nfs_path = "/nfs/{}/autoseq-output".format(project_name)

	if not os.path.isdir(nfs_path):
		print('The {}, path specified does not exist'.format(nfs_path))
		sys.exit()
	else:
		main(nfs_path, project_name, sample_id, capture_id)

