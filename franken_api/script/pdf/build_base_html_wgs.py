#!/usr/bin/env python
# coding: utf-8

# ### Fetch curated information into pdf 
# 
# weasyprint -e utf8 -vd --full-fonts <input.html> <ouput.pdf>
# 

import os
import pandas as pd
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
from bs4 import BeautifulSoup
import subprocess
from datetime import date, datetime
import math
from decimal import Decimal

def path():
	return os.path.dirname(os.path.realpath(__file__))


# ### Run command using subprocess
def subprocess_cmd(command):
	process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
	proc_stdout = process.communicate()[0].strip()
	for line in proc_stdout.decode().split('\n'):
		logging.info("Subprocess : {}".format(line))


# ### Read a DB information from yml file
def readConfig(section):
	filename = path()+"/config.yml"
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
		logging.error("Fetch SQL Exception : {}".format(str(error)))
	finally:
		if conn is not None:
			conn.close()

def json_serial(obj):
	if isinstance(obj, (datetime, date)):
		return obj.isoformat()
	elif isinstance(obj, Decimal):
		return float(obj)
	raise TypeError("Type %s not serializable" % type(obj))


# ### Fetch the sample information from ipcm referral table
def build_ipcm_sample_details(cfdna, normal_cfdna):
	result = ''
	identifier_status = False
	try:
		## tissue cfdna
		cfdna = re.sub(r'[a-zA-Z]', '', cfdna)
		sql = "select rf.pnr, rf.cdk from ipcm_referral_t as rf WHERE rf.rid like'%{}%' OR rf.blood like'%{}%' OR rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' OR rf.dna3 like'%{}%'".format(cfdna, cfdna, cfdna, cfdna, cfdna)
		res_data = fetch_sql_query('ipcmLeaderboard', sql)

		## normal cfdna
		sql2 = "select rf.pnr, rf.cdk from ipcm_referral_t as rf WHERE rf.rid like'%{}%' OR rf.blood like'%{}%' OR rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' OR rf.dna3 like'%{}%'".format(normal_cfdna, normal_cfdna, normal_cfdna, normal_cfdna, normal_cfdna)
		res_normal_data = fetch_sql_query('ipcmLeaderboard', sql2)

		t_pnr = res_data[0]['pnr'] if len(res_data) else ''
		t_cdk = res_data[0]['cdk'] if len(res_data) else ''

		n_pnr = res_normal_data[0]['pnr'] if len(res_normal_data) else''
		n_cdk = res_normal_data[0]['cdk'] if len(res_normal_data) else ''

		## Compare two pnr number
		if (t_pnr == n_pnr or t_pnr == ''):
			study_id = n_cdk
			dob = n_pnr[0:8]
			extra_cond = "" if study_id in ['None', ''] else " and ec.study_id='{}'".format(study_id)

			query_ecrf = "SELECT ec.study_id as identifier, to_date(rf.datum::text, 'YYYYMMDD') as sample_date, to_date(rf.date_birth::text, 'YYYYMMDD') as birthdate, get_hospital_name(ec.site_id) as hospital, 'oncotree' as cancer_taxonomy,  CASE WHEN ec.cancer_type_id != 0 THEN get_tissue_name(ec.cancer_type_id,ec.cancer_type_code) ELSE 'NA' END as tissue, ec.cancer_type_code as cancer_code, 'primary' as tissue_source, CASE WHEN (ec.cancer_type_code != '' AND ec.cancer_type_code != 'N/A' AND ec.cancer_type_code != 'NA') THEN get_cancer_info(ec.cancer_type_id, ec.cancer_type_code) ELSE get_tissue_name(ec.cancer_type_id,ec.cancer_type_code) END as disease_name, ec.cell_fraction as pathology_ccf, ec.germline_dna  from ipcm_referral_t as rf INNER JOIN ipcm_ecrf_t as ec ON to_date(rf.date_birth::text, 'YYYYMMDD') = ec.birth_date WHERE ec.birth_date=to_date('{}', 'YYYYMMDD') {} limit 1 ;".format(dob, extra_cond)
			res_ecrd_data = fetch_sql_query('ipcmLeaderboard', query_ecrf)
			if res_ecrd_data:
				result = res_ecrd_data[0]
				identifier_status = True      
		else:
			dob = t_pnr[0:8]
			query_ecrf_2="SELECT ec.study_id as identifier, to_date(rf.datum::text, 'YYYYMMDD') as sample_date, to_date(rf.date_birth::text, 'YYYYMMDD') as birthdate, get_hospital_name(ec.site_id) as hospital, 'oncotree' as cancer_taxonomy, CASE WHEN ec.cancer_type_id != 0 THEN get_tissue_name(ec.cancer_type_id,ec.cancer_type_code) ELSE 'NA' END as tissue, ec.cancer_type_code as cancer_code, 'primary' as tissue_source, CASE WHEN (ec.cancer_type_code != '' AND ec.cancer_type_code != 'N/A' AND ec.cancer_type_code != 'NA') THEN get_cancer_info(ec.cancer_type_id, ec.cancer_type_code) ELSE get_tissue_name(ec.cancer_type_id,ec.cancer_type_code) END as disease_name, ec.cell_fraction as pathology_ccf, ec.germline_dna  from ipcm_referral_t as rf INNER JOIN ipcm_ecrf_t as ec ON to_date(rf.date_birth::text, 'YYYYMMDD') = ec.birth_date WHERE ec.birth_date=to_date('{}', 'YYYYMMDD')".format(dob)
			res_ecrd_data_2 = fetch_sql_query('ipcmLeaderboard', query_ecrf_2)
			if res_ecrd_data_2:
				result = res_ecrd_data_2[0]
				identifier_status = True
		
	except Exception as e:
		logging.error("Build iPCM Sample Details SQL Exception : {}".format(str(e)))
	return result, identifier_status


# ### Fetch the sample information from biobank referral table
def build_sample_details(cfdna):

	sample_data = {}
	identifier_status = False
	try:
		sample_data["identifier"] = ""
		sample_data["sample_date"] = ""
		sample_data["birthdate"] = ""
		sample_data["hospital"] = ""
		sample_data["disease_name"] = ""

		sql = "SELECT pnr, datum, rid, tid, cdk from probio_bloodreferrals WHERE cf_dna1 like '%{}%' OR cf_dna2 like '%{}%' or cf_dna3 like '%{}%' or kommentar like '%{}%'".format(cfdna, cfdna, cfdna, cfdna)
		res_data = fetch_sql_query('referral', sql)

		if(res_data):
			tid = res_data[0]["tid"]
			cdk_1 = res_data[0]["cdk"]
			cdk =  cdk_1.strip() if cdk_1 != None else cdk_1
			identifier_status = True
			sample_data["identifier"] =  cdk if cdk != "" else cdk
			sample_data["sample_date"] = res_data[0]["datum"]
			pnr_1 = res_data[0]["pnr"]
			pnr = pnr_1[0:8] if pnr_1 != None else pnr_1

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
				sample_data["identifier"] = subject_id
				identifier_status = True
				if(glb_data_2):
					sample_data["birthdate"] = glb_data_2[0]["dob"]
					sample_data["hospital"] = glb_data_2[0]["site_name"]

	except Exception as e:
		logging.error("Build Sample Details Exception : {}".format(str(e)))
	return sample_data, identifier_status


def build_genomic_profile_sample_details(project_name, sample_id, capture_id):

	hospital_lookup = { "Karolinska": "KS", "Karolinska Sjukhuset": "KS", "Södersjukhuset": "SO", "St Göran": "ST" }

	sql = "SELECT study_code, study_site, dob, disease FROM genomic_profile_summary where project_name='{}' and sample_id='{}' and capture_id='{}'".format(project_name, sample_id, capture_id)
	res_data = fetch_sql_query('curation', sql)
	res_json = json.dumps(res_data, default = json_serial)
	sample_data = {}

	capture_arr = capture_id.split("-")
	seq_date_str = re.findall(r'\d+', capture_arr[5])[0][:8]
	seq_date = datetime.strptime(seq_date_str, "%Y%m%d").date().strftime("%Y-%m-%d")

	sample_data["identifier"] = "NA"
	sample_data["sample_date"] = seq_date
	sample_data["birthdate"] = "NA"
	sample_data["hospital"] = "NA"
	sample_data["disease_name"] = "NA"
	
	if(res_data):
		for key, val in enumerate(res_data):
			disease = res_data[key]["disease"]
			study_code = res_data[key]["study_code"]
			dob = res_data[key]["dob"]
			study_site = res_data[key]["study_site"]
		
			if(study_code):
				sample_data["identifier"] = study_code
			if(dob and dob !='NA'):
				pattern = re.compile(r'^\d{4}\d{2}\d{2}$')  # '%Y%m%d'
				if pattern.match(dob):
					dob_format = "%Y%m%d"
				else:
					dob_format = "%Y-%m-%d"

				sample_data["birthdate"] = datetime.strptime(dob.strip(), dob_format).date().strftime("%Y-%m-%d")

			if(study_site and study_site != 'NA'):
				sample_data["hospital"] = hospital_lookup[study_site] if (study_site in hospital_lookup) else study_site

			if(disease != ""):
				sample_data["disease_name"] = disease

	return sample_data

# ### Build a patient, specimen & assay and genome wide information
def build_basic_html(project_name, sample_id, capture_id, study_id, disease_name):

	patient_info_html = ''
	specimen_assay_html = ''
	genome_wide_html = ''
	summary_txt = ''
	sample_txt = ''

	study_code = 'None'
	disease = 'None'
	tumörcellsandel = 'None'
	tumörmutationsbörda = 'None'
	msi_status = 'None'
	potentiellt_ärftliga = 'None'

	try:
		sql = "SELECT study_code, study_site, dob, disease, specimen_assay, genome_wide, summary_txt from genomic_profile_summary WHERE sample_id='{}' and capture_id='{}' and project_name='{}' order by id desc limit 1".format(sample_id, capture_id, project_name)
		res_data = fetch_sql_query('curation', sql)

		if(res_data):

			study_code = res_data[0]['study_code'] if study_id =='' else study_id
			study_site = res_data[0]['study_site']
			dob = res_data[0]['dob']
			disease = res_data[0]['disease'] if disease_name =='' else disease_name

			summary_txt = res_data[0]['summary_txt']

			patient_info_html += '<tr><th>STUDY ID</th><td>{}</td></tr>'.format(study_code)
			patient_info_html += '<tr><th>PERSONAL NUMBER</th><td>{}</td></tr>'.format(dob)
			patient_info_html += '<tr><th>DISEASE</th><td>{}</td></tr>'.format(disease)
			patient_info_html += '<tr><th>HOSPITAL</th><td>{}</td></tr>'.format(study_site)

			specimen_assay = res_data[0]['specimen_assay'].replace("\'", "\"")

			specimen_assay_json = json.loads(specimen_assay)

			for i in specimen_assay_json:

				quality =''
				quality_str = ''
				if 'quality' in specimen_assay_json[i]:
					quality = specimen_assay_json[i]['quality'] if isinstance(specimen_assay_json[i]['quality'], str) else ''
					if quality =='Pass':
						quality_str= '<span class="pass"></span>'
					elif quality == 'Fail':
						quality_str= '<span class="fail"></span>'
				
				specimen_assay = specimen_assay_json[i]["assay"] if 'assay' in specimen_assay_json[i] else ''

				specimen_assay_html += '<tr>'
				specimen_assay_html += '<th>'+specimen_assay_json[i]["specimen"]+'</th>'
				specimen_assay_html += '<td>'+specimen_assay_json[i]["analyte"]+'</td>'
				specimen_assay_html += '<td>'+specimen_assay+'</td>'
				specimen_assay_html += '<td class="quality">'+quality_str+'</td>'
				specimen_assay_html += '</tr>'

			genome_wide = res_data[0]['genome_wide'].replace("\'", "\"")

			genome_wide_json = json.loads(genome_wide)

			for j in genome_wide_json:
				if 'CHIP' not in genome_wide_json[j]["title"]:
					result_data = genome_wide_json[j]["result"] if 'result' in genome_wide_json[j] and  isinstance(genome_wide_json[j]["result"], str) else '-'
				else:
					result_data = genome_wide_json[j]["result"] if 'result' in genome_wide_json[j] and  isinstance(genome_wide_json[j]["result"], str) else ''

				assessment =''
				if 'assessment' in genome_wide_json[j]:
					assessment = 'Yes' if genome_wide_json[j]["assessment"] == 'Possible' else ( 'No' if genome_wide_json[j]["assessment"] == 'Not possible' else '')

				genome_wide_html += '<tr>'
				genome_wide_html += '<th>'+genome_wide_json[j]["title"]+'</th>'
				genome_wide_html += '<td>'+result_data+'</td>'
				genome_wide_html += '<td>'+assessment+'</td>'


				if ('assessment' in genome_wide_json[j]) and (genome_wide_json[j]["title"] == 'OTHER GENOMIC PHENOTYPE' and genome_wide_json[j]["assessment"] == 'Possible' and result_data == 'Yes'):
					genome_wide_html += '<td>'+genome_wide_json[j]["categories"]+'</td>'
				else:
					genome_wide_html += '<td>'+genome_wide_json[j]["comment"]+'</td>'
				genome_wide_html += '</tr>'

				### Generate a clinical report txt 
				title = genome_wide_json[j]['title']
				if (title == "FRACTION OF CANCER DNA"):
					tumörcellsandel = result_data
				elif (title == "TUMOR MUTATIONAL BURDEN"):
					tumörmutationsbörda = result_data
				elif (title == "MSI STATUS"):
					msi_status = result_data
				elif (title == "PATHOGENIC GERMLINE VARIANTS"):
					potentiellt_ärftliga = result_data 
		else:
			patient_info_html += '<tr><th>STUDY ID</th><td>-</td></tr>'
			patient_info_html += '<tr><th>PERSONAL NUMBER</th><td>-</td></tr>'
			patient_info_html += '<tr><th>DISEASE</th><td>-</td></tr>'
			patient_info_html += '<tr><th>HOSPITAL</th><td>-</td></tr>'
	
	except Exception as e:
		logging.error("Build Base HTML Exception : {}".format(str(e)))

	### Generate a clinical report txt 
	sample_txt = 'Genomisk karaktärisering av solid cancer - Implementation of Personalized Cancer Medicine (iPCM)\n'
	sample_txt += 'STUDIENUMMER\t{}\nTUMÖRTYP\t{}\n'.format(study_code, disease)
	sample_txt += '\nSEKVENSERAT MATERIAL:\n'
	sample_txt += 'TUMÖRCELLSANDEL: \t{}\n'.format(tumörcellsandel)
	sample_txt += 'TUMÖRMUTATIONSBÖRDA (TMB): \t{}\n'.format(tumörmutationsbörda)
	sample_txt += 'MMR/MSI-status: \t{}\n'.format(msi_status)
	sample_txt += 'KONSTITUTIONELLA (POTENTIELLT ÄRFTLIGA) VARIANTER: \t{}\n\n'.format(potentiellt_ärftliga)

	return patient_info_html, specimen_assay_html, genome_wide_html, summary_txt, sample_txt, study_code


# ### Build a Small variant Json (Somatic & Germline)
def build_small_variants(root_path):

	smt_variant_html = ''
	smt_variant_txt = ''
	try:
		file_path = root_path + '/'
		smv_file_list = list(filter(lambda x: x.endswith('-igvnav-input.txt') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))
		regex = '^(?:(?!-(CFDNA|T)).)*igvnav-input.txt$'
		regex2 = '(.*)-(CFDNA|T)-(\w.*)(germline-igvnav-input).*txt$'

		for i in smv_file_list:
			smv_filename = file_path + i
			smv_df_data = pd.read_csv(smv_filename, delimiter = "\t")

			if 'CALL' in smv_df_data.columns:
				smv_df_call_filter = smv_df_data.loc[(smv_df_data['CALL'] == "S") | (smv_df_data['CALL'] == "G")]

				if 'include_variant_report_pdf' in smv_df_call_filter.columns:
					smv_df_call_filter = smv_df_call_filter.loc[(smv_df_call_filter['include_variant_report_pdf'] == "True") | (smv_df_call_filter['include_variant_report_pdf'] == True)]

				column_list = ['CHROM', 'START', 'END', 'REF', 'ALT', 'GENE', 'CONSEQUENCE', 'HGVSp']
				column_dict = {'CHROM': 'chr', 'START': 'start', 'END': 'end', 'REF' : 'ref', 'ALT' : 'alt'}

				clonality_boolean = False

				if 'CLONALITY' in smv_df_call_filter.columns:
					column_list.append('CLONALITY')
					column_dict['CLONALITY'] = 'clonality'
					clonality_boolean = True

				if 'SECONDHIT' in smv_df_call_filter.columns:
					column_list.append('SECONDHIT')
					column_dict['SECONDHIT'] = 'second_hit'

				if 'HGVSp_org' in smv_df_call_filter.columns:
					column_list.append('HGVSp_org')
					column_dict['HGVSp_org'] = 'HGVSp_org'

				if 'TRANSCRIPT' in smv_df_call_filter.columns:
					column_list.append('TRANSCRIPT')
					column_dict['TRANSCRIPT'] = 'TRANSCRIPT'

				if 'ASSESSMENT' in smv_df_call_filter.columns:
					column_list.append('ASSESSMENT')
					column_dict['ASSESSMENT'] = 'ASSESSMENT'

				if 'zygosity' in smv_df_call_filter.columns:
					column_list.append('zygosity')
					column_dict['zygosity'] = 'zygosity'

				if 'RSID' in smv_df_call_filter.columns:
					column_list.append('RSID')
					column_dict['RSID'] = 'RSID'

				smv_df_call_filter = smv_df_call_filter[column_list]
				smv_df_filter_data = smv_df_call_filter.rename(columns=column_dict)

				sort_arr = ["GENE"]

				if clonality_boolean : 
					sort_arr = ["clonality", "GENE"]

				smv_df_filter_data = smv_df_filter_data.sort_values(sort_arr, ascending = True)

				smv_df_filter_data.fillna('-', inplace=True)

				source_type = ''

				if(re.match(regex, i) or re.match(regex2, i)):
					source_type ='germline'
				else:
					source_type ='somatic'

				if len(smv_df_filter_data):
					smt_variant_txt += '\nGENE\tSource\tVariant-details\tConsequence\tClonality\tSecond-Hit\tHGVSP\tONCOKB-TIER\n'

					for index, row in smv_df_filter_data.iterrows():

						variant_det = "chr"+str(row['chr'])+":"+str(row['end'])+'; '+row['ref']+'>'+row['alt'].lstrip('[').rstrip(']')
						if source_type == 'somatic':
							clonality = row['clonality'] if 'clonality' in row else '-'
						else:
							clonality = row['zygosity'].upper() if 'zygosity' in row else '-'
						# clonality = row['clonality'] if 'clonality' in row else '-'
						transcript = row['TRANSCRIPT'] if 'TRANSCRIPT' in row else '-'
						assessment = row['ASSESSMENT'] if 'ASSESSMENT' in row else '-'

						if 'RSID' in row:
							rs_id_arr = ''
							if re.findall("'\s*([^']*?)\s*'", row['RSID']):
								rsId_arr = eval(row['RSID'])
								for rsId_str in rsId_arr:
									if 'rs' in rsId_str:
										rs_id_arr += rsId_str+','

							rsID = rs_id_arr[:-1]
						else: 
							rsID = '' 

						HGVSp_org_rx = row['HGVSp_org'] if 'HGVSp_org' in row else '-'
						HGVSp_org = HGVSp_org_rx.split("p.")[1] if 'p.' in HGVSp_org_rx else HGVSp_org_rx
						second_hit = str(row['second_hit']) if 'second_hit' in row else '-'
						consequence = row['CONSEQUENCE']
						HGVSp_one_cod = row['HGVSp']
										
						smt_variant_html += '<tr>'
						smt_variant_html +='<td>'+row['GENE']+'</td>'
						smt_variant_html +='<td>'+source_type+'</td>'
						smt_variant_html +='<td class="sm-var-dets"><p>'+variant_det+'</p></td>'
						smt_variant_html +='<td>'+consequence+'</td>'
						smt_variant_html +='<td>'+clonality+'</td>'
						smt_variant_html +='<td>'+second_hit+'</td>'
						smt_variant_html +='<td>'+assessment+'</td>'
						smt_variant_html +='<td>'+rsID+'</td>'
						smt_variant_html +='<td>'+transcript+'</td>'
						smt_variant_html +='<td>'+HGVSp_org+'</td>'
						smt_variant_html += '</tr>'

						smt_variant_txt += row['GENE']+'\t'
						smt_variant_txt += source_type+'\t'
						smt_variant_txt += variant_det+'\t'
						smt_variant_txt += consequence+'\t'
						smt_variant_txt += clonality+'\t'
						smt_variant_txt += second_hit+'\t'
						smt_variant_txt += HGVSp_one_cod+'\t'
						smt_variant_txt += '-\t'
						smt_variant_txt += '\n'

	
	except Exception as e:	
		smt_variant_html = ''
		smt_variant_txt = ''
		logging.error("Small Variants Exception : {}".format(str(e)))

	return smt_variant_html, smt_variant_txt


# ### Build a CNV Json (Somatic & Germline)
def build_cnv(root_path):
	cnv_html = ''
	cnv_txt = '' 

	try:
		cnv_df = pd.DataFrame()
		file_path = root_path + '/cnv/'
		cnv_file_list = list(filter(lambda x: x.endswith('_curated.cns') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))


		if len(cnv_file_list) > 0:

			regex = '^(?:(?!-(CFDNA|T)).)*_curated.cns$'
			for i in cnv_file_list:
				cnv_filename = file_path + i
				cnv_df_data = pd.read_csv(cnv_filename, delimiter = "\t")

				column_list = ['chromosome', 'start', 'end', 'gene', 'ASSESSMENT']
				column_dict = {'chromosome': 'chr'}

				if 'ASSESSMENT' in cnv_df_data.columns:
					cnv_txt += 'GENE\tSource\tVariant-details\tAssessment\n'
					cnv_df_data = cnv_df_data.loc[(cnv_df_data['ASSESSMENT'].notnull())]

					if 'include_variant_report_pdf' in cnv_df_data.columns:
						cnv_df_data = cnv_df_data.loc[(cnv_df_data['include_variant_report_pdf'] == "True") | (cnv_df_data['include_variant_report_pdf'] == True)]

					if(re.match(regex, i)):
						source_type = "germline"
					else:
						source_type = "somatic"

					if 'COPY_NUMBER' in cnv_df_data.columns:
						column_list.append('COPY_NUMBER')
						column_dict['COPY_NUMBER'] = 'copy_number'
					else:
						column_list.append('copy_number')
						cnv_df_data['copy_number'] = '-'

					cnv_df_data = cnv_df_data[column_list]
					cnv_df_data = cnv_df_data.rename(columns=column_dict)

					cnv_df_data.fillna('-', inplace=True)

					for index, row in cnv_df_data.iterrows():

						variant_det = "chr"+str(row['chr'])+":"+str(row['start'])+"-"+str(row['end'])

						# gene_det = re.sub(r'\s', ' ', row['gene'])
						gene_det = " ".join(row['gene'].split())

						copy_number = row['copy_number'] if isinstance(row['copy_number'],str) else int(row['copy_number'])

						cnv_html += '<tr>'
						cnv_html +='<td class="cp-gene">'+gene_det+'</td>'
						cnv_html +='<td>'+source_type+'</td>'
						cnv_html +='<td>'+variant_det+'</td>'
						cnv_html +='<td>'+row['ASSESSMENT']+'</td>'
						cnv_html +='<td>'+str(copy_number)+'</td>'
						cnv_html += '</tr>'

						cnv_txt +=gene_det+'\t'
						cnv_txt +=source_type+'\t'
						cnv_txt +=variant_det+'\t'
						cnv_txt +=row['ASSESSMENT']+'\t'
						cnv_txt += '\n'

	except Exception as e:	
		cnv_html = ''
		cnv_txt = ''
		logging.error("CNV Exception : {}".format(str(e)))

	return cnv_html, cnv_txt


# ### Build a SVS Json
def build_svs(root_path):
	file_path = root_path + '/svs/igv/'

	svs_html = ''
	svs_txt = ''
	try:
		svs_filename = file_path + list(filter(lambda x: (re.match('[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+-sv-annotated.txt', x) or x.endswith('_annotate_combined_SV.txt')) and not x.startswith('.') and not x.endswith('.out'),os.listdir(file_path)))[0]

		sample_list = ['germline', 'somatic', 'tumor', 'cfdna']
		svs_filter = pd.read_csv(svs_filename, delimiter = "\t")

		column_list = ['CHROM_A', 'START_A', 'END_A', 'GENE_A', 'CHROM_B', 'START_B', 'END_B', 'GENE_B', 'SAMPLE', 'IGV_COORD', 'CLONALITY', 'SECONDHIT', 'CALL']
		column_dict = {}

		if 'CALL' in svs_filter.columns:
			svs_filter = svs_filter.loc[(svs_filter['CALL'] == True ) | (svs_filter['CALL'] == 'true')]

			if 'include_variant_report_pdf' in svs_filter.columns:
				svs_filter = svs_filter.loc[(svs_filter['include_variant_report_pdf'] == "True") | (svs_filter['include_variant_report_pdf'] == True)]

			if not svs_filter.empty:
				svs_txt += 'GENE_A,GENE_B\tSource\tVariant-details\tConsequence\tClonality\tSecond-Hit\n'
				svs_filter = svs_filter[svs_filter['SAMPLE'].isin(sample_list)]

				if 'CONSEQUENCE' in svs_filter.columns:
					column_list.append('CONSEQUENCE')
					column_dict['CONSEQUENCE'] = 'consequence'
				else:
					svs_filter["consequence"] = '-'

				svs_filter = svs_filter[column_list]
				svs_filter = svs_filter.rename(columns=column_dict)

				svs_filter.fillna('-', inplace=True)

				for index, row in svs_filter.iterrows():

					chr_b = ' | chr'+str(row['CHROM_B'])+":" if row['CHROM_B'] != '-' else ' | '
					start_b = str(int(row['START_B']))+"-" if row['START_B'] != '-' else ''
					end_b = str(int(row['END_B'])) if row['END_B'] != '-' else ''

					if start_b != '':
						variant_det = "chr"+str(row['CHROM_A'])+":"+str(int(row['START_A']))+"-"+str(int(row['END_A']))+chr_b+start_b+end_b
					else:
						variant_det = "chr"+str(row['CHROM_A'])+":"+str(int(row['START_A']))+"-"+str(int(row['END_A']))

					gene_det = row['GENE_A']+","+ row['GENE_B']
					sample_type = row['SAMPLE']
					consequence = row['consequence']
					clonality = row['CLONALITY']
					secondHit = str(row['SECONDHIT'])

					svs_html += '<tr>'
					svs_html +='<td>'+row['GENE_A']+'</td>'
					svs_html +='<td>'+row['GENE_B']+'</td>'
					svs_html +='<td>'+sample_type+'</td>'
					svs_html +='<td>'+variant_det+'</td>'
					svs_html +='<td>'+consequence+'</td>'
					svs_html +='<td>'+clonality+'</td>'
					svs_html +='<td>'+secondHit+'</td>'
					svs_html += '</tr>'

					svs_txt +=gene_det+'\t'
					svs_txt +=sample_type+'\t'
					svs_txt +=variant_det+'\t'
					svs_txt +=consequence+'\t'
					svs_txt +=clonality+'\t'
					svs_txt +=secondHit+'\t'
					svs_txt += '\n'

	except Exception as e:	
		svs_html = ''
		svs_txt = ''
		logging.error("SVS Exception : {}".format(str(e)))

	return svs_html, svs_txt


# ### Build a Technical Info from QC
def build_tech_val_QC(root_path, project_name, capture_id):
	file_path = root_path + '/qc/'

	tech_html = ''

	tech_header_html = ''

	try:

		capture_arr = capture_id.split("_")
		specimen =  'cfDNA' if 'CFDNA' in capture_arr[0] else ( 'FFPE' if 'T' in capture_arr[0] else 'Tumor')
		specimen_2 = 'Whole blood' if 'N' in capture_arr[1] else ''

		tech_header_html += '<tr>'
		tech_header_html +='<th>NGS</th>'
		tech_header_html +='<th>'+specimen+'</th>'
		tech_header_html +='<th>'+specimen_2+'</th>'
		tech_header_html +='</tr>'

		regex = '^(?:(?!-(CFDNA|T)).)*$'

		tech_html += '<tr>'
		tech_html +='<th>LIBRARY PREP KIT</th>'
		tech_html +='<td>Kapa Hyperprep</td>'
		tech_html +='<td>Kapa Hyperprep</td>'
		tech_html +='</tr>'

		tech_html +='<tr>'
		tech_html +='<th>HYBRIDISATION CAPTURE</th>'
		if project_name == 'PROBIO':
			tech_html +='<td>Prostate cancer comprehensive panel v2</td>'
			tech_html +='<td>Prostate cancer comprehensive panel v2</td>'
		else:
			tech_html +='<td>ProBio comprehensive panel v2</td>'
			tech_html +='<td>ProBio comprehensive panel v2</td>'

		tech_html +='</tr>'

		qc_filename = file_path + list(filter(lambda x: x.endswith('.qc_overview.txt') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))[0]

		if len(qc_filename) > 0:
			qc_df_data = pd.read_csv(qc_filename, delimiter = "\t")

			column_list = ['SAMP', 'MEAN_TARGET_COVERAGE', 'READ_PAIRS_EXAMINED', 'PERCENT_DUPLICATION']
			column_dict = {'SAMP': 'sample', 'MEAN_TARGET_COVERAGE': 'coverage', 'READ_PAIRS_EXAMINED': 'read_pair', 'PERCENT_DUPLICATION':'duplication'}

			qc_df_data = qc_df_data[column_list]
			qc_df_data = qc_df_data.rename(columns=column_dict)

			qc_df_data.fillna('-', inplace=True)

			qc_df_json = qc_df_data.to_dict()

			df_len = len(qc_df_data.index)

			tech_html += '<tr>'
			tech_html +='<th>TOTAL NUMBER OF READS PAIRS</th>'
			tech_html +='<td>'+str(qc_df_json["read_pair"][0])+'</td>'
			if(df_len > 1):
				tech_html +='<td>'+str(qc_df_json["read_pair"][1])+'</td>'
			else:
				tech_html +='<td>-</td>'

			tech_html +='</tr>'
			tech_html +='<tr>'
			tech_html +='<th>MEAN TARGET COVERAGE</th>'
			tech_html +='<td>'+str(math.ceil(int(qc_df_json["coverage"][0])))+'</td>'
			if(df_len > 1):
				tech_html +='<td>'+str(math.ceil(int(qc_df_json["coverage"][1])))+'</td>'
			else:
				tech_html +='<td>-</td>'

			tech_html +='</tr>'
			tech_html +='<tr>'
			tech_html +='<th>FRACTION DUPLICATES</th>'
			tech_html +='<td>'+str(round(qc_df_json["duplication"][0],2))+'</td>'
			if(df_len > 1):
				tech_html +='<td>'+str(round(qc_df_json["duplication"][1],2))+'</td>'
			else:
				tech_html +='<td>-</td>'
			tech_html +='</tr>'


	except Exception as e:
		logging.error("TECHNICAL VALIDATION Exception : {}".format(str(e)))


	return tech_html, tech_header_html


# ### Build HTML from output files
def build_html(root_path, html_root_path, file_name, project_name, cfdna, capture_format, base_html_path, sample_id,capture_id, study_id, disease_name, header_txt, sub_header_txt):
	logging.info("--- MTBP Json Format Started ---")
	logging.info("Path : {} / {}".format(root_path,file_name))

	report_txt = ''

	# Patient, Specimen & assay and genome wide information 
	logging.info('--- Patient, Specimen & assay and genome wide information started ---')
	pat_html, specimen_html, genome_wide_html, summary_txt, sample_details_txt, study_id = build_basic_html(project_name, sample_id,capture_id, study_id, disease_name)

	report_txt += sample_details_txt

	# # Small Variant (Somatic & Germline)
	logging.info('--- Small Variant started ---')
	small_variant_html, small_variant_txt = build_small_variants(root_path)
	report_txt += small_variant_txt
	if(small_variant_html == ""):
		small_variant_html = '<tr><td class="no-data" colspan="9">No relevant variants detected</td></tr>'
	logging.info('--- Small Variant completed ---')

	# # SVS
	logging.info('--- SVS started ---')
	svs_html, svs_txt = build_svs(root_path)
	report_txt += "\n"+svs_txt
	if(svs_html == ""):
		svs_html = '<tr><td class="no-data" colspan="7">No relevant variants detected</td></tr>'
	logging.info('--- SVS completed ---')

	# # CNVs
	logging.info('--- CNVs started ---')
	cnv_html, cnv_txt = build_cnv(root_path)
	report_txt += "\n"+cnv_txt
	if(cnv_html == ""):
		cnv_html = '<tr><td class="no-data" colspan="5">No relevant variants detected</td></tr>'
	logging.info('--- CNVs completed ---')

	new_text = ""
	# # Read a base html and replace the curated text based on ids
	with open(base_html_path, 'r') as f:
		contents = f.read()
		soup = BeautifulSoup(contents, "html.parser")

		for tag in soup.find_all(title="custom css"):
			tag['href'] = os.path.join(html_root_path, tag['href'])

		for tag in soup.find_all(id='header_ul_data'):
			header_html = header_txt
			tag.string.replace_with(header_html)

		for tag in soup.find_all(id='header_ul_data2'):
			sub_header_html = sub_header_txt
			tag.string.replace_with(sub_header_html)

		for images in soup.find_all(id='ki_logo_img'):
			img_path = html_root_path +'/static/img/logo/KS_SE.svg'
			images['src'] = images['src'].replace("#", img_path)

		for images in soup.find_all(id='logo_img'):
			img_path = html_root_path +'/static/img/logo/KS_logo.svg'
			images['src'] = images['src'].replace("#", img_path)


		for tag in soup.find_all(id='patient_info_table_data'):
			tag.string.replace_with(pat_html)
			 
		for tag in soup.find_all(id='specimen_assay_table_data'):
			tag.string.replace_with(specimen_html)

		for tag in soup.find_all(id='genome_wide_table_data'):
			tag.string.replace_with(genome_wide_html)

		for tag in soup.find_all(id='somatic_table_data'):
			tag.string.replace_with(small_variant_html)

		for tag in soup.find_all(id='svs_table_data'):
			tag.string.replace_with(svs_html)

		for tag in soup.find_all(id='cnv_table_data'):
			tag.string.replace_with(cnv_html)

		for tag in soup.find_all(id='summary_notes_data'):
			tag.string.replace_with(summary_txt)          
		
		new_text = soup.prettify(formatter=None)


	# Create a new html based on the base template
	with open(file_name, "w", encoding = 'utf-8') as file:
		file.write(str(new_text))

	# # create a clinical report txt file and store information
	output_path = root_path+"/pdf"
	study_id = sample_id if study_id == None else study_id
	file_name = output_path+"/autoseq_clinical_report_"+study_id+".txt"

	report_txt += "\n"
	report_txt += "\nESCAT TIER I - IV: A framework to rank genomic alterations as targets for cancer precision medicine: the ESMO Scale for Clinical Actionability of molecular Targets (ESCAT). Mateo J, et al. Ann Oncol. 2018 Sep 1;29(9):1895-1902. doi: 10.1093/annonc/mdy263. PMID: 30137196; PMCID: PMC6158764."
	report_txt += "\nAMP TIER I - IV: Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists. Li et al., J Mol Diagn. 2017 Jan;19(1):4-23. doi: 10.1016/j.jmoldx.2016.10.002. PMID: 27993330; PMCID: PMC5707196."
	
	final_report_txt = report_txt

	with open(file_name, 'w') as f:
		f.write(final_report_txt)

	return 1


# ### Main Function
def main(nfs_path, project_name, sample_id, capture_id):

	html_root_path = os.path.dirname(sys.argv[0])

	### Serve Path 
	base_html_path = html_root_path+'/base_wgs.html'
	root_path = os.path.join(nfs_path,sample_id,capture_id)

	output_path = root_path+"/pdf"
	capture_arr = capture_id.split("-")


	# cfdna = capture_arr[4]
	normal_idx = capture_arr.index("N")
	normal_cfdna = capture_arr[normal_idx+1]

	cfdna_idx = capture_arr.index("T") if 'T' in capture_arr else capture_arr.index("CFDNA")
	cfdna = capture_arr[cfdna_idx+1]

	capture_format = capture_arr[0]

	file_name = output_path+"/base_"+project_name+"_"+sample_id+"_"+cfdna+".html"
	pdf_file_name = output_path+"/pdf_"+project_name+"_"+sample_id+"_"+cfdna+".pdf"
	log_name = output_path+"/log_"+project_name+"_"+sample_id+"_"+cfdna+".log"

	if(not os.path.exists(output_path)):
		os.mkdir(output_path)
	else:
		for f in os.listdir(output_path):
			os.remove(os.path.join(output_path, f))

	specimen =  'Fresh Frozen'

	# Sample Information
	sample_date = '-'
	study_id = ''
	disease_name = ''

	sample_details_json = build_genomic_profile_sample_details(project_name, sample_id, capture_id)
	sample_date = sample_details_json["sample_date"]
	study_id = sample_details_json["identifier"]
	disease_name = sample_details_json["disease_name"]

	report_date = date.today().strftime('%Y-%m-%d')

	if 'PROBIO' in project_name:
		epm_dnr_data_html = "2016/101-32"
	else:
		epm_dnr_data_html = "2021-00135"

	# Build Header Information 
	logging.info('--- Header information started ---')
	header_txt = ''
	header_txt = '<li><b>SPECIMENS</b><span>'+specimen+'</span></li>'
	header_txt += '<li><b>EPM DNR</b><span>'+epm_dnr_data_html+'</span></li>'
	header_txt += '<li><b>REPORT DATE</b><span>'+report_date+'</span></li>'
	header_txt += '<li><b>RECIEVED DATE</b><span>'+str(sample_date)+'</span></li>'


	sub_header_txt = ''
	sub_header_txt = '<li><b>REPORT DATE</b><span>'+report_date+'</span></li>'

	logging.basicConfig(format = '%(asctime)s  %(levelname)-10s %(name)s %(message)s', level=logging.INFO , filename=log_name, datefmt =  "%Y-%m-%d %H:%M:%S")
	logging.info('--- Generated Json format Started---')

	logging.info("Sample Id : {} || Capture Id : {} || Outpue File Name : {} ".format(sample_id,capture_id, file_name))

	try:
		
		html_result = build_html(root_path, html_root_path, file_name, project_name, cfdna, capture_format, base_html_path, sample_id, capture_id, study_id, disease_name, header_txt, sub_header_txt)
		if html_result:
			cmd ="weasyprint -e utf8 -vd {} {}".format(file_name, pdf_file_name)
			subprocess_cmd(cmd)
			logging.info("---- PDF Generated ----")

	except Exception as e:
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

	if not os.path.isdir(nfs_path):
		print('The {}, path specified does not exist'.format(nfs_path))
		sys.exit()
	else:
		main(nfs_path, project_name, sample_id, capture_id)

