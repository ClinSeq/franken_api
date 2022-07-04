#!/usr/bin/env python
# coding: utf-8

# ### Fetch curated information into pdf 
# 
# * wkhtmltopdf 'output.html' --header-html '/nfs/IPCM/pdf/layout/header.html' --footer-line --footer-html '/nfs/IPCM/pdf/layout/footer.html' 'appendix.html' --header-html '/nfs/IPCM/pdf/layout/appendix_header.html' --footer-line  --footer-html '/nfs/IPCM/pdf/layout/footer.html' iPCM_test.pdf
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

def path():
	return os.path.dirname(os.path.realpath(__file__))


# ### Run command using subprocess

def subprocess_cmd(command):
	process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
	proc_stdout = process.communicate()[0].strip()
	for line in proc_stdout.decode().split('\n'):
		print (line)


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
		print(error)
	finally:
		if conn is not None:
			conn.close()


# ### Fetch the sample information from ipcm referral table

def build_icpm_sample_details(cfdna):
	sql = "SELECT ec.study_id as identifier, to_date(rf.datum::text, 'YYYYMMDD') as sample_date, to_date(rf.date_birth::text, 'YYYYMMDD') as birthdate, get_hospital_name(rf.site_id) as hospital, 'oncotree' as cancer_taxonomy,  ec.cancer_type_code as cancer_code, 'primary' as tissue_source, get_tissue_name(ec.cancer_type_id, ec.cancer_type_code) as tissue_type, ec.cell_fraction as pathology_ccf, ec.germline_dna  from ipcm_referral_t as rf INNER JOIN ipcm_ecrf_t as ec ON CAST(rf.cdk as VARCHAR) = ec.study_id WHERE rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' or rf.dna3 like'%{}%'".format(cfdna, cfdna, cfdna)
	res_data = fetch_sql_query('ipcmLeaderboard', sql)
	if res_data:
		res_data = res_data[0] 
	return res_data


# ### Fetch the sample information from biobank referral table

def build_sample_details(cfdna):
	
	sample_data = {}
	
	sample_data["identifier"] = ""
	sample_data["sample_date"] = ""
	sample_data["birthdate"] = ""
	sample_data["hospital"] = ""

	sql = "SELECT pnr, datum, rid, tid from probio_bloodreferrals WHERE cf_dna1 like '%{}%' OR cf_dna2 like '%{}%' or cf_dna3 like '%{}%' or kommentar like '%{}%'".format(cfdna, cfdna, cfdna, cfdna)
	res_data = fetch_sql_query('referral', sql)

	if(res_data):
		sample_data["identifier"] = res_data[0]['tid']
		sample_data["sample_date"] = res_data[0]["datum"]
		pnr = res_data[0]["pnr"][0:8]
		
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
			
	
	return sample_data


# ### Build a patient, specimen & assay and genome wide information

def build_basic_html(sample_id, capture_id):
	
	patient_info_html = ''
	specimen_assay_html = ''
	genome_wide_html = ''
	summary_txt = ''
	
	sql = "SELECT study_code, study_site, dob, disease, specimen_assay, genome_wide, summary_txt from genomic_profile_summary WHERE sample_id='{}' and capture_id='{}'".format(sample_id, capture_id)
	res_data = fetch_sql_query('curation', sql)
	
	if(res_data):
	
		study_code = res_data[0]['study_code']
		study_site = res_data[0]['study_site']
		dob = res_data[0]['dob']
		disease = res_data[0]['disease']

		# summary_txt = '<tr><td>'+res_data[0]['summary_txt']+'</td></tr>'
		summary_txt = res_data[0]['summary_txt']


		patient_info_html += '<tr><th>STUDY ID</th><td>{}</td></tr>'.format(study_code)
		patient_info_html += '<tr><th>DATE OF BIRTH</th><td>{}</td></tr>'.format(dob)
		patient_info_html += '<tr><th>DISEASE</th><td>{}</td></tr>'.format(disease)
		patient_info_html += '<tr><th>HOSPITAL</th><td>{}</td></tr>'.format(study_site)

		specimen_assay = res_data[0]['specimen_assay'].replace("\'", "\"")

		specimen_assay_json = json.loads(specimen_assay)

		for i in specimen_assay_json:

			specimen_assay_html += '<tr>'
			specimen_assay_html += '<td>'+specimen_assay_json[i]["specimen"]+'</td>'
			specimen_assay_html += '<td>'+specimen_assay_json[i]["analyte"]+'</td>'
			specimen_assay_html += '<td>'+specimen_assay_json[i]["assay"]+'</td>'
			specimen_assay_html += '<td>'+specimen_assay_json[i]["quality"]+'</td>'
			specimen_assay_html += '</tr>'

		genome_wide = res_data[0]['genome_wide'].replace("\'", "\"")

		genome_wide_json = json.loads(genome_wide)
		
		for j in genome_wide_json: 
			
			result_data = genome_wide_json[j]["result"] if 'result' in genome_wide_json[j] else ''

			assessment = 'Yes' if genome_wide_json[j]["assessment"] == 'Possible' else ( 'No' if genome_wide_json[j]["assessment"] == 'Not possible' else '')
			genome_wide_html += '<tr>'
			genome_wide_html += '<td>'+genome_wide_json[j]["title"]+'</td>'
			genome_wide_html += '<td>'+result_data+'</td>'
			genome_wide_html += '<td>'+assessment+'</td>'


			if genome_wide_json[j]["title"] == 'OTHER GENOMIC PHENOTYPE' and genome_wide_json[j]["assessment"] == 'Possible' and result_data == 'Yes':
				genome_wide_html += '<td>'+genome_wide_json[j]["categories"]+'</td>'
			else:
				genome_wide_html += '<td>'+genome_wide_json[j]["comment"]+'</td>'
			genome_wide_html += '</tr>'
	
	return patient_info_html, specimen_assay_html, genome_wide_html, summary_txt


# ### Build a Small variant Json (Somatic & Germline)

def build_small_variants(root_path):

	file_path = root_path + '/'
	smv_file_list = list(filter(lambda x: x.endswith('-igvnav-input.txt') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))
	regex = '^(?:(?!-(CFDNA|T)).)*igvnav-input.txt$'
	
	smt_variant_html = ''
	for i in smv_file_list:
		smv_filename = file_path + i
		smv_df_data = pd.read_csv(smv_filename, delimiter = "\t")

		if 'CALL' in smv_df_data.columns:
			
			smv_df_call_filter = smv_df_data.loc[(smv_df_data['CALL'] == "S") | (smv_df_data['CALL'] == "G")]

			column_list = ['CHROM', 'START', 'END', 'REF', 'ALT', 'GENE', 'CONSEQUENCE', 'HGVSp']
			column_dict = {'CHROM': 'chr', 'START': 'start', 'END': 'end', 'REF' : 'ref', 'ALT' : 'alt'}
			
			if 'CLONALITY' in smv_df_call_filter.columns:
				column_list.append('CLONALITY')
				column_dict['CLONALITY'] = 'clonality'

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
				
			smv_df_call_filter = smv_df_call_filter[column_list]
			smv_df_filter_data = smv_df_call_filter.rename(columns=column_dict)
			
			smv_df_filter_data.fillna('-', inplace=True)
			
			source_type = ''
			
			if(re.match(regex, i)):
				source_type ='germline'
			else:
				source_type ='somatic'

			for index, row in smv_df_filter_data.iterrows():

				variant_det = "chr"+str(row['chr'])+":"+str(row['start'])+', '+row['ref']+'>'+row['alt'].lstrip('[').rstrip(']')
				clonality = row['clonality'] if 'clonality' in row else '-'
				transcript = row['TRANSCRIPT'] if 'TRANSCRIPT' in row else '-'
				assessment = row['ASSESSMENT'] if 'ASSESSMENT' in row else '-'
				HGVSp_org_rx = row['HGVSp_org'] if 'HGVSp_org' in row else '-'
				HGVSp_org = HGVSp_org_rx.split("p.")[1] if 'p.' in HGVSp_org_rx else HGVSp_org_rx
				second_hit = str(row['second_hit']) if 'second_hit' in row else '-'

				smt_variant_html += '<tr>'
				smt_variant_html +='<td>'+row['GENE']+'</td>'
				smt_variant_html +='<td>'+source_type+'</td>'
				smt_variant_html +='<td>'+variant_det+'</td>'
				smt_variant_html +='<td>'+row['CONSEQUENCE']+'</td>'
				smt_variant_html +='<td>'+clonality+'</td>'
				smt_variant_html +='<td>'+second_hit+'</td>'
				smt_variant_html +='<td>'+assessment+'</td>'
				smt_variant_html +='<td>'+transcript+'</td>'
				smt_variant_html +='<td>'+HGVSp_org+'</td>'
				smt_variant_html += '</tr>'
				
	return smt_variant_html


# ### Build a CNV Json (Somatic & Germline)

def build_cnv(root_path):
	
	cnv_df = pd.DataFrame()
	file_path = root_path + '/cnv/'
	cnv_file_list = list(filter(lambda x: x.endswith('_curated.cns') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))
	
	cnv_html = '' 
	
	if len(cnv_file_list) > 0:
		
		regex = '^(?:(?!-(CFDNA|T)).)*_curated.cns$'
		for i in cnv_file_list:
			cnv_filename = file_path + i
			cnv_df_data = pd.read_csv(cnv_filename, delimiter = "\t")
			
			column_list = ['chromosome', 'start', 'end', 'gene', 'ASSESSMENT']
			column_dict = {'chromosome': 'chr'}
		
			if 'ASSESSMENT' in cnv_df_data.columns:
				cnv_df_data = cnv_df_data.loc[(cnv_df_data['ASSESSMENT'].notnull())]

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
					
					gene_det = row['gene']
					
					copy_number = row['copy_number'] if isinstance(row['copy_number'],str) else int(row['copy_number'])

					cnv_html += '<tr>'
					cnv_html +='<td>'+gene_det+'</td>'
					cnv_html +='<td>'+source_type+'</td>'
					cnv_html +='<td>'+variant_det+'</td>'
					cnv_html +='<td>'+row['ASSESSMENT']+'</td>'
					cnv_html +='<td>'+str(copy_number)+'</td>'
					cnv_html += '</tr>'
					

	return cnv_html


# ### Build a SVS Json

def build_svs(root_path):
	file_path = root_path + '/svs/igv/'
	
	svs_html = ''
	try:
		svs_filename = file_path + list(filter(lambda x: (re.match('[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+-sv-annotated.txt', x) or x.endswith('_annotate_combined_SV.txt')) and not x.startswith('.') and not x.endswith('.out'),os.listdir(file_path)))[0]

		sample_list = ['germline', 'somatic', 'tumor', 'cfdna']
		svs_filter = pd.read_csv(svs_filename, delimiter = "\t")
		
		column_list = ['CHROM_A', 'START_A', 'END_A', 'GENE_A', 'CHROM_B', 'START_B', 'END_B', 'GENE_B', 'SAMPLE', 'IGV_COORD', 'CLONALITY', 'SECONDHIT', 'CALL']
		column_dict = {}
		
		if 'CALL' in svs_filter.columns:
			svs_filter = svs_filter.loc[(svs_filter['CALL'] == True ) | (svs_filter['CALL'] == 'true')]
			
			if not svs_filter.empty:
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


					chr_b = int(row['CHROM_B']) if row['CHROM_B'] != '-' else row['CHROM_B']
					start_b = int(row['START_B']) if row['START_B'] != '-' else row['START_B']
					end_b = int(row['END_B']) if row['END_B'] != '-' else row['END_B']

					variant_det = "chr"+str(row['CHROM_A'])+":"+str(int(row['START_A']))+"-"+str(int(row['END_A']))+','+"chr"+str(chr_b)+":"+str(start_b)+"-"+str(end_b)
					
					gene_det = row['GENE_A']+","+ row['GENE_B']
					
					svs_html += '<tr>'
					svs_html +='<td>'+row['GENE_A']+'</td>'
					svs_html +='<td>'+row['GENE_B']+'</td>'
					svs_html +='<td>'+row['SAMPLE']+'</td>'
					svs_html +='<td>'+variant_det+'</td>'
					svs_html +='<td>'+row['consequence']+'</td>'
					svs_html +='<td>'+row['CLONALITY']+'</td>'
					svs_html +='<td>'+str(row['SECONDHIT'])+'</td>'
					svs_html += '</tr>'
				
	except Exception as e:
		print(" SVS Exception", str(e))
		svs_html = ''
	
	return svs_html


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
#         tech_html +='<tr>'
#         tech_html +='<th>INPUT DNA (ng)</th>'
#         tech_html +='<td>100</td>'
#         tech_html +='<td>100</td>'
#         tech_html +='</tr>'
		tech_html +='<tr>'
		tech_html +='<th>HYBRIDISATION CAPTURE</th>'
		if project_name == 'PROBIO':
			tech_html +='<td>Prostate cancer comprehensive panel v2</td>'
			tech_html +='<td>Prostate cancer comprehensive panel v2</td>'
		else:
			tech_html +='<td>Pan-cancer comprehensive panel (GMCK) v1</td>'
			tech_html +='<td>Pan-cancer comprehensive panel (GMCK) v1</td>'
 
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
			tech_html +='<td>'+str(math.ceil(qc_df_json["duplication"][0]))+'</td>'
			if(df_len > 1):
				tech_html +='<td>'+str(math.ceil(qc_df_json["duplication"][1]))+'</td>'
			else:
				tech_html +='<td>-</td>'
			tech_html +='</tr>'
		
	
	except Exception as e:
		print("TECHNICAL VALIDATION Exception", str(e))


	return tech_html, tech_header_html


# ### Build HTML from output files

def build_html(root_path, file_name, project_name, cfdna, capture_format, base_html_path, sample_id,capture_id, appendix_page, appendix_name):
	
	print(capture_format)
	print("--- MTBP Json Format Started ---\n")
	print("Path : ", root_path, "/", file_name)
	
	project_json = {}
		
	# Patient, Specimen & assay and genome wide information 
	logging.info('--- Patient, Specimen & assay and genome wide information started ---')
	pat_html, specimen_html, genome_wide_html, summary_txt = build_basic_html(sample_id,capture_id)

   
	# Small Variant (Somatic & Germline)
	logging.info('--- Small Variant started ---')
	small_variant_html = build_small_variants(root_path)
	if(small_variant_html == ""):
		small_variant_html = '<tr><td class="text-center" colspan="9">No relevant variants detected</td></tr>'
	logging.info('--- Small Variant completed ---')
	
	# SVS
	logging.info('--- SVS started ---')
	svs_html = build_svs(root_path)
	if(svs_html == ""):
		svs_html = '<tr><td class="text-center" colspan="7">No relevant variants detected</td></tr>'
	logging.info('--- SVS completed ---')
	
	# CNVs
	logging.info('--- CNVs started ---')
	cnv_html = build_cnv(root_path)
	if(cnv_html == ""):
		cnv_html = '<tr><td class="text-center" colspan="5">No relevant variants detected</td></tr>'
	logging.info('--- CNVs completed ---')
	
	## Get Technical Validation for QC
	logging.info('--- Technical Validation started ---')
	tech_valid_html, tech_header_html = build_tech_val_QC(root_path, project_name, capture_id)
	if(tech_valid_html == ""):
		tech_valid_html = '<tr><td class="text-center" colspan="3">No relevant variants detected</td></tr>'
	logging.info('--- Technical Validation completed ---')
	

	# Read a base html and replace the curated text based on ids
	with open(base_html_path, 'r') as f:
		contents = f.read()
		soup = BeautifulSoup(contents, "html.parser")
		
		for tag in soup.find_all(id='patient_info_table_data'):
			tag.string.replace_with(pat_html)
		
		for tag in soup.find_all(id='specimen_assay_table_data'):
			tag.string.replace_with(specimen_html)
		
		for tag in soup.find_all(id='genome_wide_table_data'):
			tag.string.replace_with(genome_wide_html)
			
		for tag in soup.find_all(id='point_table_data'):
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
	
	with open(appendix_page, 'r') as f:
		contents = f.read()
		soup1 = BeautifulSoup(contents, "html.parser")

		for tag in soup1.find_all(id='tech_val_table_header'):
			tag.string.replace_with(tech_header_html)
		
		for tag in soup1.find_all(id='tech_val_table_data'):
			tag.string.replace_with(tech_valid_html)
		
		appendix_text = soup1.prettify(formatter=None)
		  
	# Create a new appendix html based on the appendix template
	with open(appendix_name, "w", encoding = 'utf-8') as file1:
		file1.write(str(appendix_text))
	
	return 1


# ### Change the header value

def change_header_logo(header_html, project_name, output_header, specimen, sample_date):
	
	report_date = date.today().strftime('%Y-%m-%d')
	#sample_date = datetime.datetime.strptime(sample_date, '%Y%m%d').date()

	with open(header_html, 'r') as f:
		contents = f.read()
		soup = BeautifulSoup(contents, "html.parser")
		
		for tag in soup.find_all(id='specimen_data'):
			specimen_html = '<span>'+specimen+'</span>'
			tag.string.replace_with(specimen_html)
		
		for tag in soup.find_all(id='recieved_date'):
			sample_date_html = '<span>'+str(sample_date)+'</span>'
			tag.string.replace_with(sample_date_html)
			
		for tag in soup.find_all(id='report_date'):
			report_date_html = '<span>'+report_date+'</span>'
			tag.string.replace_with(report_date_html)
		
		for images in soup.find_all('img'):
			
#             if 'IPCM' in project_name or 'iPCM' in project_name:
			if 'PROBIO' in project_name:
				img_path = '/nfs/IPCM/script/pdf/static/img/probio.png'
			else:
				img_path = '/nfs/IPCM/script/pdf/static/img/iPCM.png'
				
			images['src'] = images['src'].replace("", img_path)

		new_text = soup.prettify(formatter=None)
				
	with open(output_header, "w") as f_output:
		f_output.write(str(new_text))


# ### Main Function

def main(nfs_path, project_name, sample_id, capture_id):
	
   
	## PDF Tempalte Path
	### Local Path 
	#base_html_path = '/home/karthick/project/code/curator/pdf/base.html'
	#appendix_page = '/home/karthick/project/code/curator/pdf/appendix_c4.html'
		
	### Serve Path 
	base_html_path = '/nfs/IPCM/script/pdf/base.html'
	tml_header_page = '/nfs/IPCM/script/pdf/layout/header.html' 
	footer_page = '/nfs/IPCM/script/pdf/layout/footer.html'
	appendix_page = '/nfs/IPCM/script/pdf/appendix_c4.html'
	tml_appendix_header_page = '/nfs/IPCM/script/pdf/layout/appendix_header.html' 
	 
	root_path = os.path.join(nfs_path,sample_id,capture_id)
	
	output_path = root_path+"/pdf"
	
	## Check the sample is 'PN' and 'C4' design
	KN_value = capture_id.split("-")[6]
	
	if "PN" in KN_value:
		#appendix_page = '/home/karthick/project/code/curator/pdf/appendix_pn.html'
		appendix_page = '/nfs/IPCM/script/pdf/appendix_pn.html'
	
	cfdna = capture_id.split("-")[4]
	
	capture_format = capture_id.split("-")[0]
	   
	file_name = output_path+"/base_"+project_name+"_"+sample_id+"_"+cfdna+".html"
	appendix_name = output_path+"/appendix_"+project_name+"_"+sample_id+"_"+cfdna+".html" 
	header_page = output_path+"/header_"+project_name+"_"+sample_id+"_"+cfdna+".html"
	appendix_header_page = output_path+"/append_header_"+project_name+"_"+sample_id+"_"+cfdna+".html"
	pdf_file_name = output_path+"/pdf_"+project_name+"_"+sample_id+"_"+cfdna+".pdf"
	log_name = output_path+"/log_"+project_name+"_"+sample_id+"_"+cfdna+".log"
	
	if(not os.path.exists(output_path)):
		os.mkdir(output_path)
	
	capture_arr = capture_id.split("_")
	specimen =  'cfDNA' if 'CFDNA' in capture_arr[0] else ( 'FFPE' if 'T' in capture_arr[0] else '')
		
	# Sample Information
	if(project_name == "ICPM" or capture_format == "iPCM"):
		sample_details_json = build_icpm_sample_details(cfdna)
	else:
		sample_details_json = build_sample_details(cfdna)
	
	sample_date = '-'
	
	if(sample_details_json):
		sample_date = sample_details_json["sample_date"]
		
	## Change the logo based on the project  
	change_header_logo(tml_header_page, project_name, header_page, specimen, sample_date)
	change_header_logo(tml_appendix_header_page, project_name, appendix_header_page, specimen, sample_date)
	   
	logging.basicConfig(format = '%(asctime)s  %(levelname)-10s %(name)s %(message)s', level=logging.INFO , filename=log_name, datefmt =  "%Y-%m-%d %H:%M:%S")
	logging.info('--- Generated Json format Started---')
	
	logging.info("Sample Id : {} || Capture Id : {} || Outpue File Name : {} ".format(sample_id,capture_id, file_name))
		
	try:
		html_result = build_html(root_path, file_name, project_name, cfdna, capture_format, base_html_path, sample_id,capture_id, appendix_page, appendix_name)
		if html_result:
			
			cmd = 'wkhtmltopdf --enable-local-file-access {} --header-html {} --footer-line --footer-html {} {} --header-html {} --footer-line  --footer-html {} {}'.format(file_name, header_page, footer_page, appendix_name, appendix_header_page, footer_page, pdf_file_name)
			subprocess_cmd(cmd)
			logging.info("PDF Generated")

	except Exception as e:
		print(" Main Exception", str(e))  
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
	
