from franken_api.database import db
from franken_api.database.models import ProbioBloodReferral as probio
from franken_api.database.models import PSFFBloodReferral as psff
from franken_api.database.models import TableIgvGermline as igv_germline_table
from franken_api.database.models import TableIgvSomatic as igv_somatic_table
from franken_api.database.models import TableSVS as svs_table
from franken_api.database.models import TableIgvHotspotUpdate as igv_hotspot_update_table
from franken_api.database.models import TableIgvHotspot as igv_hotspot_table
from franken_api.database.models import TableIgvWarmspot as igv_warmspot_table
from franken_api.database.models import TablePsffSummary as psff_profile
from franken_api.database.models import TableProbioSummary as probio_profile
from franken_api.database.models import TableGenomicProfileSummary as genomic_profile


from sqlalchemy import and_
import os, io
#from franken_api.settings import MOUNT_POINT
from flask import current_app
from franken_api.util_notfound import Not_found
import json
import csv
import re
import ast
from flask import jsonify
import subprocess
from collections import OrderedDict
from datetime import datetime
import pandas as pd
import sys
from flask import request
import requests
from sqlalchemy import create_engine

# check the string contains special character or not 
def check_special_char(seq_str):
	result = any(not c.isalnum() for c in seq_str)
	return result

# split the sequence into three letter and convert into one letter
def get_three_to_one_amino_code(code_seq):
	amino_code_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
	 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
	 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

	one_code_res = ''
	
	# split the code based on digit and special character
	three_code = re.split('([0-9]+|[a-zA-Z \s\n\.]+)', code_seq)
	three_code = [i for i in three_code if i]
	
	for c in three_code:
		validate_seq = check_special_char(c)
		if(not validate_seq):
			if not c.isdigit() and c.upper() in amino_code_dict:
				code = amino_code_dict[c.upper()]
				c = code
		one_code_res = one_code_res + c
	
	return one_code_res

def generate_list_to_dict(result):
	d, row = {}, []
	for rowproxy in result:
		for column, value in rowproxy.items():
			if isinstance(value, datetime):
				value = datetime.strftime(value, "%Y-%m-%d")
			else:
				value = str(value)
			d = {**d, **{column: value}}
		row.append(d)
	
	return row
	
def run_cmd(cmd):
	"Run external commands"
	return subprocess.check_output(cmd, shell=True).decode('utf-8')

def check_nfs_mount(file_path=None):
	"Check anchorage is mount to /nfs"
	if os.path.exists(file_path) and len(file_path) > 0:
		return  True, 200
	else :
		return False, 400

def get_sample_design_ids(project_path, sample_id):
	"get all sample design ids for given sample id"
	capture_dir = project_path + '/' + sample_id
	status, error =  check_nfs_mount(capture_dir)
	if not status:
		return {'sample_capture': [], 'status': False}, error

	sample_capture_list = list(filter(lambda x: (x.startswith('PB-') or x.startswith('LB-') or x.startswith('AL-') or x.startswith('OT-') or x.startswith('PSFF-') or x.startswith('RB-') or x.startswith('iPCM-')),
				os.listdir(capture_dir)))

	if len(sample_capture_list) < 1:
		return {'sample_capture': [], 'status': False}, 400

	return {'sample_capture': sample_capture_list, 'status': True}, 200


def get_sample_ids(project_path):
	"Get all probio samples"
	status, error = check_nfs_mount(project_path)
	if status:
		files = list(filter(lambda x: x.startswith('P-'), os.listdir(project_path)))
		files.sort(key=lambda x : os.path.getmtime(project_path + '/' + x), reverse=True)
		return {'sidis': files, 'status': True}, error

	return {'sidis': [], 'status': False}, error

def check_franken_plot_link(url_list):
	valid_url_list = []
	for url in url_list:
		try:
			response = requests.get(url)
			status = response.status_code
			valid_url_list.append(url)
		except requests.ConnectionError:
		   print("Failed to connect : {}".format(url))
	return valid_url_list

def get_static_frankenplot(project_path, project_name, sample_id, capture_id):
	"return the static franken plots as image"

	file_path = project_path + '/' + sample_id + '/' + capture_id + '/qc/'
	temp_url_list = []
	ip_addr = 'localhost' if '5000' in request.host else request.host
	port_no = ':5000' if 'localhost' in ip_addr else ''
	status = True if os.path.exists(file_path) and len(os.listdir(file_path)) > 0 else False
	if status:
		for each_file in filter(lambda x: x.endswith('liqbio-cna.png') and not x.startswith('.'), os.listdir(file_path)):
			link = 'http://'+ip_addr+port_no+'/'+ 'api/franken/staticimage?project_name=' + project_name + '&sdid=' + sample_id + '&capture_id=' + capture_id + '&imagename=' + each_file
			temp_url_list.append(link)

		if len(temp_url_list) > 0:
			image_url_list = temp_url_list
			return {'image_url': image_url_list, 'status': True}, 200

	return {'image_url':[], 'status': False}, 400


def get_static_image(project_path, sample_id, capture_id, image_name):
	"""retrun the franken static png image"""

	file_path = project_path + '/' + sample_id + '/' + capture_id + '/qc/' + image_name
	
	if os.path.exists(file_path):
		return file_path, 200

	return file_path, 400

def get_xml_image(project_path, sample_id, capture_id, xml):
	"""retrun the franken static png image"""

	# file_path = project_path + '/' + sample_id + '/' + capture_id + '/IGVnav/' + xml
	file_path = "/home/karman/igv.js/examples/sessions/IGVnav/igv_session_sv.xml"

	if os.path.exists(file_path):
		return file_path, 200

	return file_path, 400

def get_interactive_plot(project_path, sample_id, capture_id, plotname):
	'retrun the json files to plot interative frankenplots'
	file_path = project_path + '/' + sample_id + '/' + capture_id + '/qc/'
	status = True if os.path.exists(file_path) and len(os.listdir(file_path)) > 0 else False
	if status:
		for each_file in filter(lambda x: x.endswith('_' + plotname + '.json'),
								os.listdir(file_path)):
			with open(file_path + each_file, 'r') as f:
				return {'plot_data': json.load(f), 'status': True}, 200

	return {'plot_data': {}, 'status': False}, 400

def generate_headers_table_sv(headers):
	head_obj = {}
	for each_head in headers:
		head_obj[each_head] = {
			'title':each_head,
			'type' : 'text',
			'editable': False
		}

	return head_obj

def generate_headers_ngx_table(headers):
	columns= []
	for each_head in headers:
		  columns.append({ 'key': each_head,'title':each_head.upper().replace("_"," ")})
	return columns


def get_table_qc_header(project_path, sdid, capture_id, header='true'):
	"read qc file from qc_overview.txt and return as json"
	data = []
	file_path = project_path + '/' + sdid + '/' + capture_id + '/' + 'qc/'
	qc_filename = file_path + list(filter(lambda x: x.endswith('.qc_overview.txt')
													and not x.startswith('.')
													and not x.endswith('.out'), os.listdir(file_path)))[0]
	if os.path.exists(qc_filename):
		hide_header = ["indexs",  "FOLD_80_BASE_PENALTY", "dedupped_on_bait_rate"]
		has_rows = False
		data = []
		with open(qc_filename, 'r') as f:
			reader_ponter = csv.DictReader(f, delimiter ='\t')
			for i, each_row in enumerate(reader_ponter):
				has_rows = True
				each_row = dict(each_row)
				each_row['indexs'] = i
				data.append(each_row)
			
		if(not has_rows):
			return {'header': [], 'data': [], 'status': True, 'error': 'No Data Found For QC'}, 200

		header = list(generate_headers_ngx_table(data[0].keys()))

		new_keys = {
			'CHIP': {'key': 'CHIP', 'title': 'CHIP'},    
			'PURITY': {'key': 'PURITY', 'title': 'PURITY'},
			'PLOIDY': {'key': 'PLOIDY', 'title': 'PLOIDY'},
			'Overall_QC': {'key': 'Overall_QC', 'title': 'Overall_QC'},
			'Comment': {'key': 'Comment', 'title': 'Comment'}
		}

		for idx,value in enumerate(new_keys):
			n_key = [item for item in header if item.get('key')==value]
			if(not n_key):
				header.append(new_keys[value])  

		new_header = []
		for i,value in enumerate(header):
			key_name = value['key']
			if(key_name not in hide_header):
				new_header.append(value)

		return {'header': new_header, 'data': data, 'filename': qc_filename, 'status': True}, 200

	else:
		return {'header': [], 'data': [], 'filename': '', 'status': False}, 400


def get_table_svs_header(project_path, sdid, capture_id, header='true'):
	"read structural variant file from sdid_annotate_combined_SV.txt and return as json"
	file_path = project_path + '/' + sdid + '/' + capture_id + '/svs/igv/'

	file_path = file_path + list(filter(lambda x: (re.match('[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+-sv-annotated.txt$', x) or
									   x.endswith('_annotate_combined_SV.txt'))
									  and not x.startswith('.')
									  and not x.endswith('.out'),
							os.listdir(file_path)))[0]
	data = []
	if os.path.exists(file_path):
		hide_header = ["GENE_A-GENE_B-sorted", "CAPTURE_ID", "PROJECT_ID", "SDID", "indexs"]
		df = pd.read_csv(file_path,delimiter="\t")

		if(df.empty):
			return {'header': [], 'data': [], 'status': True, 'error': 'No Data Found For Structural Variants'}, 200
		else:
			# Dataframe soted based on the below columns
			df_sorted = df.sort_values(["GENE_A-GENE_B-sorted","CHROM_A","START_A","CHROM_B","START_B","TOOL","SUPPORT_READS"], ascending = (True,False,False,False,False,False,False))
			df_filter = df_sorted.loc[(df['IN_DESIGN_A'] == 'YES') | (df['IN_DESIGN_B'] == 'YES') | (df['TOOL'] == 'svcaller')]
			
			# Add Index column in the dataframe
			df_filter['indexs'] = pd.RangeIndex(len(df_filter.index))
			

			# if "IGV_COORD" in df_filter.columns:
			# 	df_filter = df_filter[['CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B', 'END_B', 'IGV_COORD', 'SVTYPE', 'SV_LENGTH', 'SUPPORT_READS', 'TOOL', 'SDID', 'SAMPLE', 'GENE_A', 'GENE_B', 'IN_DESIGN_A', 'IN_DESIGN_B', 'GENE_A-GENE_B-sorted', 'indexs']]
			# else:
			# 	df_filter = df_filter[['CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B', 'END_B', 'SVTYPE', 'SV_LENGTH', 'SUPPORT_READS', 'TOOL', 'SDID', 'SAMPLE', 'GENE_A', 'GENE_B', 'IN_DESIGN_A', 'IN_DESIGN_B', 'GENE_A-GENE_B-sorted', 'indexs']]

			column_list = list(df_filter.columns)

			result = df_filter.to_json(orient="records")
			data = json.loads(result)

			cal_key = 'CALL'
			typ_key = 'TYPE'
			sec_key = 'SECONDHIT'
			com_key = 'COMMENT'
			asst_key = 'ASSESSMENT'
			clon_key = 'CLONALITY'
			consq_key = 'CONSEQUENCE'
			funtp_key = 'FUNCTIONAL_TYPE'
			vatStr_key = 'VARIANT_STRING'

			if cal_key in column_list:
				cal_indx = column_list.index(cal_key)
				del column_list[cal_indx]
				column_list.insert(0,cal_key)

			if typ_key in column_list:
				typ_indx = column_list.index(typ_key)
				del column_list[typ_indx]
				column_list.insert(0,typ_key)
				
			if sec_key in column_list:
				sec_indx = column_list.index(sec_key)
				del column_list[sec_indx]
				column_list.insert(0,sec_key)    

			if com_key in column_list:
				com_indx = column_list.index(com_key)
				del column_list[com_indx]
				column_list.insert(0,com_key)

			if asst_key in column_list:
				asst_indx = column_list.index(asst_key)
				del column_list[asst_indx]
				column_list.insert(0,asst_key)
			
			if clon_key in column_list:
				clon_indx = column_list.index(clon_key)
				del column_list[clon_indx]
				column_list.insert(0,clon_key)

			if funtp_key in column_list:
				funtp_indx = column_list.index(funtp_key)
				del column_list[funtp_indx]
				column_list.insert(0,funtp_key)

			if consq_key in column_list:
				consq_indx = column_list.index(consq_key)
				del column_list[consq_indx]
				column_list.insert(0,consq_key)

			if vatStr_key in column_list:
				vatStr_indx = column_list.index(vatStr_key)
				del column_list[vatStr_indx]
				column_list.insert(0,vatStr_key)

			header = list(generate_headers_ngx_table(column_list))
			
			#Add additional columns to SV  [CALL(True | False):  TYPE:(Somatic| germline) and comment columns]
			new_keys = {
				cal_key: {'key': cal_key, 'title': 'CALL'},
				typ_key: {'key': typ_key, 'title': 'TYPE'},
				sec_key :  {'key': sec_key, 'title': 'SECONDHIT'},
				com_key :  {'key': com_key, 'title': 'COMMENT'},
				asst_key :  {'key': asst_key, 'title': 'ASSESSMENT'},
				clon_key :  {'key': clon_key, 'title': 'CLONALITY'},
				consq_key :  {'key': consq_key, 'title': 'CONSEQUENCE'},
				funtp_key :  {'key': funtp_key, 'title': 'FUNCTIONAL TYPE'},
				vatStr_key :  {'key': vatStr_key, 'title': 'VARIANT STRING'}
			}

			for idx,value in enumerate(new_keys):
				n_key = [item for item in header if item.get('key')==value]
				if(not n_key):
					header.insert(0, new_keys[value])

			new_header = []
			for i,value in enumerate(header):
				key_name = value['key']
				if(key_name not in hide_header):
					new_header.append(value)

			return {'header': new_header, 'data': data, 'filename': file_path, 'status': True}, 200		
		
		#====== Start : Old code for structural variant ===========#
		'''
		with open(file_path, 'r') as f:
			reader_ponter = csv.DictReader(f, delimiter ='\t')
			for i, each_row in enumerate(reader_ponter):
				each_row = dict(each_row)
				each_row['indexs'] = i
				data.append(each_row)
			header = list(generate_headers_ngx_table(data[0].keys()))

			#Add additional columns to SV  [CALL(True | False):  TYPE:(Somatic| germline) and comment columns]
			new_keys = {
				'CALL': {'key': 'CALL', 'title': 'CALL'},
				'TYPE': {'key': 'TYPE', 'title': 'TYPE'},
				'SECONDHIT': {'key': 'SECONDHIT', 'title': 'SECONDHIT'},
				'COMMENT': {'key': 'COMMENT', 'title': 'COMMENT'}
			}
			
			for idx,value in enumerate(new_keys):
				n_key = [item for item in header if item.get('key')==value]
				if(not n_key):
					header.insert(0, new_keys[value])

			return {'header': header, 'data': data, 'filename': file_path, 'status': True}, 200
		 '''   
		#====== End : Old code for structural variant ===========#

	else:
		return {'header': [], 'data': [], 'filename': '', 'status': False}, 400


def get_HGVSp_status_update(gene, HGVSp):
	
	hotspot_status = ''

	hotspot_list = get_hotspot_update_info(gene)

	if(hotspot_list):
		hotspot_list =  [x.strip(' ') for x in hotspot_list.split(',')]
		if(HGVSp in hotspot_list):
			hotspot_status = '0'
			
	return hotspot_status

def get_HGVSp_status(gene, HGVSp):
	'''
		hotspot status - '', 0, 1, 2
			'' - empty
			0 - hotspot
			1 - warmspot
			2- both hotspot and warmspot
	'''

	HGVSp = 'p.'+HGVSp 

	hotspot_list = get_hotspot_info(gene)

	warmspot_list = get_warmspot_info(gene)

	hotspot_status = ''
	warmspot_status = ''

	if(hotspot_list):
		hotspot_list =  [x.strip(' ') for x in hotspot_list.split(',')]
		hotspot_arr = list(filter(lambda x: re.findall('^({})(.+)$'.format(x), HGVSp), hotspot_list))
		
		if(HGVSp in hotspot_list or (''.join(hotspot_arr) in hotspot_list and hotspot_arr != []) ):
			hotspot_status = 0

	if(warmspot_list):
		warmspot_list =  [x.strip(' ') for x in warmspot_list.split(',')]
		warmspot_arr = list(filter(lambda x: re.findall('^({})(.+)$'.format(x), HGVSp), warmspot_list))

		if(HGVSp in warmspot_list or (''.join(warmspot_arr) in warmspot_list and warmspot_arr != []) ):
			warmspot_status = 1 

	status = ''

	if(hotspot_status == 0 and warmspot_status == ''):
		status = '0'
	elif(hotspot_status == '' and warmspot_status == 1):
		status = '1'
	elif(hotspot_status == 0 and warmspot_status == 1):
		status = '0'

	return status

def get_hotspot_update_info(gene):
	res = igv_hotspot_update_table.query.filter(igv_hotspot_update_table.gene == '{}'.format(gene)).all()
	protmut = ''
	if(res):
		for r in res:
			protmut += r.variant_arr + ','
	return protmut

def get_hotspot_info(gene):
	res = igv_hotspot_table.query.filter(igv_hotspot_table.gene == '{}'.format(gene)).all()
	protmut = ''
	if(res):
		for r in res:
			protmut += r.protmut + ','
	return protmut.rstrip(',')

def get_warmspot_info(gene):
	res = igv_warmspot_table.query.filter(igv_warmspot_table.gene == '{}'.format(gene)).all()
	protmut = ''
	if(res):
		for r in res:
			# protmut.append(r.protmut)
			protmut += r.protmut + ','
	return protmut.rstrip(',')

def get_table_igv(variant_type, project_path, sdid, capture_id, header='true'):
	"read  variant file for given sdid and return as json"
	file_path = project_path + '/' + sdid + '/' + capture_id
	data = []
	missing_header = []
	header = []

	if variant_type == 'germline':
		missing_header = ['purecn_status', 'purecn_probability', 'purecn_tot_copies']
		regex = '^(?:(?!-(CFDNA|T)).)*igvnav-input.txt$'
		regex2 = '(.*)-(CFDNA|T)-(\w.*)(germline-igvnav-input).*txt$'
	elif variant_type == 'somatic':
		missing_header = ['GENE', 'IMPACT', 'CONSEQUENCE', 'HGVSp', 'HGVSp_org', 'T_DP', 'T_ALT', 'T_VAF', 'N_DP', 'N_ALT', 'N_VAF', 'CLIN_SIG', 'gnomAD', 'BRCAEx', 'OncoKB', 'purecn_probability', 'purecn_status', 'purecn_tot_copies']
		regex = '.*-(CFDNA|T)-.*igvnav-input.txt$'
		regex2 = '(.*)-(CFDNA|T)-(\w.*)(somatic-igvnav-input).*txt$'
	else:
		missing_header = ['GENE', 'IMPACT', 'CONSEQUENCE', 'HGVSp', 'HGVSp_org', 'N_DP', 'N_ALT', 'N_VAF', 'CLIN_SIG', 'gnomAD', 'BRCAEx', 'OncoKB']
		return {'header': {}, 'data': [], 'status': False, 'error': 'unknown variant type: ' + variant_type}, 400

	try:
		hide_header = ["HGVSp_org", "PureCN_probability", "PureCN_status", "PureCN_tot_copies", "purecn_probability", "purecn_status", "purecn_tot_copies", "indexs", "CAPTURE_ID", "PROJECT_ID", "SDID"]
		igv_nav_file = list(filter(lambda x: (re.match(regex2, x) if ('-somatic-' in x or '-germline-' in x) else re.match(regex, x) ) and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))[0]
		igv_nav_file = file_path + '/' + igv_nav_file

		has_rows = False
		data = []
		with open(igv_nav_file, 'r') as f:
			reader_pointer = csv.DictReader(f, delimiter='\t')
			for i, each_row in enumerate(reader_pointer):
				has_rows = True
				each_row = dict(each_row)
				each_row['indexs'] = i
				each_row['HGVSp_org'] = each_row['HGVSp']
				gene = each_row['GENE']
				if each_row['HGVSp'] and ':p.' in each_row['HGVSp']:
					one_amino_code = get_three_to_one_amino_code(each_row['HGVSp'].split("p.")[1])
					each_row['HGVSp'] = one_amino_code
				
				if variant_type == 'somatic' and each_row['HGVSp']:
					HGVSp_status = get_HGVSp_status(gene,each_row['HGVSp'])
					if(HGVSp_status == ''):
						HGVSp_status = get_HGVSp_status_update(gene,each_row['HGVSp'])
					each_row['HOTSPOT'] = HGVSp_status

				consequence = each_row['CONSEQUENCE'].replace('&', ' & ')
				each_row['CONSEQUENCE'] = consequence              

				if None in each_row:
					if isinstance(each_row[None], list):
						for i, each_none in enumerate(each_row[None]):
							each_row[missing_header[i]] = each_none
						#each_row['Notes'] = " ".join( each_row[None])
						del each_row[None]

				data.append(dict(each_row))

		if(not has_rows):
			return {'header': [], 'data': [], 'status': True, 'error': 'No Data Found For {} Variants'.format(variant_type.capitalize())}, 200

		header = list(data[0])
		if 'HOTSPOT' in header and variant_type == 'somatic':
			del header[header.index('HOTSPOT')]

		if 'HOTSPOT' not in header and variant_type == 'somatic':
			conseq_index = header.index('CONSEQUENCE') + 1
			header.insert(conseq_index, 'HOTSPOT')

		asec_key = 'SECONDHIT'
		ass_key = 'ASSESSMENT'
		acl_key = 'CLONALITY'

		if asec_key in header:
			asec_indx = header.index(asec_key)
			del header[asec_indx]
			header.insert(0,asec_key)

		if ass_key in header:
			ass_indx = header.index(ass_key)
			del header[ass_indx]
			header.insert(0,ass_key)

		if acl_key in header:
			acl_indx = header.index(acl_key)
			del header[acl_indx]
			header.insert(0,acl_key)
		
		header = generate_headers_ngx_table(header)

		new_header = []
		for i,value in enumerate(header):
			key_name = value['key']
			if(key_name not in hide_header):
				new_header.append(value)

		return {'header': new_header, 'data': data, 'filename' : igv_nav_file, 'status': True}, 200

	except Exception as e:
		return {'header': [], 'data': [], 'status': False, 'error': str(e)}, 400


def save_igvnav_input_file(filename, data):
	"save the igvnav input file with updated calls and tags column"
	try:
		with open(filename, 'w') as fpw:
			#headers = data[0].keys()
			headers = list(OrderedDict.fromkeys([k for i in data for k in i.keys()]))
			writer = csv.DictWriter(fpw, fieldnames=headers, delimiter='\t' )
			writer.writeheader()
			for each_row in data:
				writer.writerow(each_row)

		return str("written to file"), 200

	except Exception as e:
		return str(e), 400

def get_probio_blood_referrals():
	"Fetch the all the records from probio referral database"
	header = ['crid','pnr','rid','datum','tid','remisstyp', 'studieid', 'sign','countyletter','new','progression','follow_up','cf_dna1','cf_dna2','cf_dna3','blood','kommentar','filnamn', 'hormonkänslig', 'kastrationsresistent', 'cdk', 'endofstudy']
	try:
		return {'status': True, 'data': probio.query.filter().all(), 'header': generate_headers_ngx_table(header), 'error': '' }, 200
	except Exception as e:
		return {'status': True, 'data': [], 'header': generate_headers_ngx_table(header), 'error': str(e) }, 400

def get_psff_blood_referrals():
	"Fetch the all the records from probio referral database"
	header = ['crid','rid','datum','tid','sign','blood1','blood2','blood3','blood4','comment', 'filnamn', 'cdk']
	try:
		return {'status': True, 'data': psff.query.filter().all(), 'header': generate_headers_ngx_table(header), 'error': '' }, 200
	except Exception as e:
		return {'status': True, 'data': [], 'header': generate_headers_ngx_table(header), 'error': str(e) }, 400


def update_referrals(db_name):
	"Update the referrals data from ftp into postgres db using referral-manager tool"
	referral_conf = {
		'fetch': {
			'common': 'refman --sentry-login /nfs/PROBIO/referraldb/.sentrylogin fetch --referrals-login /nfs/PROBIO/referraldb/referral-manager_conf_files/login.json ',
			'probio': ' --local-data-dir /nfs/PROBIO/referraldb/remote_files --remote-data-dir /ProBio2/Scannade_remisser ',
			#'psff': ' --local-data-dir /nfs/CLINSEQ/PSFF/referraldb/remote_files --remote-data-dir /PSFF/Scannade_remisser ',
			'psff': ' --local-data-dir /nfs/PSFF/referraldb/remote_files --remote-data-dir /PSFF/Scannade_remisser ',
			#'psff_log': '/nfs/CLINSEQ/PSFF/referraldb/referral_db_fetch.log',
			'psff_log': '/nfs/PSFF/referraldb/referral_db_fetch.log',
			'probio_log': '/nfs/PROBIO/referraldb/referral_db_fetch.log'
		},
		'db_import':{
			'common' : 'refman --sentry-login /nfs/PROBIO/referraldb/.sentrylogin dbimport --dbcred /nfs/PROBIO/referraldb/referral-manager_conf_files/config.json ',
			'probio': ' --local-data-dir /nfs/PROBIO/referraldb/remote_files/csv --referral-type ProbioBloodReferral ',
			#'psff': '--local-data-dir /nfs/CLINSEQ/PSFF/referraldb/remote_files/csv --referral-type PsffBloodReferral ',
			'psff': '--local-data-dir /nfs/PSFF/referraldb/remote_files/csv --referral-type PsffBloodReferral ',
			#'psff_log': '/nfs/CLINSEQ/PSFF/referraldb/referral_db_dbimport.log',
			'psff_log': '/nfs/PSFF/referraldb/referral_db_dbimport.log',
			'probio_log': '/nfs/PROBIO/referraldb/referral_db_dbimport.log'

		}
	}
	try:
		logfile_fetch = open(referral_conf['fetch'][db_name+'_log'], 'w')
		logfile_dbimport = open(referral_conf['db_import'][db_name+'_log'], 'w')
		cmd_fetch = referral_conf['fetch']['common'] + referral_conf['fetch'][db_name]
		cmd_dbimport = referral_conf['db_import']['common'] + referral_conf['db_import'][db_name]
		proc = subprocess.check_call(cmd_fetch, stdout=logfile_fetch, stderr=logfile_fetch, shell=True)
		proc = subprocess.check_call(cmd_dbimport, stdout=logfile_dbimport, stderr=logfile_dbimport, shell=True)
		logfile_fetch.close()
		logfile_dbimport.close()
		return {'status': True, 'error': ''}, 200
	except subprocess.CalledProcessError as err:
		logfile_fetch.close()
		logfile_dbimport.close()
		return {'status': False, 'error': err}, 400


def pdfs_files(variant_type, project_path, sdid, capture_id):
	if variant_type not in ['qc', 'multiqc', 'purecn']:
		return '', 400

	file_path = project_path + '/' + sdid + '/' + capture_id + '/' + variant_type 

	pdf_file = list(filter(lambda x: (re.match('[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+.pdf$', x) or re.match('[-\w]+-(multiqc)+.html$', x) or x.endswith('.qc_overview.pdf')) and not x.startswith('.') and not x.endswith('.out'),
							   os.listdir(file_path)))[0]

	if os.path.exists(file_path):
		file_path = file_path + '/' + pdf_file
		return file_path, 200

	return file_path, 400

def check_curation_germline_record(table, record):
	return table.query.filter(table.PROJECT_ID==record['PROJECT_ID'],
								   table.SDID == record['SDID'],
								   table.CAPTURE_ID == record['CAPTURE_ID'],
								   table.CHROM == record['CHROM'],
								   table.START == record['START'],
								   table.END == record['END'],
								   table.REF == record['REF'],
								   table.ALT == record['ALT']
								   ).first()

def check_curation_somatic_record(table, record):
	return table.query.filter(table.PROJECT_ID==record['PROJECT_ID'],
								   table.SDID == record['SDID'],
								   table.CAPTURE_ID == record['CAPTURE_ID'],
								   table.CHROM == record['CHROM'],
								   table.START == record['START'],
								   table.END == record['END'],
								   table.REF == record['REF'],
								   table.ALT == record['ALT']
								   ).first()


def check_curation_svs_record(table, record):
	return table.query.filter(table.PROJECT_ID==record['PROJECT_ID'],
								   table.SDID == record['SDID'],
								   table.CAPTURE_ID == record['CAPTURE_ID'],
								   table.CHROM_A == record['CHROM_A'],
								   table.START_A == record['START_A'],
								   table.END_A == record['END_A'],
								   table.CHROM_B == record['CHROM_B'],
								   table.START_B == record['START_B']
								   ).first()

def check_curation_psff_profile_record(table, record):
	return table.query.filter(table.sample_id == record['sample_id'],
								   table.capture_id == record['capture_id']
								   ).first()

def check_curation_probio_profile_record(table, record):
	return table.query.filter(table.sample_id == record['sample_id'],
								   table.capture_id == record['capture_id']
								   ).first()

def check_curation_genomic_profile_record(table, record):
	return table.query.filter(table.project_name == record['project_name'],
								   table.sample_id == record['sample_id'],
								   table.capture_id == record['capture_id']
								   ).first()

def curation_update_profile(record, table_name):
	try:
		tables_dict = {
			'psff_summary':psff_profile,
			'probio_summary':probio_profile,
			'genomic_profile_summary': genomic_profile
		}
		func_dict = {
			'psff_summary':check_curation_psff_profile_record,
			'probio_summary':check_curation_probio_profile_record,
			'genomic_profile_summary':check_curation_genomic_profile_record
		}
		current_record = func_dict[table_name](tables_dict[table_name], record )
		
		if not bool(current_record):
			obj_germline = tables_dict[table_name](record)
			db.session.add(obj_germline)
			db.session.commit()
			return {'status': True, 'error': ''}, 200
		else:
			for each_col in record:
				setattr(current_record, each_col, record[each_col])
			db.session.commit()
			return {'status': True, 'error': ''}, 200
			
	except Exception as e :
		return {'status': False, 'error': str(e)}, 400

def post_curation(record, table_name):
	try:
		tables_dict = {
			'germline': igv_germline_table,
			'somatic': igv_somatic_table,
			'svs': svs_table
		}
		func_dict = {
			'germline': check_curation_germline_record,
			'somatic': check_curation_somatic_record,
			'svs': check_curation_svs_record
		}

		current_record = func_dict[table_name](tables_dict[table_name], record )

		if not bool(current_record):
			obj_germline = tables_dict[table_name](record)
			db.session.add(obj_germline)
			db.session.commit()
			return {'status': True, 'error': ''}, 200
		else:
			for each_col in record:
				setattr(current_record, each_col, record[each_col])
			db.session.commit()
			return {'status': True, 'error': ''}, 200
	except Exception as e :
		return {'status': False, 'error': str(e)}, 400


def get_curation_cancer_hotspot():
	try:
		header = ['id', 'gene', 'residue', 'res_type', 'variants', 'variant_arr']
		try:
			return {'status': True, 'data': igv_hotspot_update_table.query.filter().all(),
					'header': generate_headers_ngx_table(header),
					'error': ''}, 200
		except Exception as e:
			return {'status': True, 'data': [], 'header':  generate_headers_ngx_table(header), 'error': str(e)}, 400

	except Exception as e:
		return "Error :" + str(e), 400

def get_curation_hotspot():
	try:
		header = ['gene', 'aapos', 'nmut', 'protmut', 'prot2mut', 'dnamut', 'canmut', 'conseqmut', 'transcript', 'dn_ds', 'community_notes']
		try:
			return {'status': True, 'data': igv_hotspot_table.query.filter().all(),
					'header': generate_headers_ngx_table(header),
					'error': ''}, 200
		except Exception as e:
			return {'status': True, 'data': [], 'header':  generate_headers_ngx_table(header), 'error': str(e)}, 400

	except Exception as e:
		return "Error :" + str(e), 400

def get_curation_warmspot():
	try:
		header = ['gene', 'aapos', 'nmut', 'protmut', 'prot2mut', 'dnamut', 'canmut', 'conseqmut', 'transcript', 'dn_ds', 'community_notes']
		try:
			return {'status': True, 'data': igv_warmspot_table.query.filter().all(),
					'header': generate_headers_ngx_table(header),
					'error': ''}, 200
		except Exception as e:
			return {'status': True, 'data': [], 'header':  generate_headers_ngx_table(header), 'error': str(e)}, 400

	except Exception as e:
		return "Error :" + str(e), 400

def get_curation_igv_germline():
	try:
		header = ['PROJECT_ID', 'SDID', 'CAPTURE_ID', 'CHROM', 'START', 'END',
				  'REF', 'ALT', 'CALL', 'TAG', 'NOTES', 'GENE', 'IMPACT', 'CONSEQUENCE',
				  'HGVSp', 'N_DP', 'N_ALT', 'N_VAF', 'CLIN_SIG', 'gnomAD', 'BRCAEx', 'OncoKB', 'purecn_probability', 'purecn_status', 'purecn_tot_copies']
		try:
			return {'status': True, 'data': igv_germline_table.query.filter().all(),
					'header': generate_headers_ngx_table(header),
					'error': ''}, 200
		except Exception as e:
			return {'status': True, 'data': [], 'header':  generate_headers_ngx_table(header), 'error': str(e)}, 400

	except Exception as e:
		return "Error :" + str(e), 400

def get_curation_igv_somatic():
	try:
		header = ['PROJECT_ID', 'SDID', 'CAPTURE_ID', "CHROM", 'START', 'END',
					'REF', 'ALT', 'CALL', 'TAG', 'NOTES', 'ASSESSMENT', 'CLONALITY',  'GENE', 'IMPACT',
				  'CONSEQUENCE', 'HGVSp', 'T_DP', 'T_ALT', 'T_VAF', 'N_DP', 'N_ALT', 'N_VAF',
				  'CLIN_SIG', 'gnomAD', 'BRCAEx', 'OncoKB', 'purecn_probability', 'purecn_status', 'purecn_tot_copies']
		try:
			return {'status': True, 'data': igv_somatic_table.query.filter().all(),
					'header': generate_headers_ngx_table(header),
					'error': ''}, 200
		except Exception as e:
			return {'status': True, 'data': [], 'header': header, 'error': str(e)}, 400

	except Exception as e:
		return "Error :" + str(e), 400

def get_curation_svs():
	try:
		header = ['PROJECT_ID', 'CAPTURE_ID', 'SDID', 'CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B',
				  'END_B', 'SVTYPE', 'SV_LENGTH', 'SUPPORT_READS', 'TOOL', 'SAMPLE', 'GENE_A', 'IN_DESIGN_A', 'GENE_B',
				  'IN_DESIGN_B', 'GENE_A-GENE_B-sorted', 'CALL', 'TYPE', 'SECONDHIT', 'COMMENT', 'ASSESSMENT', 'CLONALITY', 'CONSEQUENCE', 'FUNCTIONAL_TYPE', 'VARIANT_STRING']
		try:
			return {'status': True, 'data': svs_table.query.filter().all(),
					'header': generate_headers_ngx_table(header),
					'error': ''}, 200
		except Exception as e:
			return {'status': True, 'data': [], 'header': header, 'error': str(e)}, 400

	except Exception as e:
		return "Error :" + str(e), 400

def list_curation_psff_profile():
	try:
		header = ['sample_id', 'capture_id', 'basic_qc', 'franken_plot', 'ploidy', 'ctdna_fraction', 'ctdna_param', 'ctdna_method', 'study_code', 'study_site', 'somatic_mutations', 'germline_alterations', 'structural_variants', 'cnvs', 'comment_info']
		try:
			return {'status': True, 'data': psff_profile.query.filter().all(),
					'header': generate_headers_ngx_table(header),
					'error': ''}, 200
		except Exception as e:
			return {'status': True, 'data': [], 'header': header, 'error': str(e)}, 400

	except Exception as e:
		return "Error :" + str(e), 400


def list_curation_probio_profile():
	try:
		header = ['sample_id', 'capture_id', 'basic_qc', 'franken_plot', 'ploidy', 'ctdna_fraction', 'ctdna_param', 'ctdna_method', 'ctdna_category',  'study_code', 'study_site', 'somatic_mutations', 'germline_alterations', 'structural_variants', 'cnvs', 'comment_info']
		try:
			return {'status': True, 'data': probio_profile.query.filter().all(),
					'header': generate_headers_ngx_table(header),
					'error': ''}, 200
		except Exception as e:
			return {'status': True, 'data': [], 'header': header, 'error': str(e)}, 400

	except Exception as e:
		return "Error :" + str(e), 400

def list_curation_genomic_profile():
	try:
		header = ['project_name', 'sample_id', 'capture_id', 'study_code', 'study_site', 'dob', 'disease', 'specimen_assay', 'ctdna_param', 'ctdna_method', 'genome_wide', 'somatic_mutations', 'germline_alterations', 'structural_variants', 'cnvs', 'summary_txt']
		try:
			return {'status': True, 'data': genomic_profile.query.filter().all(),
					'header': generate_headers_ngx_table(header),
					'error': ''}, 200
		except Exception as e:
			return {'status': True, 'data': [], 'header': header, 'error': str(e)}, 400

	except Exception as e:
		return "Error :" + str(e), 400

def get_curation_psff_profile(record):
	try:
		return {'status': True, 'data': psff_profile.query.filter(psff_profile.sample_id == record['sample_id']).all(), 'error': ''}, 200
	except Exception as e:
		return {'status': True, 'data': [], 'error': str(e)}, 400

def get_curation_probio_profile(record):
	try:
		return {'status': True, 'data': probio_profile.query.filter(probio_profile.sample_id == record['sample_id']).all(), 'error': ''}, 200
	except Exception as e:
		return {'status': True, 'data': [], 'error': str(e)}, 400

def get_curation_genomic_profile(record):
	try:
		return {'status': True, 'data': genomic_profile.query.filter(genomic_profile.project_name == record['project_name'] and genomic_profile.sample_id == record['sdid'] and genomic_profile.capture_id == record['capture_id']).all(), 'error': ''}, 200
	except Exception as e:
		return {'status': True, 'data': [], 'error': str(e)}, 400		

def get_table_cnv_header(project_path, sdid, capture_id, variant_type, header='true'):
	"read qc file from qc_overview.txt and return as json"
	data = []
	file_path = project_path + '/' + sdid + '/' + capture_id + '/' + 'cnv/'
	if variant_type == 'somatic':
		regex = '[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+.cns$'
		set_save_file = '_somatic_curated.cns'
	elif variant_type == 'germline':
		regex = '^(?:(?!(-CFDNA|germline_curated|-T)).)*.cns$'
		set_save_file = '_germline_curated.cns'
	else:
		return {'header': [], 'data': [], 'filename': '', 'error': 'Invalid end point', 'status': False}, 400
	
	file_list = list(filter(lambda x: (re.match(regex, x) ) and not x.startswith('.') and not x.endswith('.out'),os.listdir(file_path)))

	if file_list != []:

		cnv_filename = file_path + file_list[0]

		save_to_cnv_file  = cnv_filename.split('.cns')[0] + set_save_file

		curated_cnv_file =  list(filter(lambda x: ( x.endswith(set_save_file))
										and not x.startswith('.')
										and not x.endswith('.out'),
								os.listdir(file_path)))

		curated_file_status = True if curated_cnv_file else False

		if curated_file_status:
			cnv_filename = save_to_cnv_file
		
		if os.path.exists(cnv_filename):
			has_rows = False
			hide_header = ["ABSOLUTE_COPY_NUMBER", "COMMENT", "depth", "weight", "indexs"]
			data = []
			with open(cnv_filename, 'r') as f:
				reader_ponter = csv.DictReader(f, delimiter ='\t')
				for i, each_row in enumerate(reader_ponter):
					has_rows = True
					each_row = dict(each_row)
					each_row['indexs'] = i
					gene_list = each_row['gene'].split(",")
					glist = []
					[glist.append(x) for x in gene_list if x not in glist]
					each_row['gene'] = ', '.join(glist)
					if 'SIZE' in each_row.keys() and each_row['SIZE'].isdecimal():
						each_row['SIZE'] = '{0:.4f} Mb'.format(int(each_row['SIZE'])/1000000)
					data.append(each_row)

			if(not has_rows):
				return {'header': [], 'data': [], 'status': True, 'error': 'No Data Found For CNV {} Variants'.format(variant_type.capitalize())}, 200

			header = list(data[0])
			#compute size for cnv using start and end
			if 'SIZE' not in header:
				end_index = header.index('end') + 1
				header.insert(end_index, 'SIZE')
				for data_dict in data:
					size = int(data_dict['end']) - int(data_dict['start']) + 1
					data_dict['SIZE'] = '{0:.2f} Mb'.format(size/1000000)

			acn_key = 'ABSOLUTE_COPY_NUMBER'
			ass_key = 'ASSESSMENT'
			com_key = 'COMMENT'
			pur_key = 'PURITY'
			plo_key = 'PLOIDY'
			copy_nu_key = 'COPY_NUMBER'
			plo_tp_key = 'PLOIDY_TYPE'


			if acn_key in header:
				acn_indx = header.index(acn_key)
				del header[acn_indx]
				header.insert(0,acn_key)

			if ass_key in header:
				ass_indx = header.index(ass_key)
				del header[ass_indx]
				header.insert(0,ass_key)

			if com_key in header:
				com_indx = header.index(com_key)
				del header[com_indx]
				header.insert(0,com_key)

			if pur_key in header:
				pur_indx = header.index(pur_key)
				del header[pur_indx]
				header.insert(0,pur_key)

			if plo_key in header:
				plo_indx = header.index(plo_key)
				del header[plo_indx]
				header.insert(0,plo_key)

			if copy_nu_key in header:
				copy_nu_indx = header.index(copy_nu_key)
				del header[copy_nu_indx]
				header.insert(0,copy_nu_key)
			
			if plo_tp_key in header:
				plo_tp_indx = header.index(plo_tp_key)
				del header[plo_tp_indx]
				header.insert(0,plo_tp_key)
			
			del header[header.index('gene')]
			header.append('gene')
			header = generate_headers_ngx_table(header)

			new_keys = {
				acn_key: {'key': acn_key, 'title': 'ABSOLUTE_COPY_NUMBER'},
				ass_key: {'key': ass_key, 'title': 'ASSESSMENT'},
				com_key :  {'key': com_key, 'title': 'COMMENT'},
				pur_key :  {'key': pur_key, 'title': 'PURITY'},
				plo_key :  {'key': plo_key, 'title': 'PLOIDY'},
				copy_nu_key :  {'key': copy_nu_key, 'title': 'COPY_NUMBER'},
				plo_tp_key :  {'key': plo_tp_key, 'title': 'PLOIDY_TYPE'}
			}

			for idx,value in enumerate(new_keys):
				n_key = [item for item in header if item.get('key')==value]
				if(not n_key):
					header.insert(0, new_keys[value])
			
			new_header = []
			for i,value in enumerate(header):
				key_name = value['key']
				if(key_name not in hide_header):
					new_header.append(value)

			return {'header': new_header, 'data': data, 'filename': save_to_cnv_file, 'status': True}, 200

		else:
			return {'header': [], 'data': [], 'filename': '', 'error': 'Invalid file', 'status': False}, 400
	else:
		return {'header': [], 'data': [], 'filename': '', 'error': 'File not found', 'status': False}, 400


def get_purecn_ctdna(project_path, sample_id, capture_id):
	file_path = project_path + '/' + sample_id + '/' + capture_id + '/purecn/'
	status = True if os.path.exists(file_path) and len(os.listdir(file_path)) > 0 else False
	if status:
		regex = '[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+.csv$'
		csv_filename = file_path + list(filter(lambda x: (re.match(regex, x) ),os.listdir(file_path)))[0]
		df = pd.read_csv(csv_filename)
		df.fillna('', inplace=True)
		json_data = df.to_dict(orient='records')
		return {'data': json_data, 'status': True}, 200

	return {'data': [], 'status': False}, 400

def read_purecn_variant_txt(file_path):
	df = pd.read_csv(file_path, delimiter=',')
	df['somatic_PureCN_status'] = df[[ 'SOMATIC.M0', 'SOMATIC.M1', 'SOMATIC.M2', 'SOMATIC.M3', 'SOMATIC.M4', 'SOMATIC.M5', 'SOMATIC.M6', 'SOMATIC.M7']].idxmax(axis=1)
	df['germline_PureCN_status'] = df[[ 'GERMLINE.M0', 'GERMLINE.M1', 'GERMLINE.M2', 'GERMLINE.M3', 'GERMLINE.M4', 'GERMLINE.M5', 'GERMLINE.M6', 'GERMLINE.M7']].idxmax(axis=1)
	df['max_somatic'] = df[[ 'SOMATIC.M0', 'SOMATIC.M1', 'SOMATIC.M2', 'SOMATIC.M3', 'SOMATIC.M4', 'SOMATIC.M5', 'SOMATIC.M6', 'SOMATIC.M7']].values.max(1)
	df['max_germline'] = df[[ 'GERMLINE.M0', 'GERMLINE.M1', 'GERMLINE.M2', 'GERMLINE.M3', 'GERMLINE.M4', 'GERMLINE.M5', 'GERMLINE.M6', 'GERMLINE.M7']].values.max(1)
	
	df = df.drop(['Sampleid', 'ID', 'SOMATIC.M0', 'SOMATIC.M1', 'SOMATIC.M2', 'SOMATIC.M3', 'SOMATIC.M4', 'SOMATIC.M5', 'SOMATIC.M6', 'SOMATIC.M7', 'GERMLINE.M0', 'GERMLINE.M1', 'GERMLINE.M2', 'GERMLINE.M3', 'GERMLINE.M4', 'GERMLINE.M5', 'GERMLINE.M6', 'GERMLINE.M7', 'GERMLINE.CONTHIGH', 'GERMLINE.CONTLOW', 'GERMLINE.HOMOZYGOUS', 'ML.SOMATIC', 'POSTERIOR.SOMATIC', 'ML.M', 'ML.M.SEGMENT', 'M.SEGMENT.POSTERIOR', 'M.SEGMENT.FLAGGED', 'ML.AR', 'AR', 'AR.ADJUSTED', 'MAPPING.BIAS', 'ML.LOH', 'CN.SUBCLONAL', 'CELLFRACTION', 'FLAGGED', 'log.ratio', 'depth', 'prior.somatic', 'prior.contamination', 'on.target', 'seg.id', 'gene.symbol'],axis=1, errors='ignore')
	df.rename({'chr': 'CHROM', 'start': 'START', 'end':'END'}, axis=1, inplace=True)
	
	df['compressed']=df.apply(lambda x:'%s%s%s%s' % (x['CHROM'],x['START'],x['REF'], x['ALT']),axis=1)
	return df

def update_pureCN_txt(variant_type, df, file_path, pureCN_probability, pureCN_status):
	
	if(variant_type == 'germline'):
		regex = '^(?:(?!-(CFDNA|T)).)*igvnav-input.txt$'
	elif(variant_type == 'somatic'):
		regex = '.*-(CFDNA|T)-.*igvnav-input.txt$'
		
	igv_nav_file = list(filter(lambda x: re.match(regex, x) and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))[0]
	igv_nav_file = file_path + '/' + igv_nav_file
			
	df_modified = pd.read_csv(igv_nav_file, delimiter='\t')
	df_modified['compressed']=df_modified.apply(lambda x:'%s%s%s%s' % (x['CHROM'],x['START'],x['REF'], re.sub(r"\W", "", x['ALT'])),axis=1)
	df_modified['match'] = df_modified['compressed'].isin(df['compressed'])
	
	for index, row in df_modified.iterrows():
		if(row['match']):
			df_modified.loc[index,'purecn_probability'] = df.loc[df['compressed'] == row['compressed'], pureCN_probability].values[0]
			df_modified.loc[index,'purecn_status']  = df.loc[df['compressed'] == row['compressed'], pureCN_status].values[0]
			df_modified.loc[index,'purecn_tot_copies']  = df.loc[df['compressed'] == row['compressed'], 'ML.C'].values[0]
		else:
			df_modified.loc[index,'purecn_probability'] = '-'
			df_modified.loc[index,'purecn_status']  = '-'
			df_modified.loc[index,'purecn_tot_copies']  = '-'
	
	df_modified = df_modified.drop(['compressed', 'match'],axis=1, errors='ignore')
	df_modified.to_csv(igv_nav_file,index = False, header=True, sep='\t')

def update_pureCN_somatic_germline(project_path, sample_id, capture_id, variant_type):

	root_path = project_path + '/' + sample_id + '/' + capture_id 
	
	regex = '.*-(CFDNA|T)-.*_variants.csv$'
	
	purecn_file_path = root_path +'/purecn/'
	
	purecn_file_path = purecn_file_path  +list(filter(lambda x: re.match(regex, x) and not x.startswith('.') and not x.endswith('.out'), os.listdir(purecn_file_path)))[0]
		
	try:
		df = read_purecn_variant_txt(purecn_file_path)
		if(variant_type == 'somatic'):
			update_pureCN_txt('somatic', df, root_path, 'max_somatic', 'somatic_PureCN_status')
		elif(variant_type == 'germline'):
			update_pureCN_txt('germline', df, root_path, 'max_germline', 'germline_PureCN_status')

		return {'status': True, 'data': 'Update Successfully', 'error': ''}, 200
	except Exception as e:
		return {'status': True, 'data': [], 'error': str(e)}, 400


def get_curated_json_file(project_path, project_name, sample_id, capture_id):
	"return the json format"

	file_path = project_path + '/' + sample_id + '/' + capture_id + '/MTBP/'
	status = True if os.path.exists(file_path) and len(os.listdir(file_path)) > 0 else False
	if status:
		file_name = list(filter(lambda x: x.endswith('.json') and not x.startswith('.'), os.listdir(file_path)))[0]
		json_file = file_path + file_name

		if os.path.exists(json_file):
			f = open(json_file)
			data = json.load(f)
			return {'data':data, 'file_name': file_name, 'status': True}, 200

	return {'data':[], 'file_name': '', 'status': False}, 400

def generate_curated_json(project_path, project_name, sample_id, capture_id, script_path):

	python_cmd = "python3 {}/MTBP_samplewise_json_format.py --path {} --project {} --sample {} --capture {}".format(script_path,project_path, project_name, sample_id, capture_id)

	try:
		proc = subprocess.check_output(python_cmd,shell=True,stderr=subprocess.STDOUT)
		return {'data': 'Json File Generated', 'status': True}, 200
	except subprocess.CalledProcessError as e:
		raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
		return {'data':[], 'status': False}, 400

def build_icpm_sample_details(cfdna):
	   
	sql = "SELECT ec.study_id as identifier, to_date(rf.date_birth::text, 'YYYYMMDD') as birthdate, get_hospital_name(rf.site_id) as hospital, 'oncotree' as cancer_taxonomy,  ec.cancer_type_code as cancer_code, 'primary' as tissue_source, get_tissue_name(ec.cancer_type_id, ec.cancer_type_code) as disease_name, ec.cell_fraction as pathology_ccf, ec.germline_dna  from ipcm_referral_t as rf INNER JOIN ipcm_ecrf_t as ec ON CAST(rf.cdk as VARCHAR) = ec.study_id WHERE rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' or rf.dna3 like'%{}%'".format(cfdna, cfdna, cfdna)
	res = db.session.execute(sql, bind=db.get_engine(current_app, 'ipcmLeaderboard'))
	res_data = generate_list_to_dict(res)
	return res_data[0]

def build_sample_details(project_name, cfdna):
		
	sample_data = {}
	
	sample_data["identifier"] = ""
	sample_data["birthdate"] = ""
	sample_data["hospital"] = ""
	sample_data["disease_name"] = "Prostate Adenocarcinoma" if project_name == "PROBIO" else ''

	table_name = 'psff_bloodreferrals' if project_name != "PROBIO" else 'probio_bloodreferrals'

	if(project_name != "PROBIO"):
		sql = "SELECT datum, rid, tid, cdk from psff_bloodreferrals WHERE blood1 like '%{}%' OR blood2 like '%{}%' or blood3 like '%{}%' or blood4 like '%{}%'".format(cfdna, cfdna, cfdna, cfdna)
		res = db.session.execute(sql)
		res_data = generate_list_to_dict(res)
		if(res_data):
			tid = res_data[0]["tid"]
			cdk = res_data[0]["cdk"].strip()
			sample_data["identifier"] = tid if tid != "" else cdk
	else:
		sql = "SELECT pnr, datum, rid, tid, cdk from probio_bloodreferrals WHERE cf_dna1 like '%{}%' OR cf_dna2 like '%{}%' or cf_dna3 like '%{}%' or kommentar like '%{}%'".format(cfdna, cfdna, cfdna, cfdna)

		res = db.session.execute(sql)
		res_data = generate_list_to_dict(res)
		
		if(res_data):
			tid = res_data[0]["tid"]
			cdk = res_data[0]["cdk"].strip()
			sample_data["identifier"] = tid if tid != "" else cdk
			pnr = res_data[0]["pnr"][0:8]
			sample_data["birthdate"] = pnr
			
			sql1 = "SELECT subject_id, CAST(dob as VARCHAR), site_name from sample_status_t WHERE pnr like '%{}%'".format(pnr)
			res1 = db.session.execute(sql1, bind=db.get_engine(current_app, 'leaderboard'))
			glb_data_1 = generate_list_to_dict(res1)
			if(glb_data_1):
				sample_data["hospital"] = glb_data_1[0]["site_name"]
		else:
			cfdna_rid = re.sub("^0+(?!$)", "", cfdna)
			sql2 = "SELECT DISTINCT rid, subjectid from biobank_t WHERE regexp_replace(referenceID, '[^a-zA-Z0-9]+', '','g') like '{}' or regexp_replace(referenceID, '[^a-zA-Z0-9]+', '','g') like '{}'  or regexp_replace(referenceID, '[^a-zA-Z0-9]+', '','g') IN ('{}') or regexp_replace(rid, '[^a-zA-Z0-9]+', '','g') IN('{}')".format(cfdna, cfdna, cfdna, cfdna_rid)
			res2 = db.session.execute(sql2, bind=db.get_engine(current_app, 'leaderboard'))
			bio_data = generate_list_to_dict(res2)
			if(bio_data):
				subject_id = bio_data[0]['subjectid'] 
				sql3 = "SELECT subject_id, CAST(dob as VARCHAR), site_name from sample_status_t WHERE regexp_replace(subject_id, '[^a-zA-Z0-9]+', '','g') like regexp_replace('{}', 'P-', '','g') or regexp_replace(subject_id, '[^a-zA-Z0-9]+', '','g') like regexp_replace('{}', '-', '','g')".format(subject_id, subject_id)
				res3 = db.session.execute(sql3, bind=db.get_engine(current_app, 'leaderboard'))
				glb_data_2 = generate_list_to_dict(res3)
			
				if(glb_data_2):
					sample_data["birthdate"] = glb_data_2[0]["dob"]
					sample_data["hospital"] = glb_data_2[0]["site_name"]

	return sample_data

def fetch_patient_info(project_name, sample_id, capture_id):

	capture_id_arr = capture_id.split("-")
	cfdna1 = capture_id_arr[4]
	cfdna =  cfdna1 if cfdna1.isdigit() else capture_id_arr[10]

	capture_format = capture_id.split("-")[0]

	try:
		if(project_name == "ICPM" or capture_format == "iPCM"):
			sample_details_json = build_icpm_sample_details(cfdna)
		else:
			sample_details_json = build_sample_details(project_name, cfdna)

		return {'status': True, 'data': sample_details_json, 'error': ''}, 200

	except subprocess.CalledProcessError as e:
		raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
		return {'data':[], 'status': False}, 400

def generate_curated_pdf(project_path, project_name, sample_id, capture_id, script_path):

	python_cmd = "python3 {}/build_base_html.py --path {} --project {} --sample {} --capture {}".format(script_path,project_path, project_name, sample_id, capture_id)
	print(python_cmd)
	try:
		proc = subprocess.check_output(python_cmd,shell=True,stderr=subprocess.STDOUT)
		return {'data': 'PDF File Generated', 'status': True}, 200
	except subprocess.CalledProcessError as e:
		raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
		return {'data':[], 'status': False}, 400


def fetch_curated_pdf(project_path, project_name, sample_id, capture_id):

	file_path = project_path + '/' + sample_id + '/' + capture_id + '/pdf/'
	temp_pdf_url_list = []
	ip_addr = 'localhost' if '5000' in request.host else request.host
	port_no = ':4200' if 'localhost' in ip_addr else ''

	status = True if os.path.exists(file_path) and len(os.listdir(file_path)) > 0 else False
	if status:
		for each_file in filter(lambda x: x.endswith('.pdf') and not x.startswith('.'), os.listdir(file_path)):
			link = 'http://'+ip_addr+port_no+'/'+ 'api/franken/viewPdf?project_name=' + project_name + '&sdid=' + sample_id + '&capture_id=' + capture_id + '&pdf_name=' + each_file
			temp_pdf_url_list.append(link)

		if len(temp_pdf_url_list) > 0:
			pdf_url_list = temp_pdf_url_list
			return {'pdf_url': pdf_url_list, 'status': True}, 200

	return {'pdf_url':[], 'status': False}, 400


def get_pdf_file(project_path, sample_id, capture_id, pdf_name):
	"""retrun the pdf file"""

	file_path = project_path + '/' + sample_id + '/' + capture_id + '/pdf/' + pdf_name
	
	if os.path.exists(file_path):
		return file_path, 200

	return file_path, 400


def get_pdf_file2(project_path, sample_id, capture_id):
	"""retrun the pdf file"""

	file_path = project_path + '/' + sample_id + '/' + capture_id + '/pdf'

	pdf_name = list(filter(lambda x: x.endswith('.pdf') and not x.startswith('.'), os.listdir(file_path)))[0]
	if os.path.exists(file_path):
		file_path = file_path + '/' + pdf_name
		return file_path, 200

	return file_path, 400