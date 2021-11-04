#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import psycopg2
from sqlalchemy import create_engine
import numpy as np
import json
import os 
import glob
import re
import argparse
from configparser import ConfigParser
import yaml


def path():
	return os.path.dirname(os.path.realpath(__file__))

def readConfig(section):
	root_path = path()
	filename = "config.yml"
	with open(filename, "r") as ymlfile:
		cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
		section = cfg[section]
	return section

def readCNVOutFile(file_name):
	cnv_df = pd.read_csv(file_name, delimiter = "\t")
	ploidy_in_dataframe = "PLOIDY_TYPE" in cnv_df
	
	if(ploidy_in_dataframe):
		cnv_df = cnv_df[(cnv_df['PLOIDY_TYPE'] == "Haploid") | (cnv_df['PLOIDY_TYPE'] == "Diploid")]
		cnv_json = cnv_df.to_json(orient="records")
		cnv_res = json.loads(cnv_json)
	else:
		cnv_res = []
		
	return cnv_res

def readSVSOutFile(file_name):
	svs_df = pd.read_csv(file_name, delimiter = "\t")
	call_in_dataframe = "CALL" in svs_df
	
	if call_in_dataframe:
		svs_df = svs_df[(svs_df['CALL'] == True)]
		svs_json = svs_df.to_json(orient="records")
		svs_res = json.loads(svs_json)
	else:
		svs_res = []
		
	return  svs_res

def readOutFile(file_name, condition_txt):
	out_df = pd.read_csv(file_name, delimiter = "\t")
	call_in_dataframe = "CALL" in out_df
	
	if call_in_dataframe:
		out_df = out_df[(out_df['CALL'] == condition_txt)]
		out_df = out_df.replace(np.nan, '', regex=True)
		out_json = out_df.to_json(orient="records")
		out_res = json.loads(out_json)
	else:
		out_res = []
		
	return out_res

def main(nfs_path, project_name):

	#nfs_path = '/sdata/PROBIO/autoseq-output/'

	ignored = {"._.DS_Store", ".DS_Store", ".nohup.log"}
	sdid_string = [x for x in os.listdir(nfs_path) if x not in ignored]

	profile_list = []
	col_id = 1

	config_param = readConfig('PROJECT')
	project_pattern = config_param[project_name.capitalize()]

	for s in sdid_string:

		if(re.match('^[^.]+$', s)):
			sample_id = s
			sample_path = os.path.join(nfs_path, sample_id) 
			
			if(os.path.isdir(sample_path)):
				ignored = {"._.DS_Store", ".DS_Store", ".nohup.log"}
				folders = [x for x in os.listdir(sample_path) if x not in ignored]

				for fold in folders:
					if re.search(project_pattern, fold):
						profile_dict = {}
						capture_id = fold
						profile_dict['id'] = col_id
						profile_dict['sample_id'] = sample_id
						profile_dict['capture_id'] = capture_id
						profile_dict['basic_qc'] = ''
						profile_dict['franken_plot'] = ''
						profile_dict['ploidy'] = ''
						profile_dict['ctdna_fraction'] = ''
						profile_dict['ctdna_param'] = []
						profile_dict['ctdna_method'] = ''
						profile_dict['ctdna_category'] = ''
						profile_dict['study_code'] = ''
						profile_dict['study_site'] = ''
						profile_dict['comment_info'] = ''

						print("============== {} =======================".format(sample_id))

						src_path = os.path.join(sample_path, capture_id)

						files_arr = glob.glob(src_path + "/*-igvnav-input.txt", recursive = True)
						svs_files_arr = glob.glob(src_path + "/svs/igv/*sv-annotated.txt", recursive = True)
						cnv_files_arr = glob.glob(src_path + "/cnv/*_curated.cns", recursive = True)

						svs_variants_txt = ''
						if(svs_files_arr):
							svs_res = readSVSOutFile(svs_files_arr[0])
							if svs_res:
								for ss in svs_res:
									if (ss['CALL'] == True or ss['CALL'] == 'True'):
										svs_variants_txt += ss['SVTYPE'] + ', ' + ss['IGV_COORD'] + ', sv length:' + ss['SV_LENGTH'] + ', supporting reads: ' + ss['SUPPORT_READS'] + ', tool: ' + ss['TOOL'] + ', sample: ' + ss['SAMPLE'] + ', GeneA: ' + ss['GeneA'] + ', GeneB: ' + ss['GeneB'] + '; \n';
								print("SVS Done")
							else:
								svs_variants_txt = "NA"
								print("Empty SVS")
						else:
							print("File Not Found SVS ")

						profile_dict['structural_variants'] = svs_variants_txt

						germline_mut_txt = ''
						somatic_arr = {}
						if(files_arr):
							for f in files_arr:
								file_name = f.split('/')[-1]
								if re.search('-CFDNA-', file_name):
									somatic_res = readOutFile(f, 'S')
									if somatic_res:
										
										for sm in somatic_res:
											hotspot_cont = (sm['HGVSp'] if(sm['HOTSPOT'] == '') else 'hotspot mutation ('+sm['HGVSp']+') ') if 'HOTSPOT' in sm else sm['HGVSp']

											chrom = sm['Chromosome'] if 'Chromosome' in sm else sm['CHROM']
											start = sm['Start'] if 'Start' in sm else sm['START']
											end = sm['Stop'] if 'Stop' in sm else sm['END']

											somatic_mut_txt = sm['GENE'] + ' (chr' + str(chrom) + ':' + str(start) + '-' + str(end) + '), ' + hotspot_cont + ', ' + sm['IMPACT'] + ' impact (' + sm['CONSEQUENCE'] + ')' + '; \n ';

											if(sm['CLONALITY'] not in somatic_arr):
												somatic_arr[sm['CLONALITY']] = somatic_mut_txt;
											else:
												somatic_arr[sm['CLONALITY']] = somatic_arr[sm['CLONALITY']] + somatic_mut_txt;

										print("Somatic Done")
									else:
										somatic_arr = "NA"
										print("Empty Somatic")                            
								else:
									germline_res = readOutFile(f, 'G')
									if germline_res:
										for gm in germline_res:
											rsid = gm['RSID'] if 'RSID' in gm else ''
											chrom = gm['Chromosome'] if 'Chromosome' in gm else gm['CHROM']
											start = gm['Start'] if 'Start' in gm else gm['START']
											end = gm['Stop'] if 'Stop' in gm else gm['END']

											germline_mut_txt += gm['GENE'] + ' (chr' + str(chrom) + ':' + str(start) + '-' + str(end) + '), ' + ', ' + gm['HGVSp'] + ', ' + gm['IMPACT'] + ' impact (' + gm['CONSEQUENCE'] + '),  VAF ' + str(gm['N_VAF']) + ', GNOMAD ' + str(gm['gnomAD']) + ',' + gm['CLIN_SIG'] + ',' + rsid + '; \n';
										print("Germline Done")
									else:
										germline_mut_txt = "NA"
										print("Empty Germline")
						else:
							print("File Not Found for Germline and Somatic")

						profile_dict['somatic_mutations'] =  "NA" if somatic_arr == "" else str(somatic_arr);
						profile_dict['germline_alterations'] = germline_mut_txt

						cnvs_txt = ""
						cnvs_arr = {}
						cnvs_type = ""
						cnvs_arr["somatic"] = "NA"
						cnvs_arr["germline"] = "NA"
						regex_germline = '^(\W.*)(_germline_curated.cns)$'
						if(cnv_files_arr):
							for cnv_f in cnv_files_arr:
								if(re.match(regex_germline, cnv_f)):
									cnvs_type = "germline"
								else:
									cnvs_type = "somatic"

								cnvs_arr[cnvs_type] = ""

								cnv_res = readCNVOutFile(cnv_f)
								if cnv_res:
									for cvs in cnv_res:

										assessment = "";
										if ("ASSESSMENT" in cvs):
											assessment = cvs['ASSESSMENT'] +', ' if cvs['ASSESSMENT'] != 'undefined' or cvs['ASSESSMENT'] != '' else ''

										cnvs_txt = assessment +cvs['PLOIDY_TYPE']+', Ploidy : '+str(cvs['PLOIDY'])+',- PURITY : '+str(cvs['PURITY'])+ ', COPY_NUMBER : ' + str(cvs['COPY_NUMBER'])+ ', genes in segment : ' + cvs['gene']+ '; \n';
										cnvs_arr[cnvs_type] = cnvs_arr[cnvs_type] + cnvs_txt;

									print("CNV {} Done".format(cnvs_type))

								else:
									cnvs_arr[cnvs_type] = "NA"
									print("Empty CNV {}".format(cnvs_type))
						else:
							print("File Not Found CNVs ")

						profile_dict['cnvs'] = str(cnvs_arr)


						profile_list.append(profile_dict)
						col_id += 1
				print('\n')
			else:
				print("Not found : ",sample_id)

	#json_data = json.dumps(profile_list)
	#print(json_data, "\n")

	df = pd.DataFrame.from_dict(profile_list, orient='columns')

	dbconfig_param = readConfig('DB')
	db_engine = dbconfig_param["curation"]

	alchemyEngine           = create_engine(db_engine, pool_recycle=3600)
	postgreSQLConnection    = alchemyEngine.connect()
	postgreSQLTable         = project_name.lower()+"_summary"

	try:
		frame = df.to_sql(postgreSQLTable, postgreSQLConnection, if_exists='replace',  index = False)
	except ValueError as vx:
		print("ValueError:",vx)
	except Exception as ex:  
		print("Execption:",ex)
	else:
		print("PostgreSQL Table %s has been created successfully."%postgreSQLTable)
	finally:
		postgreSQLConnection.close()


if __name__ == "__main__":

	# Create the parser
	profile_parser = argparse.ArgumentParser(description='Generate Profile Summary for PROBIO and PSFF')

	# Add the arguments
	profile_parser.add_argument('folder_path', metavar='nfs root path', type=str, help='define the sample autoseq-output path')
	profile_parser.add_argument('project_name', metavar='project name', type=str, help='define the project name (eg: PROBIO, PSFF )')

	args = profile_parser.parse_args()

	nfs_path = args.folder_path
	project_name = args.project_name

	if not os.path.isdir(nfs_path):
		print('The path specified does not exist')
		sys.exit()
	else:
		main(nfs_path,project_name)