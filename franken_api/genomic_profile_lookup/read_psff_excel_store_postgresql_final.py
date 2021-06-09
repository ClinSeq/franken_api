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

def main():

	nfs_path = '/nfs/PSFF/autoseq-output/'

	empty_file_list = []

	profile_list = []
	col_id = 1

	ignored = {"._.DS_Store", ".DS_Store", ".nohup.log"}
	sdid_string = [x for x in os.listdir(nfs_path) if x not in ignored]

	for s in sdid_string:

		if(re.match('^[^.]+$', s)):
			sample_id = s
			sample_path = os.path.join(nfs_path, sample_id) 
			
			if(os.path.isdir(sample_path)):
				ignored = {"._.DS_Store", ".DS_Store", ".nohup.log"}
				folders = [x for x in os.listdir(sample_path) if x not in ignored]

				for fold in folders:
					if re.search('PSFF-P', fold):
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
									svs_variants_txt += ss['SVTYPE'] + ', '+ss['IGV_COORD']+', '+str(ss['SV_LENGTH'])+' sv length, with '+str(ss['SUPPORT_READS'])+' reads, '+ss['TOOL']+' tool, with '+ss['SAMPLE']+' sample '
								print("SVS Done")
							else:
								print("Empty SVS")
						else:
							print("File Not Found SVS ")

						profile_dict['structural_variants'] = svs_variants_txt

						somatic_mut_txt = ''
						germline_mut_txt = ''
						if(files_arr):
							for f in files_arr:
								file_name = f.split('/')[-1]
								if re.search('-CFDNA-', file_name):
									somatic_res = readOutFile(f, 'S')

									if somatic_res:
										for sm in somatic_res:
											second_hit = ('' if(sm['SECONDHIT'] == '' or sm['SECONDHIT'] == '-') else ', Second hit with'+sm['SECONDHIT']) if 'SECONDHIT' in sm else ''
											hotspot = (sm['HGVSp'] if(sm['HOTSPOT'] == '') else 'hotspot mutation ('+sm['HGVSp']+') ') if 'HOTSPOT' in sm else sm['HGVSp']
											CLONALITY = sm['CLONALITY'] if 'CLONALITY' in sm else ''
											somatic_mut_txt += CLONALITY + ', '+sm['GENE']+' ('+str(sm['CHROM'])+':'+str(sm['START'])+'-'+str(sm['END'])+'), '+hotspot+', '+sm['IMPACT']+' impact ('+sm['CONSEQUENCE']+')'+second_hit+' '
										print("Somatic Done")
									else:
										print("Empty Somatic")                            
								else:
									germline_res = readOutFile(f, 'G')
									if germline_res:
										for gm in germline_res:
											rsid = gm['RSID'] if 'RSID' in gm else ''
											germline_mut_txt += gm['GENE']+', '+gm['HGVSp']+', '+gm['IMPACT']+' impact ('+gm['CONSEQUENCE']+'), '+rsid+' and VAF='+str(gm['N_VAF'])+' '
										print("Germline Done")
									else:
										print("Empty Germline")
						else:
							print("File Not Found for Germline and Somatic")

						profile_dict['somatic_mutations'] = somatic_mut_txt
						profile_dict['germline_alterations'] = germline_mut_txt

						cnvs_txt = ''
						if(cnv_files_arr):
							for cnv_f in cnv_files_arr:
								cnv_res = readCNVOutFile(cnv_f)
								if cnv_res:
									for cvs in cnv_res:
										cnvs_txt += 'Genomic region: ["'+ cvs['gene']+'"] ('+str(cvs['chromosome'])+':'+str(cvs['start'])+'-'+str(cvs['end'])+'), '+cvs['PLOIDY_TYPE'] + ', ' +str(cvs['COPY_NUMBER'])+', region size: '+cvs['SIZE']+','+str(cvs['log2'])+' '
									print("CNV Done")
								else:
									print("Empty CNV")
						else:
							print("File Not Found CNVs ")

						profile_dict['cnvs'] = cnvs_txt


						profile_list.append(profile_dict)
						col_id += 1
				print('\n')
			else:
				print("Not found : ",sample_id)

	df = pd.DataFrame.from_dict(profile_list, orient='columns')

	alchemyEngine           = create_engine('postgresql+psycopg2://referral_writer:ProbioWriter@127.0.0.1:5432/curation', pool_recycle=3600)
	postgreSQLConnection    = alchemyEngine.connect()
	postgreSQLTable         = "psff_summary"

	try:
		frame = df.to_sql(postgreSQLTable, postgreSQLConnection, if_exists='fail',  index = False)
	except ValueError as vx:
		print(vx)
	except Exception as ex:  
		print(ex)
	else:
		print("PostgreSQL Table %s has been created successfully."%postgreSQLTable)
	finally:
		postgreSQLConnection.close()


if __name__ == "__main__":
	main()