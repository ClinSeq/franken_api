import os 
import pandas as pd
import json
import re
import logging
from flask import current_app
from datetime import date, datetime
from bs4 import BeautifulSoup
import subprocess
import math
from decimal import Decimal

from franken_api.database import db
from tenacity import retry, stop_after_attempt, wait_fixed # type: ignore
from sqlalchemy.exc import SQLAlchemyError # type: ignore

# ### Run command using subprocess
def subprocess_cmd(command):
	process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
	proc_stdout = process.communicate()[0].strip()
	for line in proc_stdout.decode().split('\n'):
		logging.info("Subprocess : {}".format(line))

def json_serial(obj):
	if isinstance(obj, (datetime, date)):
		return obj.isoformat()
	elif isinstance(obj, Decimal):
		return float(obj)
	raise TypeError("Type %s not serializable" % type(obj))

@retry(stop=stop_after_attempt(3), wait=wait_fixed(2)) 
def create_db_session(db_name, query):
	session = None
	try:
		engine = db.get_engine(current_app, db_name)
		session = db.create_scoped_session(options={'bind': engine}) 
		res = session.execute(query)
		session.commit()
		return res
	except SQLAlchemyError as e:
		if session:
			session.rollback()
		raise
	finally:
		if session:
			session.remove()

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
def build_genomic_profile_sample_details(project_name, sample_id, capture_id):

	hospital_lookup = { "Karolinska": "KS", "Karolinska Sjukhuset": "KS", "Södersjukhuset": "SO", "St Göran": "ST" }

	sql = "SELECT study_code, study_site, dob, disease FROM genomic_profile_summary where project_name='{}' and sample_id='{}' and capture_id='{}'".format(project_name, sample_id, capture_id)
	res_data = create_db_session('curation', sql)
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

def build_specimen(specimen_assay):

	specimen_assay_html = ''
	specimen_assay = specimen_assay.replace("\'", "\"")
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

	return 	specimen_assay_html

def build_genomic_feature(genomic_wide):

	genome_wide_html = ''
	genomic_wide = genomic_wide.replace("\'", "\"")
	genome_wide_json = json.loads(genomic_wide)

	genome_wide_txt = ''

	tumörcellsandel = 'None'
	tumörmutationsbörda = 'None'
	msi_status = 'None'
	potentiellt_ärftliga = 'None'

	for j in genome_wide_json:
		
		if 'title' in genome_wide_json[j]:
			genome_title = genome_wide_json[j]["title"] if 'title' in genome_wide_json[j] else genome_wide_json[j]
			genome_result = genome_wide_json[j]["result"] if 'result' in  genome_wide_json[j] else ''
			result_data = genome_result if genome_result != [] else '-'

			assessment =''
			if 'assessment' in genome_wide_json[j]:
				assessment = 'Yes' if genome_wide_json[j]["assessment"]== 'Possible' else ( 'No' if genome_wide_json[j]["assessment"] == 'Not possible' else '')

			genome_wide_html += '<tr>'
			genome_wide_html += '<th>'+genome_title+'</th>'
			genome_wide_html += '<td>'+str(result_data)+'</td>'	
			genome_wide_html += '<td>'+assessment+'</td>'

			if ('assessment' in genome_wide_json[j]) and (genome_wide_json[j]["title"] == 'OTHER GENOMIC PHENOTYPE' and genome_wide_json[j]["assessment"] == 'Possible' and result_data == 'Yes'):
				genome_wide_html += '<td>'+genome_wide_json[j]["categories"]+'</td>'
			else:
				genome_wide_html += '<td>'+genome_wide_json[j]["comment"]+'</td>'

			genome_wide_html += '</tr>'

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
			genome_wide_html += '<tr>'
			genome_wide_html += '<th>'+j.upper()+'</th>'
			genome_wide_html += '<td colspan="3">'+str(genome_wide_json[j])+'</td>'	
			genome_wide_html += '</tr>'

	genome_wide_txt += '\nSEKVENSERAT MATERIAL:\n'
	genome_wide_txt += 'TUMÖRCELLSANDEL: \t{}\n'.format(tumörcellsandel)
	genome_wide_txt += 'TUMÖRMUTATIONSBÖRDA (TMB): \t{}\n'.format(tumörmutationsbörda)
	genome_wide_txt += 'MMR/MSI-status: \t{}\n'.format(msi_status)
	genome_wide_txt += 'KONSTITUTIONELLA (POTENTIELLT ÄRFTLIGA) VARIANTER: \t{}\n\n'.format(potentiellt_ärftliga)

	return genome_wide_html, genome_wide_txt


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


def build_clinical_txt(project_path, study_code, report_txt):

	file_name = os.path.join(project_path, "autoseq_clinical_report_"+study_code+".txt")

	report_txt += "\n"
	report_txt += "\nESCAT TIER I - IV: A framework to rank genomic alterations as targets for cancer precision medicine: the ESMO Scale for Clinical Actionability of molecular Targets (ESCAT). Mateo J, et al. Ann Oncol. 2018 Sep 1;29(9):1895-1902. doi: 10.1093/annonc/mdy263. PMID: 30137196; PMCID: PMC6158764."
	report_txt += "\nAMP TIER I - IV: Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists. Li et al., J Mol Diagn. 2017 Jan;19(1):4-23. doi: 10.1016/j.jmoldx.2016.10.002. PMID: 27993330; PMCID: PMC5707196."
	
	final_report_txt = report_txt

	with open(file_name, "w", encoding = 'utf-8') as f:
		f.write(final_report_txt)

	return 1


def build_html(root_path, output_pdf_path, base_html_path, html_root_path, project_name, sample_id, capture_id, specimen, seq_date, report_date, epm_dnr_data_html):
	
	report_txt = ''
	
	sql = "SELECT study_code, study_site, dob, disease, specimen_assay, genome_wide, summary_txt from genomic_profile_summary WHERE sample_id='{}' and capture_id='{}' and project_name='{}' order by id desc limit 1".format(sample_id, capture_id, project_name)
	res_sql = create_db_session('curation', sql)
	genomic_json = generate_list_to_dict(res_sql)
	
	patient_info_html = ''
	specimen_html = ''
	genome_wide_html=''

	study_code = genomic_json[0]['study_code']
	study_site = genomic_json[0]['study_site']
	dob = genomic_json[0]['dob']
	disease = genomic_json[0]['disease']

	summary_txt = genomic_json[0]['summary_txt']

	patient_info_html += '<tr><th>STUDY ID</th><td>{}</td></tr>'.format(study_code)
	patient_info_html += '<tr><th>PERSONAL NUMBER</th><td>{}</td></tr>'.format(dob)
	patient_info_html += '<tr><th>DISEASE</th><td>{}</td></tr>'.format(disease)
	patient_info_html += '<tr><th>HOSPITAL</th><td>{}</td></tr>'.format(study_site)

	### Generate a clinical report txt 
	report_txt = 'Genomisk karaktärisering av solid cancer - Implementation of Personalized Cancer Medicine (iPCM)\n'
	report_txt += 'STUDIENUMMER\t{}\nTUMÖRTYP\t{}\n'.format(study_code, disease)
	
	specimen_assay = genomic_json[0]['specimen_assay']
	specimen_html = build_specimen(specimen_assay)
	if(specimen_html == ""):
		specimen_html = '<tr><td class="no-data" colspan="3">No relevant variants detected</td></tr>'

	genome_wide = genomic_json[0]['genome_wide']
	genome_wide_html, genome_wide_txt = build_genomic_feature(genome_wide)
	report_txt += genome_wide_txt
	if(genome_wide_html == ""):
		genome_wide_html = '<tr><td class="no-data" colspan="3">No relevant variants detected</td></tr>'

	small_variant_html, smt_variant_txt= build_small_variants(root_path)
	report_txt += smt_variant_txt
	if(small_variant_html == ""):
		small_variant_html = '<tr><td class="no-data" colspan="10">No relevant variants detected</td></tr>'

	svs_html, svs_variant_txt= build_svs(root_path)
	report_txt += svs_variant_txt
	if(svs_html == ""):
		svs_html = '<tr><td class="no-data" colspan="7">No relevant variants detected</td></tr>'
	
	cnv_html, cnv_variant_txt= build_cnv(root_path)
	report_txt += cnv_variant_txt
	if(cnv_html == ""):
		cnv_html = '<tr><td class="no-data" colspan="5">No relevant variants detected</td></tr>'
	
	tech_valid_html, tech_header_html = build_tech_val_QC(root_path, project_name, capture_id)
	if(tech_valid_html == ""):
		tech_valid_html = '<tr><td class="no-data" colspan="3">No relevant variants detected</td></tr>'

	build_clinical_txt(output_pdf_path, study_code, report_txt)

	header_txt = '<li><b>SPECIMENS</b><span>'+specimen+'</span></li>'
	header_txt += '<li><b>EPM DNR</b><span>'+epm_dnr_data_html+'</span></li>'
	header_txt += '<li><b>REPORT DATE</b><span>'+report_date+'</span></li>'
	header_txt += '<li><b>RECIEVED DATE</b><span>'+str(seq_date)+'</span></li>'

	sub_header_txt = '<li><b>REPORT DATE</b><span>'+report_date+'</span></li>'

	new_text = ""
	## Read a base html and replace the curated text based on ids
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

		if 'PROBIO' in project_name or 'CLINPROST' in project_name:
			for images in soup.find_all(id='probio_logo_img'):
				img_path = html_root_path +'/static/img/logo/probio_logo.png'
			images['src'] = images['src'].replace("#", img_path)
		else:
			for images in soup.find_all(id='ki_logo_img'):
				img_path = html_root_path +'/static/img/logo/KS_SE.svg'
			images['src'] = images['src'].replace("#", img_path)

		for images in soup.find_all(id='logo_img'):
			if 'PROBIO' in project_name or 'CLINPROST' in project_name:
				img_path = html_root_path +'/static/img/logo/probio_logo.png'
			else:
				img_path = html_root_path +'/static/img/logo/KS_logo.svg'
			images['src'] = images['src'].replace("#", img_path)

		for tag in soup.find_all(id='patient_info_table_data'):
			tag.string.replace_with(patient_info_html)

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

		for tag in soup.find_all(id='tech_val_table_header'):
			tag.string.replace_with(tech_header_html)

		for tag in soup.find_all(id='tech_val_table_data'):
			tag.string.replace_with(tech_valid_html)

		for tag in soup.find_all(id='summary_notes_data'):
			tag.string.replace_with(summary_txt)

		new_text = soup.prettify(formatter=None)


	return new_text
	
def generate_pdf(project_path, project_name, sample_id, capture_id):
	logging.info("--- PDF Started ---")
	
	root_path = os.path.join(project_path,sample_id,capture_id)

	output_pdf_path =  os.path.join(root_path, "pdf")

	capture_arr = capture_id.split("-")

	html_root_path = os.path.join(current_app.root_path, 'templates/pdf')
	
	## Check the sample is 'PN' and 'C4' design
	KN_value = capture_arr[6]

	if project_name == 'WGS' or project_name == 'SARCOMA_PROSP':
			base_html_path = os.path.join(html_root_path,'base_wgs.html')
	else:
		if "PN" in KN_value:
			base_html_path = os.path.join(html_root_path,'base_pn.html')
		else:
			base_html_path = os.path.join(html_root_path,'base.html')
	
	## Get the Normal value from the capture ID 
	normal_idx = capture_arr.index("N")
	normal_cfdna = capture_arr[normal_idx+1]

	## Get the CFDNA Value from the capture ID
	cfdna_idx = capture_arr.index("T") if 'T' in capture_arr else capture_arr.index("CFDNA")
	cfdna = capture_arr[cfdna_idx+1]

	## Create a PDF Folder and copy pdf template 
	base_html_file = output_pdf_path+"/base_"+project_name+"_"+sample_id+"_"+cfdna+".html"
	pdf_file_name = output_pdf_path+"/report_"+project_name+"_"+sample_id+"_"+cfdna+".pdf"
	log_name = output_pdf_path+"/log_"+project_name+"_"+sample_id+"_"+cfdna+".log"

	if(not os.path.exists(output_pdf_path)):
		os.mkdir(output_pdf_path)
	else:
		for f in os.listdir(output_pdf_path):
			os.remove(os.path.join(output_pdf_path, f))

	report_date = date.today().strftime('%Y-%m-%d')

	if 'PROBIO' in project_name:
		epm_dnr_data_html = "2016/101-32"
	else:
		epm_dnr_data_html = "2021-00135"

	if project_name == 'WGS' or project_name == 'SARCOMA_PROSP':
		specimen =  'Fresh Frozen'
	else:
		specimen =  'CFDNA' if 'CFDNA' in capture_arr else ( 'FFPE' if 'T' in capture_arr else '')

	seq_date_str = re.findall(r'\d+', capture_arr[5])[0][:8]
	seq_date = datetime.strptime(seq_date_str, "%Y%m%d").date().strftime("%Y-%m-%d")

	
	logging.basicConfig(format = '%(asctime)s  %(levelname)-10s %(name)s %(message)s', level=logging.INFO , filename=log_name, datefmt =  "%Y-%m-%d %H:%M:%S")
	logging.info('--- Generated PDF format Started---')

	logging.info("Sample Id : {} || Capture Id : {} || Outpue File Name : {} ".format(sample_id,capture_id, base_html_file))

	try:
		
		html_text = build_html(root_path, output_pdf_path, base_html_path, html_root_path, project_name, sample_id, capture_id, specimen, seq_date, report_date, epm_dnr_data_html)

		if html_text:
		 	# Create a new html based on the base template
			with open(base_html_file, "w", encoding = 'utf-8') as file:
				file.write(str(html_text))
			logging.info("---- PDF Generated ----")
			return html_text, pdf_file_name
		

	except Exception as e:
		logging.error("Failed : {}".format(str(e)))
		logging.error('--- Generated PDF format Failed ---\n')
		raise