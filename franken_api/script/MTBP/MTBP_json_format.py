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


def readConfig(section,filename="config.yml"):
    with open(filename, "r") as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
        section = cfg[section]
    return section


### Fetch SQL Query 
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


def build_icpm_sample_details(cfdna):
       
    sql = "SELECT ec.study_id as identifier, to_date(rf.datum::text, 'YYYYMMDD') as sample_date, to_date(rf.date_birth::text, 'YYYYMMDD') as birthdate, get_hospital_name(rf.site_id) as hospital, 'oncotree' as cancer_taxonomy,  ec.cancer_type_code as cancer_code, 'primary' as tissue_source, get_tissue_name(ec.cancer_type_id, ec.cancer_type_code) as tissue_type, ec.cell_fraction as pathology_ccf, ec.germline_dna  from ipcm_referral_t as rf INNER JOIN ipcm_ecrf_t as ec ON CAST(rf.cdk as VARCHAR) = ec.study_id WHERE rf.dna1 like'%{}%' OR rf.dna2 like'%{}%' or rf.dna3 like'%{}%'".format(cfdna, cfdna, cfdna)
    res_data = fetch_sql_query('ipcmLeaderboard', sql)
    
    return res_data


def build_sample_details(cfdna):
        
    sample_data = {}
    
    sample_data["identifier"] = "NA"
    sample_data["sample_date"] = "NA"
    sample_data["birthdate"] = "NA"
    sample_data["hospital"] = "NA"

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
            
    sample_data["seq_date"] = "NA"
    sample_data["cancer_taxonomy"] = "NA"
    sample_data["cancer_code"] = "NA"
    sample_data["tissue_source"] = "NA"
    sample_data["tissue_type"] = "NA"
    sample_data["pathology_ccf"] = "NA"
    sample_data["bioinf_ccf"] = "NA"
    sample_data["germline_dna"] = "NA"
    
    return sample_data


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

            msi_list = "NA"
            if 'msing_score' in qc_df_data.columns:
                msi_list = qc_df_data['msing_score'].dropna().tolist()[0]

            qc_df_data = qc_df_data[column_list]
            qc_df_data = qc_df_data.rename(columns=column_dict)

            qc_df_data.fillna('NA', inplace=True)
            return msi_list, qc_df_data.to_json(orient = 'index')
    
    except Exception as e:
        print("Exception", str(e))

        
    return msi_list, qc_df_data


def build_ploidy(root_path):
    file_path = root_path + '/purecn/'
    
    regex = '[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+.csv'
    
    ploidy_list = "NA"
    purity_list = "NA"
    
    try:
        purecn_filename = file_path + list(filter(lambda x: (re.match(regex, x) ),os.listdir(file_path)))[0]    
    
        if len(purecn_filename) > 0:
            purecn_df_data = pd.read_csv(purecn_filename, delimiter = ",")
    
            purecn_df_data = purecn_df_data[['Sampleid', 'Purity', 'Ploidy', 'Sex', 'Contamination']]
            ploidy_list = round(purecn_df_data['Ploidy'].tolist()[0],4)
            purity_list = round(purecn_df_data['Purity'].tolist()[0],4)
            
    except Exception as e:
        print("Exception", str(e))
        
    return purity_list, ploidy_list


def build_small_variants(root_path):

    file_path = root_path + '/'
    smv_file_list = list(filter(lambda x: x.endswith('-igvnav-input.txt') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))
    regex = '^(?:(?!-(CFDNA|T)).)*igvnav-input.txt$'
    
    smv_df = pd.DataFrame()
    
    for i in smv_file_list:
        smv_filename = file_path + i
        smv_df_data = pd.read_csv(smv_filename, delimiter = "\t")
        
        if 'CALL' in smv_df_data.columns:

            column_list = ['CHROM', 'START', 'END', 'REF', 'ALT', 'N_DP', 'N_ALT', 'N_VAF']
            column_dict = {'CHROM': 'chr', 'START': 'start', 'END': 'end', 'REF' : 'ref', 'ALT' : 'alt', 'ALT' : 'alt', 'N_DP' : 'ref_reads',  'N_ALT' : 'alt_reads', 'N_VAF' : 'vaf'}
            
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
                smv_df_data['origin'] = "germline"
            else:
                smv_df_data['origin'] = "somatic"

            smv_df_data['alt'] = smv_df_data['alt'].map(lambda x: x.lstrip('[').rstrip(']'))

            smv_df_data['strand'] = '+'

        smv_df = smv_df.append(smv_df_data)

    smv_df.fillna('NA', inplace=True)
    smv_df.reset_index(drop=True, inplace=True)
    return smv_df.to_json(orient = 'index') # records -> get array dict format


def build_cnv(root_path):
    
    cnv_df = pd.DataFrame()
    file_path = root_path + '/cnv/'
    cnv_file_list = list(filter(lambda x: x.endswith('_curated.cns') and not x.startswith('.') and not x.endswith('.out'), os.listdir(file_path)))
    if len(cnv_file_list) > 0:
        regex = '^(?:(?!-(CFDNA|T)).)*_curated.cns$'
        for i in cnv_file_list:
            cnv_filename = file_path + i
            cnv_df_data = pd.read_csv(cnv_filename, delimiter = "\t")
            
            column_list = ['chromosome', 'start', 'end', 'gene', 'ASSESSMENT', 'origin']
            column_dict = {'chromosome': 'chr', 'gene': 'genes', 'ASSESSMENT': 'type'}
        
#             if 'PLOIDY_TYPE' in cnv_df_data.columns:
#                 cnv_df_data = cnv_df_data.loc[(cnv_df_data['PLOIDY_TYPE'] == "Diploid") | (cnv_df_data['PLOIDY_TYPE'] == "Haploid")]

            if 'ASSESSMENT' in cnv_df_data.columns:
                cnv_df_data = cnv_df_data.loc[(cnv_df_data['ASSESSMENT'].notnull())]

                if(re.match(regex, i)):
                    cnv_df_data['origin'] = "germline"
                else:
                    cnv_df_data['origin'] = "somatic"
                
                if 'COPY_NUMBER' in cnv_df_data.columns:
                    column_list.append('COPY_NUMBER')
                    column_dict['COPY_NUMBER'] = 'copy_number'
                else:
                    column_list.append('copy_number')
                    cnv_df_data['copy_number'] = 'NA'
                    
                cnv_df_data = cnv_df_data[column_list]
                cnv_df_data = cnv_df_data.rename(columns=column_dict)
                    
                cnv_df_data['genes'] = cnv_df_data.genes.apply(lambda x: x.split(', '))
                cnv_df = cnv_df.append(cnv_df_data)
    
    cnv_df.fillna('NA', inplace=True)
    cnv_df.reset_index(drop=True, inplace=True)
    return cnv_df.to_json(orient = 'index')
    

def build_svs(root_path):
    file_path = root_path + '/svs/igv/'
    
    svs_filter = pd.DataFrame()
    
    try:
        svs_filename = file_path + list(filter(lambda x: (re.match('[-\w]+-(CFDNA|T)-[A-Za-z0-9-]+-sv-annotated.txt', x) or x.endswith('_annotate_combined_SV.txt')) and not x.startswith('.') and not x.endswith('.out'),os.listdir(file_path)))[0]

        sample_list = ['germline', 'somatic']
        svs_filter = pd.read_csv(svs_filename, delimiter = "\t")
        
        column_list = ['SAMPLE', 'SVTYPE', 'strand', 'variant', 'vaf', 'GENE_A', 'GENE_B', 'CLONALITY', 'SECONDHIT', 'consequence', 'variant_string']
        column_dict = {'SAMPLE': 'origin', 'SVTYPE': 'sv_type', 'GENE_A': 'gene_A' , 'GENE_B' : 'gene_B', 'CLONALITY': 'clonality', 'SECONDHIT' : 'secondhit'}
        
        if 'CALL' in svs_filter.columns:
            svs_filter = svs_filter.loc[(svs_filter['CALL'] == True ) | (svs_filter['CALL'] == 'true')]
                      
            if not svs_filter.empty:
                svs_filter["variant"] = svs_filter.apply(lambda x:'%s:%s_%s' % (x['CHROM_A'],x['START_A'],x['END_A']),axis=1)
                svs_filter = svs_filter[['SAMPLE', 'SVTYPE', 'variant', 'GENE_A', 'GENE_B','SECONDHIT', 'CLONALITY']]

                svs_filter = svs_filter[svs_filter['SAMPLE'].isin(sample_list)]

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
        else:
            svs_filter = pd.DataFrame()

    except Exception as e:
        print("Exception", str(e))
        svs_filter = pd.DataFrame()
        
    return svs_filter.to_json(orient = 'index')


def build_json(root_path, file_name, project_name, cfdna):
    
    print("--- MTBP Json format Started ---\n")
    print("Path : ", root_path, "/", file_name)

    project_json = {}
    
    # Sample Information
    if(project_name == "ICPM"):
        sample_details_json = build_icpm_sample_details(cfdna)
    else:
        sample_details_json = build_sample_details(cfdna)
       
    if(sample_details_json):
        project_json["sample"] = sample_details_json
    else:
        project_json["sample"] = {"identifier": "NA",  "sample_date": "NA", "seq_date": "NA", "birthdate": "NA", "hospital": "NA", "cancer_taxonomy": "NA", "cancer_code": "NA", "tissue_source": "NA", "tissue_type": "NA", "pathology_ccf": "NA", "bioinf_ccf": "NA", "germline_dna": "NA"}
    
    # Pipeline 
    project_json["pipeline"] = { "genome_reference": "hg19", "version": "1.0"}
    
    # QC
    logging.info('--- Build QC json started ---')
    msi_val, qc_json = build_qc(root_path)
    project_json["qc"] =  "NA" if qc_json == "" else json.loads(qc_json)
    logging.info('--- QC completed ---')
    
    # Get the Ploidy value
    purity_val, ploidy_val = build_ploidy(root_path)
    
    # Phenotype
    project_json["phenotypes"] = {"purity" : purity_val ,"ploidy": ploidy_val, "msi": "LOW",  "hrrd":"NA",  "tmb": "NA"}
    
   
    # Small Variant (Somatic & Germline)
    logging.info('--- Build Small Variant json started ---')
    small_variant_json = build_small_variants(root_path)
    project_json["small_variants"] = json.loads(small_variant_json)
    logging.info('--- Small Variant completed ---')
    
    # CNVs
    logging.info('--- Build CNvs json started ---')
    cnv_json = build_cnv(root_path)
    project_json["cnas"] = json.loads(cnv_json)
    logging.info('--- CNVs completed ---')
    
    # SVS
    logging.info('--- Build SVS json started ---')
    svs_json = build_svs(root_path)
    project_json["svs"] = json.loads(svs_json)
    logging.info('--- SVS completed ---')
    
    final_json = json.dumps(project_json)
    
    with open(file_name, 'w') as f:
        json.dump(project_json, f, indent=4)
    
    logging.info('--- Generated Json format successfuly ---')
    
    print("\n----  MTBP Json format Completed -----\n")


def main(nfs_path,project_name, capture_format):
    
    regx_pat = '^P-.*[0-9]$'
    regx_capture_pat = '^{}-P-.*[0-9]$'.format(capture_format)
    
    for root, dirs, files in os.walk(nfs_path):
        for d in dirs:
            pat_check = re.match(regx_pat, d)
            if pat_check:
                dir_path = os.path.join(nfs_path,d)
                inner_dir_list = os.listdir(dir_path)
                if(len(inner_dir_list) > 0):
                    for indir in inner_dir_list:
                        capture_check = re.match(regx_capture_pat, indir)
                        if(capture_check):
                            cfdna = indir.split("-")[4]
                            root_path = os.path.join(dir_path,indir)
                            output_path = root_path+"/MTBP";
                            #project_name = indir.split("-")[0]
                            file_name = output_path+"/MTBP_"+project_name+"_"+d+"_CFDNA"+cfdna+".json"
                            log_name = output_path+"/MTBP_"+project_name+"_"+d+"_CFDNA"+cfdna+".log"
                            if(not os.path.exists(output_path)):
                                os.mkdir(output_path);

                            logging.basicConfig(format = '%(asctime)s  %(levelname)-10s %(name)s %(message)s', level=logging.INFO , filename=log_name,  filemode='w', datefmt =  "%Y-%m-%d %H:%M:%S")
                            logging.info('--- Generated Json format ---')

                            logging.info("Sample Id : {} || Capture Id : {} || Outpue File Name : {} ".format(d,indir, file_name))
                            try:
                                build_json(root_path, file_name, project_name, cfdna)
                            except Exception as e:
                                print("Exception", str(e))
                                logging.error("Failed : {}".format(str(e)))
                                logging.error('--- Generated Json format Failed ---')
                                raise
                                #break

if __name__ == "__main__":
        
    # Create the parser
    profile_parser = argparse.ArgumentParser(description='Generate MTBP Json')

    # Add the arguments
    profile_parser.add_argument('input_path', metavar='nfs root path', type=str, help='define the sample autoseq-output path')
    profile_parser.add_argument('project_name', metavar='project name', type=str, help='define the project name (eg: PROBIO, PSFF, IPCM )')
    profile_parser.add_argument('capture_format', metavar='project name', type=str, help='define the capture name (eg: PB, PSFF, iPCM )')
    
    args = profile_parser.parse_args()

    nfs_path = args.input_path
    project_name = args.project_name
    capture_format = args.capture_format

    if not os.path.isdir(nfs_path):
        print('The path specified does not exist')
        sys.exit()
    else:
        main(nfs_path,project_name, capture_format)
    