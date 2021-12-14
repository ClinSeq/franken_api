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


def main(root_path, file_name):
    
    print("--- MTBP Json format Started ---\n")
    print("Path : ", root_path, "/", file_name)

        
    project_json = {}
    
    # Sample Information
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


if __name__ == "__main__":
    
    project_path = "/sdata/RESBIO/autoseq-output/"
    regx_pat = '^P-.*[0-9]$'
    regx_capture_pat = '^RB-P-.*[0-9]$'
    
    for root, dirs, files in os.walk(project_path):
        for d in dirs:
            pat_check = re.match(regx_pat, d)
            if pat_check:
                dir_path = os.path.join(project_path,d)
                inner_dir_list = os.listdir(dir_path)
                if(len(inner_dir_list) > 0):
                    for indir in inner_dir_list:
                        capture_check = re.match(regx_capture_pat, indir)
                        if(capture_check):
                            cfdna = indir.split("-")[4]
                            root_path = os.path.join(dir_path,indir)
                            output_path = root_path+"/MTBP";
                            project_name = indir.split("-")[0]
                            file_name = output_path+"/MTBP_"+project_name+"_"+d+"_CFDNA"+cfdna+".json"
                            log_name = output_path+"/MTBP_"+project_name+"_"+d+"_CFDNA"+cfdna+".log"
                            if(not os.path.exists(output_path)):
                                os.mkdir(output_path);
                            
                            logging.basicConfig(format = '%(asctime)s  %(levelname)-10s %(name)s %(message)s', level=logging.INFO , filename=log_name,  filemode='w', datefmt =  "%Y-%m-%d %H:%M:%S")
                            logging.info('--- Generated Json format ---')

                            logging.info("Sample Id : {} || Capture Id : {} || Outpue File Name : {} ".format(d,indir, file_name))
                            try:
                                main(root_path, file_name)
                            except Exception as e:
                                print("Exception", str(e))
                                logging.error("Failed : {}".format(str(e)))
                                logging.error('--- Generated Json format Failed ---')
                                raise
                                #break
