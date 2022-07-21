#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import psycopg2
from sqlalchemy import create_engine
from sqlalchemy.types import VARCHAR, Integer
import os
import yaml


def path():
	return os.path.dirname(os.path.realpath(__file__))

def readConfig(section):

	filename = path()+"/config.yml"

	with open(filename, "r") as ymlfile:
		cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
		section = cfg[section]
	return section

def read_txt_convert_dataframe(file_name, column_arr):
    df_txt = pd.read_csv(file_name, sep='\t')
    df_dunique = df_txt.drop_duplicates()
    df_final = df_dunique[column_arr]
    return df_final

def merge_hotspot_data():
    
    indel_maf_file = 'hotspotsINDELS_mod.txt'
    indel_column = ['Hugo_Symbol', 'HGVSp_Short', 'Protein_position', 'start_AA', 'end_AA']
    indel_maf_data = read_txt_convert_dataframe(indel_maf_file, indel_column)
    
    snp_maf_file = 'hotspotsMAF_mod.txt'
    snp_column = ['Hugo_Symbol', 'HGVSp_Short']
    snp_maf_data = read_txt_convert_dataframe(snp_maf_file, snp_column)
    
    df_merge = pd.merge(indel_maf_data, snp_maf_data, on=['Hugo_Symbol', 'HGVSp_Short'], how='outer')
    
    column_dict = {"Hugo_Symbol" : "gene", "HGVSp_Short" : "hgvsp", "Protein_position": "protein_position", "start_AA": "start_aa", "end_AA": "end_aa"}
    df_merge.index = np.arange(1, len(df_merge)+1)
    df_merge = df_merge.rename(columns=column_dict)
    df_merge.index.name='hs_id'
    return df_merge


def main():
    final_dump = merge_hotspot_data()
    
    ### Create the table based on the final dump
    
    if(len(final_dump) > 0 ):
        dbconfig_param = readConfig('DB')
        db_engine = dbconfig_param["curation"]

        dtypedict = {
          'hs_id': Integer,
          'gene': VARCHAR(length=255),
          'hgvsp' : VARCHAR(length=255),
          'protein_position' : VARCHAR(length=255),
          'start_aa' : VARCHAR(length=255),
          'end_aa' : VARCHAR(length=255)
        }

        alchemyEngine           = create_engine(db_engine, pool_recycle=3600)
        postgreSQLConnection    = alchemyEngine.connect()
        postgreSQLTable         = "cancer_hotspot_summary"

        frame = final_dump.to_sql(postgreSQLTable, postgreSQLConnection, if_exists='replace',  index= True, dtype=dtypedict)

        with alchemyEngine.connect() as con:
            con.execute('ALTER TABLE {} ADD PRIMARY KEY (hs_id);'.format(postgreSQLTable))

        print("Cancer Hotspot Table Created Successfully")
    else:
        print("Cancer Hotspot Dataframe was empty")


if __name__ == "__main__":
    main()