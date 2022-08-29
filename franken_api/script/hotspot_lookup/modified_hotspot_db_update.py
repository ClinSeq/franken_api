#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import psycopg2
from flask_sqlalchemy import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.types import VARCHAR, Integer


def read_txt_convert_dataframe(file_name):
    df_txt = pd.read_csv(file_name, sep=',')
    df_final = df_txt.drop_duplicates()
    return df_final


def merge_hotspot_data():
    indel_maf_file = 'Indel.csv'
    indel_maf_data = read_txt_convert_dataframe(indel_maf_file)
    indel_count = len(indel_maf_data)
    
    snp_maf_file = 'SNV.csv'
    snp_maf_data = read_txt_convert_dataframe(snp_maf_file)
    snp_count = len(snp_maf_data)
    
    df_concat = pd.concat([snp_maf_data, indel_maf_data],axis=0)
    
    df_concat.index = np.arange(1, len(df_concat)+1)
    df_concat.index.name='h_id'
    
    print(df_concat)
    tot_count = len(df_concat)
    
    print(" SNP Count '{}' + '{}' INDEL Count = {}".format(snp_count, indel_count,tot_count))
    
    return df_concat


def main():
    final_dump = merge_hotspot_data()
    
    ### Create the table based on the final dump
    
    if(len(final_dump) > 0 ):
        db_engine = "postgresql+psycopg2://postgres:postgres@127.0.0.1:5432/curation"

        dtypedict = {
          'h_id': Integer,
          'gene': VARCHAR(length=255),
          'hgvsp' : VARCHAR(length=255),
          'amino_acid_position' : VARCHAR(length=255),
          'start_aa' : Integer,
          'end_aa' : Integer
        }

        alchemyEngine           = create_engine(db_engine, pool_recycle=3600)
        postgreSQLConnection    = alchemyEngine.connect()
        postgreSQLTable         = "cancer_hotspot_summary"

        frame = final_dump.to_sql(postgreSQLTable, postgreSQLConnection, if_exists='replace',  index= True, dtype=dtypedict)

        with alchemyEngine.connect() as con:
            con.execute('ALTER TABLE {} ADD PRIMARY KEY (h_id);'.format(postgreSQLTable))

        print("Table Created Successfully")
    else:
        print("Cancer Hotspot Dataframe was empty")### Create the table based on the final dump


if __name__ == "__main__":
    main()
