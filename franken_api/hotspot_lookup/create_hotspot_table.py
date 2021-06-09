#!/usr/bin/env python
# coding: utf-8

### Lookup Hotspot

import time
import pandas as pd
from sqlalchemy import create_engine
import psycopg2
import os

def build_dataframe_to_table(dir_path, src_file_path, table_name, csv_file_name):

    df_data = pd.read_csv(src_file_path, sep="\t")
    df_data = df_data.reset_index()
    df_data = df_data.drop(columns =['p', 'q', 'Hotness', 'Stem_Strength', 'SS_Loop_Pos', 'SS_Loop_Len', 'SynSites', 'NonSynSites', 'APOBEC3A_Hairpin', 'MC.n_syn', 'MC.n_mis', 'MC.n_non', 'MC.n_spl', 'MC.n_ind', 'MC.wmis_cv', 'MC.wnon_cv', 'MC.wspl_cv', 'MC.wind_cv', 'MC.pmis_cv', 'MC.ptrunc_cv', 'MC.pallsubs_cv', 'MC.pind_cv', 'MC.qmis_cv', 'MC.qtrunc_cv', 'MC.qallsubs_cv', 'MC.pglobal_cv', 'MC.qglobal_cv'])
    df_data = df_data.rename(columns={'index': 'id'})

    #INPUT YOUR OWN CONNECTION STRING HERE
    conn_string = 'postgres://referral_writer:ProbioWriter@127.0.0.1:5432/curation'

    csv_path = dir_path+'/'+csv_file_name+'.csv'

    #perform COPY test and print result
    sql = "COPY {} FROM '{}' DELIMITER ',' CSV;".format(table_name, csv_path)

    table_create_sql = "CREATE TABLE IF NOT EXISTS {} (id SERIAL PRIMARY KEY, Gene TEXT, aaPos TEXT, nMut TEXT, ProtMut TEXT, Prot2Mut TEXT, DNAMut TEXT, CanMut TEXT, ConseqMut TEXT, Transcript TEXT, dN_dS TEXT, Community_Notes TEXT)".format(table_name)

    pg_conn = psycopg2.connect(conn_string)
    cur = pg_conn.cursor()
    cur.execute('DROP TABLE IF EXISTS {}'.format(table_name))
    cur.execute(table_create_sql)

    start_time = time.time()
    df_data.to_csv(csv_file_name+'.csv', index=False, header=False)
    cur.execute(sql)
    pg_conn.commit()
    cur.close()
    print("COPY {} duration: {} seconds".format(csv_file_name, time.time() - start_time))

def main():
    dir_path = os.getcwd()
    hotspot_file_path = dir_path +'/'+'hotspots.txt'
    warmspot_file_path = dir_path +'/'+'warmspots.txt'
    build_dataframe_to_table(dir_path, hotspot_file_path, 'hotspot_table', 'hotspot')
    build_dataframe_to_table(dir_path, warmspot_file_path, 'warmspots_table', 'warmspots')
    
if __name__ == "__main__":
    main()
    