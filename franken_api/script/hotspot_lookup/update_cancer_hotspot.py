#!/usr/bin/env python
# coding: utf-8

# ## Fetch the Cancer hotspot https://www.cancerhotspots.org/#/home
# 

import pandas as pd
import os
import re
import psycopg2
from sqlalchemy import create_engine


def build_hotspot(file_path):
    df = pd.read_csv(file_path, sep="\t", usecols=["Gene", "Residue", "Type", "Variants"])
    df = df.rename(columns={'Gene': 'gene', "Residue" : "residue", "Type" : "res_type", "Variants" : "variants"})
    df["varaint_arr"] = ""
    df["id"] = ""
    
    for index, row in df.iterrows():
        res_type = row["res_type"]
        residue = row["residue"]
        variants = row["variants"].split("|")
        hgsvp_str = ''
        for v in variants:
            if(res_type != "in-frame indel"):
                hgsvp = residue + v.split(":")[0]
            else:
                hgsvp = v.split(":")[0]

            hgsvp_str += hgsvp + ", "
        row["id"] = index + 1
        row["varaint_arr"] = hgsvp_str.rstrip(", ")
    
    df = df[["id", "gene", "residue", "res_type", "variants", "varaint_arr"]]
    
    return df

def main():
    
    dir_path = os.getcwd()

    file_path = dir_path+"/cancer_hotspots.txt"
    
    df_hotspot = build_hotspot(file_path)
    
    if(len(df_hotspot) > 0 ):
        db_engine = "postgresql+psycopg2://referral_writer:ProbioWriter@127.0.0.1:5432/curation"

        alchemyEngine           = create_engine(db_engine, pool_recycle=3600)
        postgreSQLConnection    = alchemyEngine.connect()
        postgreSQLTable         = "hotspot_summary"

        frame = df_hotspot.to_sql(postgreSQLTable, postgreSQLConnection, if_exists='replace',  index= False)
    else:
        print("hotspot empty")


# In[17]:


if __name__ == "__main__":
    main()

