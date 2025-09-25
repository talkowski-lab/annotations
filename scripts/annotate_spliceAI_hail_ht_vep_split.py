###
# Pulled from annotateSpliceAI task in vepAnnotateHailExtra_dev.wdl on 3/20/2025.

## CHANGE LOG:
'''
3/20/2025:
- limit SpliceAI annotations to only sites with splice variants
'''
###

from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os
import json
import argparse
import datetime

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('-i', dest='ht_uri', help='Input HT URI')
parser.add_argument('--bucket_id', dest='bucket_id', help='Bucket ID')
parser.add_argument('--cores', dest='cores', help='CPU cores')
parser.add_argument('--mem', dest='mem', help='Memory')
parser.add_argument('--build', dest='build', help='Genome build')
parser.add_argument('--spliceAI-uri', dest='spliceAI_uri', help='SpliceAI scores SNV/Indel HT')

args = parser.parse_args()

ht_uri = args.ht_uri
bucket_id = args.bucket_id
cores = args.cores  # string
mem = int(np.floor(float(args.mem)))
build = args.build
spliceAI_uri = args.spliceAI_uri

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

ht = hl.read_table(ht_uri)
spliceAI_ht = hl.read_table(spliceAI_uri)

spliceAI_fields = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'spliceAI_score']

ht = ht.annotate(**{field: spliceAI_ht[ht.key][field] for field in spliceAI_fields})

prefix = os.path.basename(ht_uri).split('.ht')[0]
filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{prefix}.vep.ht"
pd.Series([filename]).to_csv('ht_uri.txt', index=False, header=None)
ht.write(filename)