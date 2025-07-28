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

# annotate SpliceAI scores
ht_by_transcript = ht.explode(ht.vep.transcript_consequences)
ht_by_locus_and_gene = ht_by_transcript.key_by('locus', 'alleles', ht_by_transcript.vep.transcript_consequences.SYMBOL)

spliceAI_ht = hl.read_table(spliceAI_uri)

# NEW 3/20/2025: limit SpliceAI annotations to only sites with splice variants
splice_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant']

has_splice_var = (
    hl.set(splice_vars)
    .intersection(hl.set(ht_by_locus_and_gene.vep.transcript_consequences.Consequence))
    .size() > 0
)

# Leave out ALLELE/SYMBOL because redundant
fields = 'ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL'.split('|')[2:]

ht_by_locus_and_gene = ht_by_locus_and_gene.annotate(
    SpliceAI_raw=hl.or_missing(
        has_splice_var, 
        spliceAI_ht[ht_by_locus_and_gene.row_key].SpliceAI
    )
)

ht_by_locus_and_gene = ht_by_locus_and_gene.annotate(
    vep=ht_by_locus_and_gene.vep.annotate(
        transcript_consequences=(
            ht_by_locus_and_gene.vep.transcript_consequences.annotate(
                **{
                    field: hl.if_else(
                        hl.is_defined(ht_by_locus_and_gene.SpliceAI_raw), 
                        ht_by_locus_and_gene.SpliceAI_raw.split('=')[1].split('\|')[i+2], 
                        ''
                    )
                    for i, field in enumerate(fields)
                }
            )
        )
    )
)

# Overall SpliceAI score
score_fields = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']

ht_by_locus_and_gene = ht_by_locus_and_gene.annotate(
    vep=ht_by_locus_and_gene.vep.annotate(
        transcript_consequences=(
            ht_by_locus_and_gene.vep.transcript_consequences.annotate(
                spliceAI_score=hl.str(
                    hl.max([
                        hl.or_missing(
                            ht_by_locus_and_gene.vep.transcript_consequences[field] != '',
                            hl.float(ht_by_locus_and_gene.vep.transcript_consequences[field])
                        )
                        for field in score_fields
                    ])
                )
            )
        )
    )
)


ht_by_gene = ht_by_locus_and_gene
ht_by_gene = (ht_by_gene.group_by(ht_by_gene.locus, ht_by_gene.alleles)
    .aggregate(transcript_consequences = hl.agg.collect(ht_by_gene.vep.transcript_consequences)))

ht = ht.annotate(vep=hl.Struct(**{'transcript_consequences': ht_by_gene[ht.key].transcript_consequences}))

# only adds new CSQ fields to header, overwrites if already present
fields = list(ht.vep.transcript_consequences[0])

prefix = os.path.basename(ht_uri).split('.ht')[0]
filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{prefix}.vep.ht"
pd.Series([filename]).to_csv('ht_uri.txt', index=False, header=None)
ht.write(filename)