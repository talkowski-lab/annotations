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

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('-i', dest='vcf_file', help='Input VCF file')
parser.add_argument('-o', dest='vep_annotated_vcf_name', help='Output filename')
parser.add_argument('--cores', dest='cores', help='CPU cores')
parser.add_argument('--mem', dest='mem', help='Memory')
parser.add_argument('--build', dest='build', help='Genome build')
parser.add_argument('--spliceAI-uri', dest='spliceAI_uri', help='SpliceAI scores SNV/Indel HT')

args = parser.parse_args()

vcf_file = args.vcf_file
vep_annotated_vcf_name = args.vep_annotated_vcf_name
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

header = hl.get_vcf_metadata(vcf_file) 
mt = hl.import_vcf(vcf_file, drop_samples=True, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build)

# expect first field to be Allele, regardless of 'Format: ' (or anything before) being present
csq_columns = ('Allele' + header['info']['CSQ']['Description'].split('Allele')[-1]).split('|')

# split VEP CSQ string
mt = mt.annotate_rows(vep=mt.info)
transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                    {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                        for i, col in enumerate(csq_columns)}), 
                                                        hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                        for i, col in enumerate(csq_columns)})))

mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))

# annotate SpliceAI scores
mt_by_transcript = mt.explode_rows(mt.vep.transcript_consequences)
mt_by_locus_and_gene = mt_by_transcript.key_rows_by('locus', 'alleles', mt_by_transcript.vep.transcript_consequences.SYMBOL)

spliceAI_ht = hl.read_table(spliceAI_uri)

# NEW 3/20/2025: limit SpliceAI annotations to only sites with splice variants
splice_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant']

has_splice_var = (
    hl.set(splice_vars)
    .intersection(hl.set(mt_by_locus_and_gene.vep.transcript_consequences.Consequence))
    .size() > 0
)

# Leave out ALLELE/SYMBOL because redundant
fields = 'ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL'.split('|')[2:]

mt_by_locus_and_gene = mt_by_locus_and_gene.annotate_rows(
    SpliceAI_raw=hl.or_missing(
        has_splice_var, 
        spliceAI_ht[mt_by_locus_and_gene.row_key].SpliceAI
    )
)

mt_by_locus_and_gene = mt_by_locus_and_gene.annotate_rows(
    vep=mt_by_locus_and_gene.vep.annotate(
        transcript_consequences=(
            mt_by_locus_and_gene.vep.transcript_consequences.annotate(
                **{
                    field: hl.if_else(
                        hl.is_defined(mt_by_locus_and_gene.SpliceAI_raw), 
                        mt_by_locus_and_gene.SpliceAI_raw.split('=')[1].split('\|')[i+2], 
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

mt_by_locus_and_gene = mt_by_locus_and_gene.annotate_rows(
    vep=mt_by_locus_and_gene.vep.annotate(
        transcript_consequences=(
            mt_by_locus_and_gene.vep.transcript_consequences.annotate(
                spliceAI_score=hl.str(
                    hl.max([
                        hl.or_missing(
                            mt_by_locus_and_gene.vep.transcript_consequences[field] != '',
                            hl.float(mt_by_locus_and_gene.vep.transcript_consequences[field])
                        )
                        for field in score_fields
                    ])
                )
            )
        )
    )
)


mt_by_gene = mt_by_locus_and_gene
mt_by_gene = (mt_by_gene.group_rows_by(mt_by_gene.locus, mt_by_gene.alleles)
    .aggregate_rows(vep = hl.agg.collect(mt_by_gene.vep))).result()

fields = list(mt_by_gene.vep.transcript_consequences[0])
new_csq = mt_by_gene.vep.transcript_consequences.scan(lambda i, j: 
                                    hl.str('|').join(hl.array([i]))
                                    +','+hl.str('|').join(hl.array([j[col] if col!='Consequence' else 
                                                                hl.str('&').join(j[col]) 
                                                                for col in list(fields)])), '')[-1][1:]
mt_by_gene = mt_by_gene.annotate_rows(CSQ=new_csq)
mt = mt.annotate_rows(info=mt.info.annotate(CSQ=mt_by_gene.rows()[mt.row_key].CSQ))
mt = mt.drop('vep')

# only adds new CSQ fields to header, overwrites if already present
csq_fields_str = 'Format: ' + '|'.join(fields)

header['info']['CSQ'] = {'Description': csq_fields_str, 'Number': '.', 'Type': 'String'}
hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header, tabix=True)
