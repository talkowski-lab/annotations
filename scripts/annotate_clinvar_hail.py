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
parser.add_argument('--clinvar', dest='clinvar_vcf_uri', help='ClinVar VCF')
parser.add_argument('--clinvar-fields', dest='clinvar_fields', help='Fields to annotate from INFO in the ClinVar VCF')

args = parser.parse_args()

vcf_file = args.vcf_file
vep_annotated_vcf_name = args.vep_annotated_vcf_name
cores = args.cores  # string
mem = int(np.floor(float(args.mem)))
build = args.build
clinvar_vcf_uri = args.clinvar_vcf_uri
clinvar_fields = args.clinvar_fields.split(',')

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

# annotate ClinVar
if build=='GRCh38':
    clinvar_vcf = hl.import_vcf(clinvar_vcf_uri,
                            reference_genome='GRCh38',
                            force_bgz=clinvar_vcf_uri.split('.')[-1] in ['gz', 'bgz'],
                            skip_invalid_loci=True)
    # Grab ClinVar header
    clinvar_header = hl.get_vcf_metadata(clinvar_vcf_uri)
    mt = mt.annotate_rows(info = mt.info.annotate(**{clinvar_field: 
                                                     clinvar_vcf.rows()[mt.row_key].info[clinvar_field]
                                                     for clinvar_field in clinvar_fields}))

    for clinvar_field in clinvar_fields:
        header['info'][clinvar_field] = clinvar_header['info'][clinvar_field]

hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header, tabix=True)