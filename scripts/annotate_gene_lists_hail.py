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
parser.add_argument('--inheritance', dest='inheritance_uri', help='File with inheritance codes (expected columns: approvedGeneSymbol, inheritance_code, genCC_classification)')
parser.add_argument('--genes', dest='gene_list_tsv', help='OPTIONAL: Gene list file, tab-separated "gene_list_name"\t"gene_list_uri"')

args = parser.parse_args()

vcf_file = args.vcf_file
vep_annotated_vcf_name = args.vep_annotated_vcf_name
cores = args.cores  # string
mem = int(np.floor(float(args.mem)))
build = args.build
inheritance_uri = args.inheritance_uri
gene_list_tsv = args.gene_list_tsv

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

# Explode rows and key by transcript
mt_by_transcript = mt.explode_rows(mt.vep.transcript_consequences)
mt_by_transcript = mt_by_transcript.key_rows_by(mt_by_transcript.vep.transcript_consequences.Feature)

# Annotate inheritance_code
inheritance_ht = hl.import_table(inheritance_uri).key_by('approvedGeneSymbol')
mt_by_gene = mt_by_transcript.key_rows_by(mt_by_transcript.vep.transcript_consequences.SYMBOL)
mt_by_gene = mt_by_gene.annotate_rows(vep=mt_by_gene.vep.annotate(
    transcript_consequences=mt_by_gene.vep.transcript_consequences.annotate(
        inheritance_ht_inheritance_code=hl.if_else(hl.is_defined(inheritance_ht[mt_by_gene.row_key]), inheritance_ht[mt_by_gene.row_key].inheritance_code, '')
        )
    )
)

# OPTIONAL: Annotate with gene list(s)
if gene_list_tsv!='NA':
    gene_list_uris = pd.read_csv(gene_list_tsv, sep='\t', header=None).set_index(0)[1].to_dict()
    gene_lists = {gene_list_name: pd.read_csv(uri, header=None)[0].tolist() 
                for gene_list_name, uri in gene_list_uris.items()}

    mt_by_gene = mt_by_gene.annotate_rows(vep=mt_by_gene.vep.annotate(
        transcript_consequences=mt_by_gene.vep.transcript_consequences.annotate(
            gene_list=hl.str('&').join(hl.array([hl.or_missing(hl.array(gene_list).contains(mt_by_gene.vep.transcript_consequences.SYMBOL), gene_list_name) 
                for gene_list_name, gene_list in gene_lists.items()]).filter(hl.is_defined)))))

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