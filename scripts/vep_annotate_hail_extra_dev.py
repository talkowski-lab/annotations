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
parser.add_argument('--mpc', dest='mpc_ht_uri', help='MPC scores HT')
parser.add_argument('--clinvar', dest='clinvar_vcf_uri', help='ClinVar VCF')
parser.add_argument('--inheritance', dest='inheritance_uri', help='File with inheritance codes (expected columns: approvedGeneSymbol, inheritance_code, genCC_classification)')
parser.add_argument('--revel', dest='revel_file', help='REVEL file')
parser.add_argument('--loeuf-v2', dest='loeuf_v2_uri', help='LOEUF scores from gnomAD v2.1.1')
parser.add_argument('--loeuf-v4', dest='loeuf_v4_uri', help='LOEUF scores from gnomAD v4.1')
parser.add_argument('--spliceAI-snv', dest='spliceAI_snv_uri', help='SpliceAI scores SNV HT')
parser.add_argument('--spliceAI-indel', dest='spliceAI_indel_uri', help='SpliceAI scores Indel HT')
parser.add_argument('--genes', dest='gene_list_tsv', help='OPTIONAL: Gene list file, tab-separated "gene_list_name"\t"gene_list_uri"')
parser.add_argument('--project-id', dest='project_id', help='Google Project ID')

args = parser.parse_args()

vcf_file = args.vcf_file
vep_annotated_vcf_name = args.vep_annotated_vcf_name
cores = args.cores  # string
mem = int(np.floor(float(args.mem)))
build = args.build
mpc_ht_uri = args.mpc_ht_uri
clinvar_vcf_uri = args.clinvar_vcf_uri
inheritance_uri = args.inheritance_uri
revel_file = args.revel_file
loeuf_v2_uri = args.loeuf_v2_uri
loeuf_v4_uri = args.loeuf_v4_uri
spliceAI_snv_uri = args.spliceAI_snv_uri
spliceAI_indel_uri = args.spliceAI_indel_uri
gene_list_tsv = args.gene_list_tsv
gcp_project = args.project_id

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

# annotate MPC
mpc = hl.read_table(mpc_ht_uri).key_by('locus','alleles')
mt = mt.annotate_rows(info = mt.info.annotate(MPC=mpc[mt.locus, mt.alleles].mpc))
        
# annotate ClinVar
if build=='GRCh38':
    clinvar_vcf = hl.import_vcf(clinvar_vcf_uri,
                            reference_genome='GRCh38',
                            force_bgz=clinvar_vcf_uri.split('.')[-1] in ['gz', 'bgz'],
                            skip_invalid_loci=True)
    # Grab ClinVar header
    clinvar_header = hl.get_vcf_metadata(clinvar_vcf_uri)
    mt = mt.annotate_rows(info = mt.info.annotate(CLNSIG=clinvar_vcf.rows()[mt.row_key].info.CLNSIG,
                                                  CLNREVSTAT=clinvar_vcf.rows()[mt.row_key].info.CLNREVSTAT,
                                                  CLNSIGCONF=clinvar_vcf.rows()[mt.row_key].info.CLNSIGCONF)
                                                  )
    # TODO: add ClinVar fields to annotate as input to workflow instead of being hardcoded?
    for clinvar_field in ['CLNSIG', 'CLNREVSTAT', 'CLNSIGCONF']:
        header['info'][clinvar_field] = clinvar_header['info'][clinvar_field]

# annotate REVEL
revel_ht = hl.import_table(revel_file, force_bgz=True)
revel_ht = revel_ht.annotate(chr='chr'+revel_ht['#chr']) 
build_chr = 'chr' if build=='GRCh38' else '#chr'
build_pos = 'grch38_pos' if build=='GRCh38' else 'hg19_pos'
revel_ht = revel_ht.annotate(locus=hl.locus(revel_ht[build_chr], hl.int(revel_ht[build_pos]), build),
                 alleles=hl.array([revel_ht.ref, revel_ht.alt]))
revel_ht = revel_ht.key_by('locus', 'alleles')
mt = mt.annotate_rows(info=mt.info.annotate(REVEL=revel_ht[mt.row_key].REVEL))

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

# explode rows and key by transcript
mt_by_transcript = mt.explode_rows(mt.vep.transcript_consequences)
mt_by_transcript = mt_by_transcript.key_rows_by(mt_by_transcript.vep.transcript_consequences.Feature)

# annotate LOEUF from gnomAD
loeuf_v2_ht = hl.read_table(loeuf_v2_uri).key_by('transcript')
loeuf_v4_ht = hl.read_table(loeuf_v4_uri).key_by('transcript')
mt_by_transcript = mt_by_transcript.annotate_rows(vep=mt_by_transcript.vep.annotate(
    transcript_consequences=mt_by_transcript.vep.transcript_consequences.annotate(
        LOEUF_v2=hl.if_else(hl.is_defined(loeuf_v2_ht[mt_by_transcript.row_key]), loeuf_v2_ht[mt_by_transcript.row_key]['oe_lof_upper'], ''),
        LOEUF_v2_decile=hl.if_else(hl.is_defined(loeuf_v2_ht[mt_by_transcript.row_key]), loeuf_v2_ht[mt_by_transcript.row_key]['oe_lof_upper_bin'], ''),
        LOEUF_v4=hl.if_else(hl.is_defined(loeuf_v4_ht[mt_by_transcript.row_key]), loeuf_v4_ht[mt_by_transcript.row_key]['lof.oe_ci.upper'], ''),
        LOEUF_v4_decile=hl.if_else(hl.is_defined(loeuf_v4_ht[mt_by_transcript.row_key]), loeuf_v4_ht[mt_by_transcript.row_key]['lof.oe_ci.upper_bin_decile'], '')
        )
    )
)

# annotate inheritance_code
inheritance_ht = hl.import_table(inheritance_uri).key_by('approvedGeneSymbol')
mt_by_gene = mt_by_transcript.key_rows_by(mt_by_transcript.vep.transcript_consequences.SYMBOL)
mt_by_gene = mt_by_gene.annotate_rows(vep=mt_by_gene.vep.annotate(
    transcript_consequences=mt_by_gene.vep.transcript_consequences.annotate(
        inheritance_code=hl.if_else(hl.is_defined(inheritance_ht[mt_by_gene.row_key]), inheritance_ht[mt_by_gene.row_key].inheritance_code, ''),
        genCC_classification=hl.if_else(hl.is_defined(inheritance_ht[mt_by_gene.row_key]), inheritance_ht[mt_by_gene.row_key].genCC_classification, '')
        )
    )
)

# OPTIONAL: annotate with gene list(s)
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
header['info']['REVEL'] = {'Description': 'REVEL scores.', 'Number': '.', 'Type': 'String'}

hl.export_vcf(dataset=mt, output=vep_annotated_vcf_name, metadata=header, tabix=True)