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
parser.add_argument('-i', dest='ht_uri', help='Input HT')
parser.add_argument('-g', dest='gnomADg_ht_uri', help='URI for gnomAD genomes HT')    
parser.add_argument('-e', dest='gnomADe_ht_uri', help='URI for gnomAD exomes HT')    
parser.add_argument('--bucket_id', dest='bucket_id', help='Bucket ID')
parser.add_argument('--cores', dest='cores', help='CPU cores')
parser.add_argument('--mem', dest='mem', help='Memory')
parser.add_argument('--build', dest='build', help='Genome build')
parser.add_argument('--BILLING_PROJECT_ID', dest='BILLING_PROJECT_ID', help='BILLING_PROJECT_ID')

args = parser.parse_args()

ht_uri = args.ht_uri
bucket_id = args.bucket_id
gnomADg_ht_uri = args.gnomADg_ht_uri
gnomADe_ht_uri = args.gnomADe_ht_uri
cores = args.cores  # string
mem = int(np.floor(float(args.mem)))
build = args.build
BILLING_PROJECT_ID = args.BILLING_PROJECT_ID

hl.init(default_reference=build,
        min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
        gcs_requester_pays_configuration=BILLING_PROJECT_ID
)

# Extract version strings
gnomade_version = os.path.basename(gnomADe_ht_uri).split('.exomes.')[1].split('.sites.')[0]
gnomadg_version = os.path.basename(gnomADg_ht_uri).split('.genomes.')[1].split('.sites.')[0]

# Start with small ht
ht = hl.read_table(ht_uri).repartition(200)

# First annotate gnomADe
gnomade_ht = hl.read_table(gnomADe_ht_uri).select('freq')
gnomade_annot = gnomade_ht.annotate(tmp=ht[gnomade_ht.key])
gnomade_annot = gnomade_annot.filter(hl.is_defined(gnomade_annot.tmp))
gnomade_annot = gnomade_annot.select(
    key=gnomade_annot.key,
    gnomADe_AC_tmp=gnomade_annot.freq.AC,
    gnomADe_AF_tmp=gnomade_annot.freq.AF,
    gnomADe_AN_tmp=gnomade_annot.freq.AN
)

# Join back with original ht
ht = ht.annotate(**gnomade_annot[ht.key])

# Now do the same for gnomADg
gnomadg_ht = hl.read_table(gnomADg_ht_uri).select('freq')
gnomadg_annot = gnomadg_ht.annotate(tmp=ht[gnomadg_ht.key])
gnomadg_annot = gnomadg_annot.filter(hl.is_defined(gnomadg_annot.tmp))
gnomadg_annot = gnomadg_annot.select(
    key=gnomadg_annot.key,
    gnomADg_AC_tmp=gnomadg_annot.freq.AC,
    gnomADg_AF_tmp=gnomadg_annot.freq.AF,
    gnomADg_AN_tmp=gnomadg_annot.freq.AN
)

# Final annotate to add both sets
ht = ht.annotate(**gnomadg_annot[ht.key])

# Collect frequency index maps
gnomade_freq_idx_map = pd.Series(gnomade_ht.freq_index_dict.collect()[0])
gnomadg_freq_idx_map = pd.Series(gnomadg_ht.freq_index_dict.collect()[0])
gnomadg_freq_idx_map.index = gnomadg_freq_idx_map.index.str.replace('-', '_')

# Define relevant population tags
populations = ['afr', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'oth', 'remaining', 'sas']

exome_populations = [
    pop for pop in populations
    if gnomade_freq_idx_map.index.str.contains(pop).sum() > 0
]

genome_populations = [
    pop for pop in populations
    if gnomadg_freq_idx_map.index.str.contains(pop).sum() > 0
]

# Map field names based on version
if gnomade_freq_idx_map.index.str.contains('gnomad').sum() != 0:
    gnomade_annot_field_map = {
        f"gnomADe_{gnomade_version}": ('gnomad', gnomade_freq_idx_map['gnomad']),
        **{
            f"gnomADe_{gnomade_version}_{pop.upper()}": (f"gnomad_{pop}", gnomade_freq_idx_map[f"gnomad_{pop}"])
            for pop in exome_populations
        }
    }
else:
    gnomade_annot_field_map = {
        f"gnomADe_{gnomade_version}": ('adj', gnomade_freq_idx_map['adj']),
        **{
            f"gnomADe_{gnomade_version}_{pop.upper()}": (f"{pop}_adj", gnomade_freq_idx_map[f"{pop}_adj"])
            for pop in exome_populations
        }
    }

gnomadg_annot_field_map = {
    f"gnomADg_{gnomadg_version}": ('adj', gnomadg_freq_idx_map['adj']),
    **{
        f"gnomADg_{gnomadg_version}_{pop.upper()}": (f"{pop}_adj", gnomadg_freq_idx_map[f"{pop}_adj"])
        for pop in genome_populations
    }
}

# Annotate AC, AF, AN from gnomADe and gnomADg
ht = ht.annotate(
    **{
        f"{new_field}_AC": ht.gnomADe_AC_tmp[hl.int32(idx)]
        for new_field, (old_field, idx) in gnomade_annot_field_map.items()
    },
    **{
        f"{new_field}_AF": ht.gnomADe_AF_tmp[hl.int32(idx)]
        for new_field, (old_field, idx) in gnomade_annot_field_map.items()
    },
    **{
        f"{new_field}_AN": ht.gnomADe_AN_tmp[hl.int32(idx)]
        for new_field, (old_field, idx) in gnomade_annot_field_map.items()
    },
    **{
        f"{new_field}_AC": ht.gnomADg_AC_tmp[hl.int32(idx)]
        for new_field, (old_field, idx) in gnomadg_annot_field_map.items()
    },
    **{
        f"{new_field}_AF": ht.gnomADg_AF_tmp[hl.int32(idx)]
        for new_field, (old_field, idx) in gnomadg_annot_field_map.items()
    },
    **{
        f"{new_field}_AN": ht.gnomADg_AN_tmp[hl.int32(idx)]
        for new_field, (old_field, idx) in gnomadg_annot_field_map.items()
    }
)

# Drop temporary fields
row_fields = pd.Series(list(ht.row))
tmp_row_fields = row_fields[row_fields.astype(str).str.contains('_tmp')].tolist()
ht = ht.drop(*tmp_row_fields)

prefix = os.path.basename(ht_uri).split('.ht')[0]
output_uri = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{prefix}.gnomAD_AC_AN.ht"
pd.Series([output_uri]).to_csv('ht_uri.txt', index=False, header=None)
ht.write(output_uri, overwrite=True)

