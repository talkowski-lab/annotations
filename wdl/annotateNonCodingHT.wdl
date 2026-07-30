version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow annotateNonCoding {
    input {
        Array[String] ht_uris
        String bucket_id
        File noncoding_bed

        String genome_build='GRCh38'
        String hail_docker        
    }
    
    scatter (ht_uri in ht_uris) {
        call helpers.getHailMTSize as getInputHTSize {
            input:
                mt_uri=ht_uri,
                hail_docker=hail_docker
        }
        
        call annotateHTFromBed {
            input:
            ht_uri=ht_uri,
            bucket_id=bucket_id,
            genome_build=genome_build,
            noncoding_bed=noncoding_bed,
            hail_docker=hail_docker,
            ht_size=getInputHTSize.mt_size
        }
    }

    output {
        Array[String] output_ht = annotateHTFromBed.output_ht
    }
}

task annotateHTFromBed {
    input {
        String ht_uri
        String bucket_id
        String genome_build
        String hail_docker

        File noncoding_bed
        Float ht_size
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = ht_size
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
   
    command <<<
    cat <<EOF > annotate_noncoding.py
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
    parser.add_argument('--bucket-id', dest='bucket_id', help='Bucket ID')
    parser.add_argument('--cores', dest='cores', help='CPU cores')
    parser.add_argument('--mem', dest='mem', help='Memory')
    parser.add_argument('--noncoding', dest='noncoding_bed', help='Noncoding BED file')
    parser.add_argument('--build', dest='build', help='Genome build')

    args = parser.parse_args()

    ht_uri = args.ht_uri
    bucket_id = args.bucket_id
    cores = args.cores  # string
    mem = int(np.floor(float(args.mem)))
    build = args.build
    noncoding_bed = args.noncoding_bed

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    bed = hl.import_bed(
        noncoding_bed,
        reference_genome=build,
        skip_invalid_intervals=True
    )

    ht = hl.read_table(ht_uri)
    
    ht = ht.annotate(
        PREDICTED_NONCODING = bed.index(ht.locus, all_matches=True).target
    )

    prefix = os.path.basename(ht_uri).split('.ht')[0]
    filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{prefix}.noncoding.ht"
    pd.Series([filename]).to_csv('ht_uri.txt', index=False, header=None)
    ht.write(filename)    
    EOF
    python3 annotate_noncoding.py -i ~{ht_uri} --bucket-id ~{bucket_id} --cores ~{cpu_cores} --mem ~{memory} \
        --noncoding ~{noncoding_bed} --build ~{genome_build}
    >>>

    output {
        String output_ht = read_lines('ht_uri.txt')[0]
    }
}