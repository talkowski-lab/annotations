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

workflow AnnotateEigenHT {
    input {
        String ht_uri
        String bucket_id
        String eigen_uri

        String genome_build='GRCh38'
        String hail_docker        
    }

    call helpers.getHailMTSize as getInputHTSize {
        input:
            mt_uri=ht_uri,
            hail_docker=hail_docker
    }
    
    call annotateEigenHT {
        input:
        ht_uri=ht_uri,
        eigen_uri=eigen_uri,
        bucket_id=bucket_id,
        genome_build=genome_build,
        hail_docker=hail_docker,
        ht_size=getInputHTSize.mt_size
    }

    output {
        String output_ht = annotateEigenHT.output_ht
    }
}

task annotateEigenHT {
    input {
        String ht_uri
        String eigen_uri
        String bucket_id
        String genome_build
        String hail_docker

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
    parser.add_argument('--eigen-uri', dest='eigen_uri', help='Noncoding BED file')
    parser.add_argument('--build', dest='build', help='Genome build')

    args = parser.parse_args()

    ht_uri = args.ht_uri
    bucket_id = args.bucket_id
    build = args.build
    eigen_uri = args.eigen_uri

    hl.init(min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                "spark.driver.memory": f"{int(np.floor(mem*0.8))}g"
            }, 
            tmp_dir="tmp", local_tmpdir="tmp", default_reference=build)

    eigen_ht = hl.read_table(eigen_uri)
    ht = hl.read_table(ht_uri)
    
    eigen_fields = ['Eigen-raw', 'Eigen-phred', 'Eigen-PC-raw', 'Eigen-PC-phred']
    ht = ht.annotate(**{eigen_field: hl.float(eigen_ht[ht.key][eigen_field]) for eigen_field in eigen_fields})

    prefix = os.path.basename(ht_uri).split('.ht')[0]
    filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{prefix}.Eigen.ht"
    pd.Series([filename]).to_csv('ht_uri.txt', index=False, header=None)
    ht.write(filename)    
    EOF
    python3 annotate_noncoding.py -i ~{ht_uri} --bucket-id ~{bucket_id} \
        --eigen-uri ~{eigen_uri} --build ~{genome_build}
    >>>

    output {
        String output_ht = read_lines('ht_uri.txt')[0]
    }
}