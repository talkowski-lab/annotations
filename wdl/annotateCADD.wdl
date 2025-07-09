version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow AnnotateCADD {
    input {
        String ht_uri
        String output_uri
        String BILLING_PROJECT_ID

        String hail_docker
        String cadd_ht_uri
        String genome_build='GRCh38'
    }

    call annotateCADD {
        input:
        ht_uri=ht_uri,
        output_uri=output_uri,
        cadd_ht_uri=cadd_ht_uri,
        hail_docker=hail_docker,
        genome_build=genome_build,
        BILLING_PROJECT_ID=BILLING_PROJECT_ID
    }

    output {
        String output_ht = annotateCADD.output_ht
    }
}

task annotateCADD {
    input {
        String ht_uri
        String output_uri
        String BILLING_PROJECT_ID

        String hail_docker
        String cadd_ht_uri
        String genome_build='GRCh38'
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
        disk_gb: ceil(base_disk_gb),
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
    cat <<EOF > annotateCADD.py
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
    parser.add_argument('-i', dest='ht_uri', help='Input HT')
    parser.add_argument('-o', dest='output_uri', help='Output URI')
    parser.add_argument('-c', dest='cadd_ht_uri', help='URI for CADD HT')    
    parser.add_argument('--cores', dest='cores', help='CPU cores')
    parser.add_argument('--mem', dest='mem', help='Memory')
    parser.add_argument('--build', dest='build', help='Genome build')
    parser.add_argument('--BILLING_PROJECT_ID', dest='BILLING_PROJECT_ID', help='BILLING_PROJECT_ID')

    args = parser.parse_args()

    ht_uri = args.ht_uri
    output_uri = args.output_uri
    cadd_ht_uri = args.cadd_ht_uri
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
    cadd_ht = hl.read_table(cadd_ht_uri)

    # Annotate CADD
    ht = hl.read_table(ht_uri)
    ht = ht.annotate(CADD_raw_score=cadd_ht[ht.locus, ht.alleles].raw_score,
                    CADD_PHRED_score=cadd_ht[ht.locus, ht.alleles].PHRED_score)
    ht.write(output_uri, overwrite=True)

    EOF

    python3 annotateCADD.py -i ~{ht_uri} -o ~{output_uri} -c ~{cadd_ht_uri} --cores ~{cpu_cores} --mem ~{memory} \
        --build ~{genome_build} --BILLING_PROJECT_ID ~{BILLING_PROJECT_ID}
    >>>

    output {
        String output_ht = output_uri
    }
}