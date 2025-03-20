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

        String hail_docker
        String genome_build='GRCh38'
    }

    call annotateCADD {
        input:
        ht_uri=ht_uri,
        output_uri=output_uri,
        hail_docker=hail_docker,
        genome_build=genome_build
    }

    output {
        String output_ht = annotateCADD.output_ht
    }
}

task annotateCADD {
    input {
        String ht_uri
        String output_uri

        String hail_docker
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
    import hail as hl
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('-i', dest='ht_uri', help='Input HT')
    parser.add_argument('-o', dest='output_uri', help='Output URI')
    parser.add_argument('--cores', dest='cores', help='CPU cores')
    parser.add_argument('--mem', dest='mem', help='Memory')
    parser.add_argument('--build', dest='build', help='Genome build')

    args = parser.parse_args()

    ht_uri = args.ht_uri
    output_uri = args.output_uri
    cores = args.cores  # string
    mem = int(np.floor(float(args.mem)))
    build = args.build

    hl.init(default_reference=build,
            min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.speculation": 'true'
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
    )
    cadd_version = '1.6' if build=='GRCh38' else '1.4'
    cadd_ht = hl.experimental.load_dataset(name='CADD', version=cadd_version, reference_genome=build,
                                region='us-central1', cloud='gcp')
    # Annotate CADD
    ht = hl.read_table(ht_uri)
    ht = ht.annotate(CADD_raw_score=cadd_ht[ht.locus, ht.alleles].raw_score,
                              CADD_PHRED_score=cadd_ht[ht.locus, ht.alleles].PHRED_score)
    ht.write(output_uri, overwrite=True)

    EOF

    python3 annotateCADD.py -i ~{ht_uri} -o ~{output_uri} --cores ~{cpu_cores} --mem ~{memory} --build ~{genome_build}
    >>>

    output {
        String output_ht = output_uri
    }
}