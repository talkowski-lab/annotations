version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow AnnotateGnomAD {
    input {
        String ht_uri
        String bucket_id
        String BILLING_PROJECT_ID

        String hail_docker
        String gnomADg_ht_uri
        String gnomADe_ht_uri
        String annotate_gnomad_script
        String genome_build='GRCh38'
    }

    call annotateGnomAD {
        input:
        ht_uri=ht_uri,
        bucket_id=bucket_id,
        gnomADg_ht_uri=gnomADg_ht_uri,
        gnomADe_ht_uri=gnomADe_ht_uri,
        annotate_gnomad_script=annotate_gnomad_script,
        hail_docker=hail_docker,
        genome_build=genome_build,
        BILLING_PROJECT_ID=BILLING_PROJECT_ID
    }

    output {
        String output_ht = annotateGnomAD.output_ht
    }
}

task annotateGnomAD {
    input {
        String ht_uri
        String bucket_id
        String BILLING_PROJECT_ID

        String hail_docker
        String gnomADg_ht_uri
        String gnomADe_ht_uri
        String annotate_gnomad_script
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
    curl ~{annotate_gnomad_script} > annotate_gnomAD.py

    python3 annotate_gnomAD.py -i ~{ht_uri} --bucket_id ~{bucket_id} -g ~{gnomADg_ht_uri} -e ~{gnomADe_ht_uri} \
        --cores ~{cpu_cores} --mem ~{memory} --build ~{genome_build} --BILLING_PROJECT_ID ~{BILLING_PROJECT_ID}
    >>>

    output {
        String output_ht = read_lines('ht_uri.txt')[0]
    }
}