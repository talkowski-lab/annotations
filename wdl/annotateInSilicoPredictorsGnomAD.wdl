version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow AnnotateInSilicoPredictorsGnomAD {

    input {
        File vcf_file

        String cadd_ht
        String pangolin_ht
        String phylop_ht
        String revel_ht
        String spliceai_ht

        String docker
        String annotate_in_silico_predictors_script
        String genome_build='GRCh38'
    }

    call annotateInSilicoPredictorsGnomAD {
        input:
        vcf_file=vcf_file,
        cadd_ht=cadd_ht,
        pangolin_ht=pangolin_ht,
        phylop_ht=phylop_ht,
        revel_ht=revel_ht,
        spliceai_ht=spliceai_ht,
        docker=docker,
        annotate_in_silico_predictors_script=annotate_in_silico_predictors_script,
        genome_build=genome_build
    }

    output {
        File output_vcf_file = annotateInSilicoPredictorsGnomAD.output_vcf_file
        File output_vcf_idx = annotateInSilicoPredictorsGnomAD.output_vcf_idx
    }
}   

task annotateInSilicoPredictorsGnomAD {
    input {
        File vcf_file

        String cadd_ht
        String pangolin_ht
        String phylop_ht
        String revel_ht
        String spliceai_ht

        String docker
        String annotate_in_silico_predictors_script
        String genome_build='GRCh38'
        
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
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
        docker: docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String annotated_vcf_name = "~{basename(vcf_file, file_ext)}.annot.in_silico_predictors.vcf.bgz"

    command <<<
        set -eou pipefail
        curl ~{annotate_in_silico_predictors_script} > annotate.py
        python3 annotate.py \
            --build ~{genome_build} \
            --cadd_ht ~{cadd_ht} \
            --pangolin_ht ~{pangolin_ht} \
            --phylop_ht ~{phylop_ht} \
            --revel_ht ~{revel_ht} \
            --spliceai_ht ~{spliceai_ht} \
            --vcf ~{vcf_file} \
            --output_vcf ~{annotated_vcf_name}
        tabix ~{annotated_vcf_name}
    >>>

    output {
        File output_vcf_file = annotated_vcf_name
        File output_vcf_idx = annotated_vcf_name + '.tbi'
    }
}