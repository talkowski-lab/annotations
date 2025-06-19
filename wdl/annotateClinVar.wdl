version 1.0

# CHANGE LOG: 
# 12/23/2024: added multiple gene list annotations (gene_list_tsv_tsv input) and overall spliceAI_score
# 3/18/2025: comment out SpliceAI score parsing below after updating input SpliceAI HT

import "scatterVCF.wdl" as scatterVCF
import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs
import "helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow AnnotateClinVar {

    input {
        Array[File] vep_vcf_files
        
        File clinvar_uri
        Array[String] clinvar_fields = ['CLNSIG', 'CLNREVSTAT', 'CLNSIGCONF', 'GENEINFO']

        String cohort_prefix
        String hail_docker
        
        String annotate_clinvar_script = "https://raw.githubusercontent.com/talkowski-lab/annotations/refs/heads/eren_dev/scripts/annotate_clinvar_hail.py"
        
        String genome_build='GRCh38'

        String? prefix  # optional prefix, for if filenames get too long

        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_annotate_add_genotypes
    }

    scatter (vcf_shard in vep_vcf_files) {
        if (!defined(prefix)) {
            String filename = basename(vcf_shard)
            String prefix_ = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_shard, ".vcf.gz") else basename(vcf_shard, ".vcf.bgz")
        }

        call annotateClinVar {
            input:
                vcf_file=vcf_shard,
                annotate_clinvar_script=annotate_clinvar_script,
                clinvar_uri=clinvar_uri,
                clinvar_fields=clinvar_fields,
                hail_docker=hail_docker,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_annotate,
                prefix=select_first([prefix_, prefix])
        }

        call helpers.addGenotypes as addGenotypes {
            input:
            annot_vcf_file=annotateClinVar.annot_vcf_file,
            vcf_file=vcf_shard,
            genome_build=genome_build,
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_annotate_add_genotypes
        }
    }

    output {
        Array[File] annot_vcf_files = addGenotypes.combined_vcf_file
        Array[File] annot_vcf_idx = addGenotypes.combined_vcf_idx
    }
}   

task annotateClinVar {
    input {
        File vcf_file
        File clinvar_uri
        Array[String] clinvar_fields

        String hail_docker
        String genome_build
        String annotate_clinvar_script

        String prefix
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String vep_annotated_vcf_name = "~{prefix}.annot.ClinVar.vcf.bgz"

    command <<<
        curl ~{annotate_clinvar_script} > annotate.py
        python3 annotate.py -i ~{vcf_file} -o ~{vep_annotated_vcf_name} --cores ~{cpu_cores} --mem ~{memory} \
        --build ~{genome_build} --clinvar ~{clinvar_uri} --clinvar-fields ~{sep=',' clinvar_fields} 
        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File annot_vcf_file = vep_annotated_vcf_name
        File annot_vcf_idx = vep_annotated_vcf_name + '.tbi'
        File hail_log = "hail_log.txt"
    }
}