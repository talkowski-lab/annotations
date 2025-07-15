version 1.0

# CHANGE LOG: 
# 12/23/2024: added multiple gene list annotations (gene_list_tsv_tsv input) and overall spliceAI_score
# 3/18/2025: comment out SpliceAI score parsing below after updating input SpliceAI HT

import "scatterVCF.wdl" as scatterVCF
import "mergeSplitVCF.wdl" as mergeSplitVCF
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/mergeVCFs.wdl" as mergeVCFs
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow vepAnnotateHailExtra {

    input {
        Array[File] vep_vcf_files

        File revel_file
        File clinvar_vcf_uri
        File inheritance_uri
        String mpc_ht_uri
        String loeuf_v2_uri
        String loeuf_v4_uri

        String cohort_prefix
        String hail_docker
        
        String vep_annotate_hail_extra_python_script = "https://raw.githubusercontent.com/talkowski-lab/annotations/refs/heads/eren_dev/scripts/vep_annotate_hail_extra_dev.py"
        String split_vcf_hail_script = "https://raw.githubusercontent.com/talkowski-lab/annotations/refs/heads/main/scripts/split_vcf_hail.py"
        String annotate_spliceAI_script = "https://raw.githubusercontent.com/talkowski-lab/annotations/refs/heads/eren_dev/scripts/annotate_spliceAI_hail.py"
        
        String genome_build='GRCh38'

        String spliceAI_uri='NA'
        String noncoding_bed='NA'
        String gene_list_tsv='NA'

        RuntimeAttr? runtime_attr_annotate_noncoding      
        RuntimeAttr? runtime_attr_annotate_extra
        RuntimeAttr? runtime_attr_annotate_spliceAI
        RuntimeAttr? runtime_attr_annotate_add_genotypes
    }

    scatter (vcf_shard in vep_vcf_files) {
        if (noncoding_bed!='NA') {
            call annotateFromBed as annotateNonCoding {
                input:
                vcf_file=vcf_shard,
                noncoding_bed=select_first([noncoding_bed]),
                hail_docker=hail_docker,
                genome_build=genome_build,
                filter=false,
                runtime_attr_override=runtime_attr_annotate_noncoding
            }
        }

        call annotateExtra {
            input:
                vcf_file=select_first([annotateNonCoding.noncoding_vcf, vcf_shard]),
                vep_annotate_hail_extra_python_script=vep_annotate_hail_extra_python_script,
                loeuf_v2_uri=loeuf_v2_uri,
                loeuf_v4_uri=loeuf_v4_uri,
                revel_file=revel_file,
                revel_file_idx=revel_file+'.tbi',
                clinvar_vcf_uri=clinvar_vcf_uri,
                inheritance_uri=inheritance_uri,
                gene_list_tsv=select_first([gene_list_tsv, 'NA']),
                mpc_ht_uri=mpc_ht_uri,
                hail_docker=hail_docker,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_annotate_extra
        }

        if (spliceAI_uri!='NA') {
            call annotateSpliceAI {
                input:
                vcf_file=annotateExtra.annot_vcf_file,
                spliceAI_uri=spliceAI_uri,
                genome_build=genome_build,
                hail_docker=hail_docker,
                annotate_spliceAI_script=annotate_spliceAI_script,
                runtime_attr_override=runtime_attr_annotate_spliceAI
            }
        }
        File annot_vcf_file = select_first([annotateSpliceAI.annot_vcf_file, annotateExtra.annot_vcf_file])

        call helpers.addGenotypes as addGenotypes {
            input:
            annot_vcf_file=annot_vcf_file,
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

task annotateFromBed {
    input {
        File vcf_file
        String noncoding_bed 
        String hail_docker
        String genome_build
        Boolean filter
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, 'GB')
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

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_filename = basename(vcf_file, file_ext) + '_noncoding_annot.vcf.bgz' 
   
    command <<<
    cat <<EOF > annotate_noncoding.py
    from pyspark.sql import SparkSession
    import hail as hl
    import numpy as np
    import sys
    import ast
    import os

    vcf_file = sys.argv[1]
    noncoding_bed = sys.argv[2]
    cores = sys.argv[3]  # string
    mem = int(np.floor(float(sys.argv[4])))
    build = sys.argv[5]
    output_filename = sys.argv[6]
    filter = ast.literal_eval(sys.argv[7].capitalize())

    hl.init(min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.speculation": 'true'
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
                        )

    bed = hl.import_bed(noncoding_bed, reference_genome=build, skip_invalid_intervals=True)
    mt = hl.import_vcf(vcf_file, drop_samples=True, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build)
    mt = mt.annotate_rows(info=mt.info.annotate(PREDICTED_NONCODING=bed[mt.locus].target))

    # filter only annotated
    if filter:
        mt = mt.filter_rows(hl.is_defined(mt.info.PREDICTED_NONCODING))

    header = hl.get_vcf_metadata(vcf_file)
    header['info']['PREDICTED_NONCODING'] = {'Description': "Class(es) of noncoding elements disrupted by SNV/Indel.", 
                                            'Number': '.', 'Type': 'String'}
    hl.export_vcf(mt, output_filename, metadata=header, tabix=True)
    EOF
    python3 annotate_noncoding.py ~{vcf_file} ~{noncoding_bed} ~{cpu_cores} ~{memory} ~{genome_build} ~{output_filename} ~{filter}
    >>>

    output {
        File noncoding_vcf = output_filename
        File noncoding_vcf_idx = output_filename + '.tbi'
    }
}

task annotateExtra {
    input {
        File vcf_file
        String loeuf_v2_uri
        String loeuf_v4_uri
        File revel_file
        File revel_file_idx
        File clinvar_vcf_uri
        File inheritance_uri
        
        String gene_list_tsv
        String mpc_ht_uri

        String hail_docker
        String genome_build
        String vep_annotate_hail_extra_python_script
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

    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")
    String vep_annotated_vcf_name = "~{prefix}.annot.vcf.bgz"

    command <<<
        curl ~{vep_annotate_hail_extra_python_script} > annotate.py
        python3 annotate.py -i ~{vcf_file} -o ~{vep_annotated_vcf_name} --cores ~{cpu_cores} --mem ~{memory} \
        --build ~{genome_build} --loeuf-v2 ~{loeuf_v2_uri} --loeuf-v4 ~{loeuf_v4_uri} \
        --mpc ~{mpc_ht_uri} --clinvar ~{clinvar_vcf_uri} --inheritance ~{inheritance_uri} \
        --revel ~{revel_file} --genes ~{gene_list_tsv} 
        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File annot_vcf_file = vep_annotated_vcf_name
        File annot_vcf_idx = vep_annotated_vcf_name + '.tbi'
        File hail_log = "hail_log.txt"
    }
}

task annotateSpliceAI {
    input {
        File vcf_file
        String spliceAI_uri

        String hail_docker
        String genome_build
        String annotate_spliceAI_script
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

    String filename = basename(vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_file, ".vcf.gz") else basename(vcf_file, ".vcf.bgz")
    String vep_annotated_vcf_name = "~{prefix}.SpliceAI.annot.vcf.bgz"

    command <<<
    curl ~{annotate_spliceAI_script} > annotateSpliceAI.py
    python3 annotateSpliceAI.py -i ~{vcf_file} -o ~{vep_annotated_vcf_name} --cores ~{cpu_cores} --mem ~{memory} \
    --build ~{genome_build} --spliceAI-uri ~{spliceAI_uri} 
    cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File annot_vcf_file = vep_annotated_vcf_name
        File annot_vcf_idx = vep_annotated_vcf_name + '.tbi'
        File hail_log = "hail_log.txt"
    }
}