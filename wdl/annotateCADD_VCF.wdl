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
        File vcf_file
        String BILLING_PROJECT_ID

        String hail_docker
        String cadd_ht_uri
        String genome_build='GRCh38'
    }

    call annotateCADD {
        input:
        vcf_file=vcf_file,
        cadd_ht_uri=cadd_ht_uri,
        hail_docker=hail_docker,
        genome_build=genome_build,
        BILLING_PROJECT_ID=BILLING_PROJECT_ID
    }

    output {
        File annotated_CADD_vcf_file = annotateCADD.output_vcf_file
        File annotated_CADD_vcf_idx = annotateCADD.output_vcf_idx
    }
}

task annotateCADD {
    input {
        File vcf_file
        String BILLING_PROJECT_ID

        String hail_docker
        String cadd_ht_uri
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
    
    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_uri = basename(vcf_file, file_ext) + '.CADD.vcf.bgz'
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
    parser.add_argument('-i', dest='vcf_file', help='Input VCF')
    parser.add_argument('-o', dest='output_uri', help='Output VCF URI')
    parser.add_argument('-c', dest='cadd_ht_uri', help='URI for CADD HT')    
    parser.add_argument('--cores', dest='cores', help='CPU cores')
    parser.add_argument('--mem', dest='mem', help='Memory')
    parser.add_argument('--build', dest='build', help='Genome build')
    parser.add_argument('--BILLING_PROJECT_ID', dest='BILLING_PROJECT_ID', help='BILLING_PROJECT_ID')

    args = parser.parse_args()

    vcf_file = args.vcf_file
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

    #split-multi
    def split_multi_ssc(mt):
        mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
        # only split variants that aren't already split
        bi = mt.filter_rows(hl.len(mt.alleles) == 2)
        bi = bi.annotate_rows(a_index=1, was_split=False, old_locus=bi.locus, old_alleles=bi.alleles)
        multi = mt.filter_rows(hl.len(mt.alleles) > 2)
        # Now split
        split = hl.split_multi(multi, permit_shuffle=True)
        sm = split.union_rows(bi)
        # sm = hl.split_multi(mt, permit_shuffle=True)
        if 'PL' in list(mt.entry.keys()):
            pl = hl.or_missing(hl.is_defined(sm.PL),
                            (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
            .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
            .map(lambda j: sm.PL[j])))))
            sm = sm.annotate_entries(PL = pl)
        split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index),
                                    AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]])
                                    ) 
            #GQ = hl.cond(hl.is_defined(pl[0]) & hl.is_defined(pl[1]) & hl.is_defined(pl[2]), hl.gq_from_pl(pl), sm.GQ) )
        mt = split_ds.drop('old_locus', 'old_alleles')
        return mt

    # Load VCF
    header = hl.get_vcf_metadata(vcf_file) 
    mt = hl.import_vcf(vcf_file, drop_samples=True, force_bgz=True, array_elements_required=False, call_fields=[])
    # Split multiallelics
    mt = split_multi_ssc(mt)

    # Annotate CADD
    mt = mt.annotate_rows(info=mt.info.annotate(
        CADD_raw_score=cadd_ht[mt.locus, mt.alleles].raw_score,
        CADD_PHRED_score=cadd_ht[mt.locus, mt.alleles].PHRED_score
        )
    )

    # Add to header
    header['info']['CADD_raw_score'] = {'Description': 'raw_score field from CADD HT.', 'Number': '.', 'Type': 'Float'}
    header['info']['CADD_PHRED_score'] = {'Description': 'PHRED_score field from CADD HT.', 'Number': '.', 'Type': 'Float'}
    hl.export_vcf(dataset=mt, output=output_uri, metadata=header, tabix=True)

    EOF

    python3 annotateCADD.py -i ~{vcf_file} -o ~{output_uri} -c ~{cadd_ht_uri} --cores ~{cpu_cores} --mem ~{memory} \
        --build ~{genome_build} --BILLING_PROJECT_ID ~{BILLING_PROJECT_ID}
    >>>

    output {
        File output_vcf_file = output_uri
        File output_vcf_idx = output_uri + '.tbi'
    }
}