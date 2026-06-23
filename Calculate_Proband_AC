version 1.0

## ProbandCohortAC.wdl
##
## Recalculates AC/AN/AF using only proband samples from a PED file,
## and annotates these values into the INFO field of the full VCF as
## PROBAND_AC, PROBAND_AN, PROBAND_AF.
##
## Inputs:
##   vcf          - Input VCF (indexed, bgzipped)
##   vcf_index    - Corresponding .tbi index
##   ped_file     - PED file (cols: family_id, sample_id, paternal_id, maternal_id, sex, phenotype)
##                  phenotype: 2 = proband, 1 = unaffected, 0 = unknown
##
## Outputs:
##   annotated_vcf             - VCF with PROBAND_AC, PROBAND_AN, PROBAND_AF in INFO
##   annotated_vcf_index       - Corresponding .tbi index
##   proband_list              - All proband sample IDs from PED
##   proband_list_in_vcf       - Proband sample IDs that were actually found in VCF header

workflow ProbandCohortAC {

    input {
        File   vcf
        File   vcf_index
        File   ped_file

        String docker      = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
        Int    mem_gb      = 16
        Int    disk_gb     = 200
        Int    cpu         = 4
        Int    preemptible = 3
    }

    call ExtractProbands {
        input:
            ped_file    = ped_file,
            docker      = docker,
            mem_gb      = mem_gb,
            disk_gb     = disk_gb,
            cpu         = cpu,
            preemptible = preemptible
    }

    call CalculateAndAnnotateProbandAC {
        input:
            vcf          = vcf,
            vcf_index    = vcf_index,
            proband_list = ExtractProbands.proband_list,
            docker       = docker,
            mem_gb       = mem_gb,
            disk_gb      = disk_gb,
            cpu          = cpu,
            preemptible  = preemptible
    }

    output {
        File annotated_vcf       = CalculateAndAnnotateProbandAC.annotated_vcf
        File annotated_vcf_index = CalculateAndAnnotateProbandAC.annotated_vcf_index
        File proband_list        = ExtractProbands.proband_list
        File proband_list_in_vcf = CalculateAndAnnotateProbandAC.proband_list_in_vcf
    }
}

# ---------------------------------------------------------------------------
# Task 1: Parse PED file to extract proband sample IDs (phenotype == 2)
#         Strips Windows-style carriage returns (\r) before parsing
# ---------------------------------------------------------------------------
task ExtractProbands {

    input {
        File   ped_file
        String docker
        Int    mem_gb
        Int    disk_gb
        Int    cpu
        Int    preemptible
    }

    command <<<
        set -euo pipefail

        # Strip carriage returns (\r) to handle Windows-formatted PED files,
        # then skip header (NR>1) and extract samples where phenotype == 2
        tr -d '\r' < ~{ped_file} \
            | awk 'NR>1 && $6 == 2 {print $2}' \
            | sort -u > proband_list.txt

        echo "Found $(wc -l < proband_list.txt) probands in PED"
        head -5 proband_list.txt
    >>>

    output {
        File proband_list = "proband_list.txt"
    }

    runtime {
        docker:      docker
        memory:      "~{mem_gb} GB"
        disks:       "local-disk ~{disk_gb} HDD"
        cpu:         cpu
        preemptible: preemptible
    }
}

# ---------------------------------------------------------------------------
# Task 2: Calculate proband-only AC/AN/AF and annotate into full VCF
# ---------------------------------------------------------------------------
task CalculateAndAnnotateProbandAC {

    input {
        File   vcf
        File   vcf_index
        File   proband_list
        String docker
        Int    mem_gb
        Int    disk_gb
        Int    cpu
        Int    preemptible
    }

    String vcf_basename = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail

        which bcftools || apt-get install -y bcftools

        # Step 1: Filter proband list to only samples present in VCF header
        bcftools query -l ~{vcf} > vcf_samples.txt

        comm -12 \
            <(sort ~{proband_list}) \
            <(sort vcf_samples.txt) \
            > proband_list_in_vcf.txt

        echo "Probands in PED:      $(wc -l < ~{proband_list})"
        echo "Probands in VCF:      $(wc -l < proband_list_in_vcf.txt)"
        echo "Probands not in VCF (skipped):"
        comm -23 \
            <(sort ~{proband_list}) \
            <(sort vcf_samples.txt)

        if [[ $(wc -l < proband_list_in_vcf.txt) -eq 0 ]]; then
            echo "ERROR: No proband samples found in VCF header. Check sample name formatting."
            exit 1
        fi

        # Step 2: Subset to probands present in VCF and compute AC/AN/AF,
        # then extract CHROM/POS/ID/REF/ALT + computed tags for annotation file
        bcftools view \
            --samples-file proband_list_in_vcf.txt \
            ~{vcf} | \
        bcftools +fill-tags -- -t AC,AN,AF | \
        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' \
            > proband_annot_raw.txt

        echo "Computed AC/AN/AF for $(wc -l < proband_annot_raw.txt) variants"

        # Step 3: Sort, bgzip, and index the annotation file
        sort -k1,1V -k2,2n proband_annot_raw.txt \
            | bgzip -c > proband_annot.txt.gz
        tabix -s1 -b2 -e2 proband_annot.txt.gz

        # Step 4: Create header lines for the new INFO fields
        cat > proband_ac_header.txt << 'EOF'
##INFO=<ID=PROBAND_AC,Number=A,Type=Integer,Description="Allele count in proband samples only (phenotype=2 in PED)">
##INFO=<ID=PROBAND_AN,Number=1,Type=Integer,Description="Total allele number in proband samples only (phenotype=2 in PED)">
##INFO=<ID=PROBAND_AF,Number=A,Type=Float,Description="Allele frequency in proband samples only (phenotype=2 in PED)">
EOF

        # Step 5: Annotate the full VCF with proband AC/AN/AF
        bcftools annotate \
            -a proband_annot.txt.gz \
            -h proband_ac_header.txt \
            -c CHROM,POS,ID,REF,ALT,PROBAND_AC,PROBAND_AN,PROBAND_AF \
            ~{vcf} \
            -Oz -o ~{vcf_basename}.proband_ac.vcf.gz

        bcftools index -t ~{vcf_basename}.proband_ac.vcf.gz

        # Sanity check
        echo "Spot check first 5 annotated variants:"
        bcftools query \
            -f '%ID\t%INFO/PROBAND_AC\t%INFO/PROBAND_AN\t%INFO/PROBAND_AF\n' \
            -r chr1 \
            ~{vcf_basename}.proband_ac.vcf.gz | head -5 || true

    >>>

    output {
        File annotated_vcf       = "~{vcf_basename}.proband_ac.vcf.gz"
        File annotated_vcf_index = "~{vcf_basename}.proband_ac.vcf.gz.tbi"
        File proband_list_in_vcf = "proband_list_in_vcf.txt"
    }

    runtime {
        docker:      docker
        memory:      "~{mem_gb} GB"
        disks:       "local-disk ~{disk_gb} HDD"
        cpu:         cpu
        preemptible: preemptible
    }
}
