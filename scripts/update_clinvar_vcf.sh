#!/bin/sh
<<<<<<< HEAD
'''
This script pulls the VCF from the latest GRCh38 ClinVar release at the time of running.
- Assumes the CHROM field does NOT have the "chr" prefix, even for GRCh38 build
- Assumes the second line of the VCF follows the format "##fileDate=%Y-%m-%d" for grabbing the release date
- Requires tabix, gsutil
'''
=======
# '''
# This script pulls the VCF from the latest GRCh38 ClinVar release at the time of running.
# - Assumes the CHROM field does NOT have the "chr" prefix, even for GRCh38 build
# - Assumes the second line of the VCF follows the format "##fileDate=%Y-%m-%d" for grabbing the release date
# - Requires tabix, gsutil
# '''
>>>>>>> eren_dev

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

# Get date of release
date_str=$(date -j -f "%Y-%m-%d" $(cat clinvar.vcf.gz | zcat | head -n 2 | tail -n 1 | cut -d '=' -f 2) "+%Y%m%d")

<<<<<<< HEAD
# Add chr
echo "Adding 'chr' prefix..."
cat clinvar.vcf.gz | zcat | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | bgzip > clinvar_${date_str}_chr_rename_GRCh38.vcf.gz
=======
# Add chr, replace MT with M
echo "Adding 'chr' prefix..."
cat clinvar.vcf.gz | zcat | sed 's/MT/M/g' | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | bgzip > clinvar_${date_str}_chr_rename_GRCh38.vcf.gz
>>>>>>> eren_dev

# Index
echo "Indexing using tabix..."
tabix clinvar_${date_str}_chr_rename_GRCh38.vcf.gz

# Copy to gsutil
gsutil -m cp clinvar_${date_str}_chr_rename_GRCh38.vcf.gz clinvar_${date_str}_chr_rename_GRCh38.vcf.gz.tbi gs://fc-0bc12741-801b-4c10-8d3c-92075b188d3c/resources/ClinVar/
echo "Copied (with index) to gs://fc-0bc12741-801b-4c10-8d3c-92075b188d3c/resources/ClinVar/clinvar_${date_str}_chr_rename_GRCh38.vcf.gz !"

# Remove intermediate files
rm clinvar.vcf.gz clinvar_${date_str}_chr_rename_GRCh38.vcf.gz clinvar_${date_str}_chr_rename_GRCh38.vcf.gz.tbi