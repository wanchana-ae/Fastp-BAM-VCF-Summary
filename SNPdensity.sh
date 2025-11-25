#!/bin/bash
set -euo pipefail
ml load vcftools
VCF_FILE="$1"

vcftools --gzvcf ${VCF_FILE} \
--out snp_index_SNPdensity \
--SNPdensity 1000000

cat snp_index_SNPdensity.snpden| \
grep -v "CHROM"|awk '{print "Chr"$1"\t"$2"\t"$2+1000"\t"".""\t"$3}' > snp_index_SNPdensity.bed