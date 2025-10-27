#!/bin/bash
ml load bcftools
# File name VCF.gz
VCF_FILE="FF68_CN_170ea.vcf.gz"

# Check index have exit?
if [ ! -f "${VCF_FILE}.tbi" ]; then
    echo "create index VCF"
    bcftools index $VCF_FILE
fi

# Get name sample
samples=( $(bcftools query -l $VCF_FILE) )

# create output
OUTPUT_FILE="genotype_summary.tsv"
echo -e "Sample\tSNP_1/1\tINDEL_1/1" > $OUTPUT_FILE

# Count 1/1 SNP and INDEL each sample
for s in "${samples[@]}"; do
    echo "Counting on : $s"
    bcftools view -s $s $VCF_FILE | \
    bcftools query -f '%REF\t%ALT\t[%GT]\n' | \
    awk -v sample="$s" '
    {
        ref=$1; alt=$2; gt=$3
        type = (length(ref)==1 && length(alt)==1) ? "SNP" : "INDEL"
        if (gt=="1/1" || gt=="1|1") count[type]++
    }
    END { print sample "\t" (count["SNP"]+0) "\t" (count["INDEL"]+0) }
    ' >> $OUTPUT_FILE
    echo "Processed $s"

done

echo "All samples processed. Results in $OUTPUT_FILE"
