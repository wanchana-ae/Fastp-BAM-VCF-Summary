#!/bin/bash

# Load samtools module
ml load samtools

INFO_FILE="/fs4/Mining_SNP/Borassus_flabellifer_Linn/Germplasm/name_bam.txt"
OUT_count="/fs4/Mining_SNP/Borassus_flabellifer_Linn/Germplasm/Report/"
Number_CPU=8        # threads per sample
Parallel_Sample=4   # number of samples to run in parallel

mkdir -p ${OUT_count}Count_Mapped/

process_sample() {
    BAM_input=$1
    Number_CPU=$2
    OUT_count=$3

    INPUT_ID=${BAM_input##*/}
    SampleID=${INPUT_ID%.*}

    mkdir -p ${OUT_count}Count_Mapped/${SampleID}

    echo "===================================================================================="
    echo "Processing sample: ${SampleID}"
    echo "===================================================================================="

    # samtools flagstat
    echo "samtools flagstat : ${BAM_input}"
    samtools flagstat --threads ${Number_CPU} ${BAM_input} \
        > ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_flagstat.txt

    # Count mapped read
    echo "Count mapped read : ${BAM_input}"
    samtools view -c -F 0x04 -@ ${Number_CPU} ${BAM_input} \
        | awk -v a="${SampleID}" '{print a"\t"$0}' \
        > ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Count_mapped_read.txt

    # Count unmapped read
    echo "Count unmapped read : ${BAM_input}"
    samtools view -c -f 4 -@ ${Number_CPU} ${BAM_input} \
        | awk -v a="${SampleID}" '{print a"\t"$0}' \
        > ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Count_unmapped_read.txt

    # Count mapped bases, Mean depth, Breadth of coverage
    echo "Calculating depth and mapped bases : ${BAM_input}"
    samtools depth -a -@ ${Number_CPU} ${BAM_input} > ${OUT_count}Count_Mapped/${SampleID}/tmp_depth.txt

    # Mapped bases
    awk -v a="${SampleID}" '{sum+=length($0)} END {print a"\t"sum}' ${OUT_count}Count_Mapped/${SampleID}/tmp_depth.txt \
        > ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Count_map_base.txt

    # Mean depth and Breadth of coverage
    awk -v a="${SampleID}" '{sum+=$3; if($3>0) total+=1; c++} END {
        print a"\t"sum/c > "'${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Mean_read_depth.txt";
        print a"\t"(total/c)*100 > "'${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Breadth_of_Coverage.txt"
    }' ${OUT_count}Count_Mapped/${SampleID}/tmp_depth.txt

    rm ${OUT_count}Count_Mapped/${SampleID}/tmp_depth.txt

    echo "Sample ${SampleID} done!"
}

export -f process_sample
export OUT_count
export Number_CPU

# Run in parallel
cat ${INFO_FILE} | parallel -j ${Parallel_Sample} --colsep ',' process_sample {1} ${Number_CPU} ${OUT_count}
