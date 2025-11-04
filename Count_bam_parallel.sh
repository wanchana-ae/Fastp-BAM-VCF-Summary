#!/usr/bin/env bash
set -euo pipefail
# Load samtools module

ml load samtools

INFO_FILE="/fs4/Mining_SNP/Borassus_flabellifer_Linn/Germplasm/name_bam.txt"
OUT_count="/fs4/Mining_SNP/Borassus_flabellifer_Linn/Germplasm/Report/"
Number_CPU=8        # threads per sample
Parallel_Sample=4   # number of samples to run in parallel

mkdir -p ${OUT_count}Count_Mapped/

# Example: list of files to check
FILES=(
    "${OUT_count}Count_Mapped/All_Count_map_base.txt"
    "${OUT_count}Count_Mapped/All_Count_mapped_read.txt"
    "${OUT_count}Count_Mapped/All_Count_unmapped_read.txt"
    "${OUT_count}Count_Mapped/All_Count_Mean_read_depth.txt"
    "${OUT_count}Count_Mapped/All_Breadth_of_Coverage.txt"
)

# Flag to track if any file exists
found_existing=false

# Check each file
for FILE in "${FILES[@]}"; do
    if [[ -f "$FILE" ]]; then
        echo "⚠️ File already exists: $FILE"
        found_existing=true
    fi
done

# If any existing files were found, stop the process
if $found_existing; then
    echo "❌ Please remove the above files before proceeding."
    exit 1
else
    echo "✅ No existing files found — safe to continue."
fi

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

	###############################################################################################	
	mkdir -p ${OUT_count}Count_Mapped/${SampleID}
	#Count Mapped read
	echo "===================================================================================="
	echo "Count Mapped read / Depth average / Zero depth / Genome coverage ${SampleID} start"
	echo "===================================================================================="
	
	# Report overall each sample for check output
	echo "samtools flagstat : ${BAM_input}"
	samtools flagstat --threads ${Number_CPU} ${BAM_input} > ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_flagstat.txt & PID_BAM1=$!
	
	# Count Maped base
	echo "Count mapped base : ${BAM_input}"
	samtools view -@ ${Number_CPU} ${BAM_input} |cut -f 10 |tr -d \n |wc -c| awk -v a="${SampleID}" '{print a"\t"$0}' > ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Count_map_base.txt & PID_BAM2=$!
	
	# Count Mapped read 
	echo "Count mapped read : ${BAM_input}"
	samtools view -c -F 0xD04 -@ ${Number_CPU} ${BAM_input} | awk -v a="${SampleID}" '{print a"\t"$0}' > ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Count_mapped_read.txt & PID_BAM3=$!
	
	# Count Unmapped read
	echo "Count unmapped read : ${BAM_input}"
	samtools view -c -f 4 -@ ${Number_CPU} ${BAM_input} | awk -v a="${SampleID}" '{print a"\t"$0}' > ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Count_unmapped_read.txt & PID_BAM4=$!
		
	#Mean Read Depth
	echo "Mean read depth : ${BAM_input}"
	samtools depth -a ${BAM_input} | awk -v a="${SampleID}" '{sum+=$3} END {print a"\t"sum/NR}' > ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Mean_read_depth.txt & PID_BAM5=$!
		
	#Breadth of Coverage
	echo "Breadth of Coverage : ${BAM_input}"
	samtools depth -a ${BAM_input} | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'| awk -v a="${SampleID}" '{print a"\t"$0}' > ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Breadth_of_Coverage.txt & PID_BAM6=$!
	
	wait $PID_BAM1 $PID_BAM2 $PID_BAM3 $PID_BAM4 $PID_BAM5 $PID_BAM6

    #Merge file
    cat ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Count_map_base.txt >> ${OUT_count}Count_Mapped/All_Count_map_base.txt
    cat ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Count_mapped_read.txt >> ${OUT_count}Count_Mapped/All_Count_mapped_read.txt
    cat ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Count_unmapped_read.txt >> ${OUT_count}Count_Mapped/All_Count_unmapped_read.txt
    cat ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Mean_read_depth.txt >> ${OUT_count}Count_Mapped/All_Count_Mean_read_depth.txt
    cat ${OUT_count}Count_Mapped/${SampleID}/${SampleID}_Breadth_of_Coverage.txt >> ${OUT_count}Count_Mapped/All_Breadth_of_Coverage.txt

    echo "Sample ${SampleID} done!"
}

export -f process_sample
export OUT_count
export Number_CPU

# Run in parallel
cat ${INFO_FILE} | parallel -j ${Parallel_Sample} --colsep ',' process_sample {1} ${Number_CPU} ${OUT_count}

