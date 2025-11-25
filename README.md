# NGS-QC-Mapping-Stats

This repository contains scripts to generate quality control (QC) and mapping statistics for NGS (Next Generation Sequencing) data. It is designed for paired-end sequencing data and allows parallel processing to speed up computations.

## Contents

### 1. `extract_fastp_summary_v3.py`

- A Python script to summarize the QC results from fastp JSON output (*_fastp.json) for multiple samples.

- Generates a table with:
  - Raw/Trimmed bases and reads (in total, Giga/Mega units)
  - Effective rate (%)
  - Q20/Q30 (%)
  - Error rate (%)
  - GC content (%)

#### Usage

```bash
python extract_fastp_summary_v3.py
```

Make sure to set the correct folder path containing *_fastp.json files in the script.

### 2. `Count_bam_parallel.sh`

- A Bash script to process BAM files in parallel using GNU Parallel and Samtools.
- Computes for each sample:
  - Flagstat statistics
  - Mapped and unmapped read counts
  - Mapped bases
  - Mean read depth
  - Breadth of coverage

#### Usage

```bash
# Make the script executable
chmod +x Count_bam_parallel.sh

# Run the script
./Count_bam_parallel.sh
```

Make sure to update the paths inside the script:

- `INFO_FILE` : path to the list of BAM files
- `OUT_count` : output directory
- `Number_CPU` and Parallel_Sample for parallelization

### 3. `genotype_summary.sh`

**VCF 1/1 Genotype Summary Script**
This script summarizes the number of homozygous alternate (`1/1` or `1|1`) genotypes for both **SNPs** and **INDELs** in a compressed VCF file (`.vcf.gz`).

### 📁 File

`genotype_summary.sh`

### ⚙️ Requirements

- `bcftools`
- `awk`
- A valid `.vcf.gz` file (with or without `.tbi` index)

### 🧩 Usage

Change input file name VCF.gz

```bash
VCF_FILE="FF68_CN_170ea.vcf.gz"
```

Save and Run

```bash
# Make the script executable
chmod +x genotype_summary.sh

# Run the script
./genotype_summary.sh
```

🧠 What it does

1. Checks if the VCF index (.tbi) exists. If not, creates it automatically.

2. Extracts all sample names from the VCF.

3. Counts the number of homozygous alternate genotypes (1/1 or 1|1) separately for:

    - SNPs (length(REF)==1 && length(ALT)==1)

    - INDELs (other cases)

4. Saves the summary to genotype_summary.tsv in tab-separated format:

|Sample|SNP_1/1|INDEL_1/1|
|:-------|:-----:|:-----------:|
|SampleA|1532|213|
|SampleB| 1450|192|

🧾 Example output

```bash
Sample  SNP_1/1 INDEL_1/1
FF68_01  1532    210
FF68_02  1498    198
```

🧩 Notes

- Supports both unindexed and indexed `.vcf.gz` files.
- Works with multi-sample VCF files.
- Requires minimal system resources.

#### Requirements

- Bash >= 4.0
- GNU Parallel
- Samtools
- Python 3.8+
- pandas library (pip install pandas)

### 4. `SNP Density Pipeline`

This is a simple two-step pipeline for calculating and visualizing SNP density across a genome using VCF files and R-based plotting.
The workflow includes:

**1. Calculating SNP density using vcftools**  
**2. Plotting SNP density along each chromosome**  

#### 🔧 1. SNPdensity.sh — SNP Density Calculation

This script uses vcftools to compute SNP density in fixed-size windows (default: 1 Mb).
It then converts the output into a BED-like file for downstream visualization.

`Usage`

```bash
chmod +x SNPdensity.sh
./SNPdensity.sh input.vcf.gz
```

- $1: Input gzipped VCF file

**Workflow Overview**  
Based on the uploaded script:

- Loads the vcftools module

- Runs SNP density calculation using:  

```bash
vcftools --SNPdensity 1000000
```

change this if you need other window size.  

- Converts the .snpden output into a 5-column BED file:

```bash
cat snp_index_SNPdensity.snpden | \
grep -v "CHROM" | awk '{print "Chr"$1"\t"$2"\t"$2+1000000"\t"".""\t"$3}' > snp_index_SNPdensity.bed
```

*change 1000000 to same number with `vcftools --SNPdensity`

**Generated Output Files**  
`*.snpden — Raw SNP density from vcftools`  
`*.bed — BED-style file used for plotting SNP density`  

#### 🎨 2. Plot_SNPdensity.R — Visualization of SNP Density

This R script generates a chromosome-by-chromosome heatmap showing SNP density across the genome.

**Run on R and Rstudio**  

```bash  
Plot_SNPdensity.R 
```

**Plot Features**  

- Faceted heatmap (one panel per chromosome)
- SNP density represented as color intensity
- X-axis labeled in Mb
- Custom gradient color scheme using RColorBrewer

**Dependencies**  

```bash  
library(ggplot2)  
library(scales)  
library(RColorBrewer)  
```

📌 **Requirements Software**  

- vcftools
- bgzip and tabix (recommended for VCF handling)  
- R (>= 4.0)
