# NGS-QC-Mapping-Stats

This repository contains scripts to generate quality control (QC) and mapping statistics for NGS (Next Generation Sequencing) data. It is designed for paired-end sequencing data and allows parallel processing to speed up computations.

## Contents

### 1 . `extract_fastp_summary_v3.py`

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

#### Requirements

- Bash >= 4.0
- GNU Parallel
- Samtools
- Python 3.8+
- pandas library (pip install pandas)
