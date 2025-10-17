import json
import pandas as pd
import glob
import os

# Path to folder containing fastp.json files
folder_path = "/fs4/Mining_SNP/Borassus_flabellifer_Linn/Germplasm/Trimed/"

# Automatically find all *_fastp.json files in the folder
json_files = glob.glob(os.path.join(folder_path, "*_fastp.json"))

# Store extracted information from each sample
rows = []

for json_file in json_files:
    with open(json_file) as f:
        data = json.load(f)

    # Extract sample name from file name
    sample_name = os.path.basename(json_file).replace("_fastp.json", "")

    # Get summary statistics before and after filtering
    before = data['summary']['before_filtering']
    after = data['summary']['after_filtering']

    # Raw bases and reads (เต็ม)
    raw_bases = before['total_bases']
    raw_reads = before['total_reads']

    # Raw bases (G) and reads (M)
    raw_bases_g = raw_bases / 1e9
    raw_reads_m = raw_reads / 1e6

    # Trimmed bases and reads (เต็ม)
    trimmed_bases = after['total_bases']
    trimmed_reads = after['total_reads']

    # Trimmed bases (G) and reads (M)
    trimmed_bases_g = trimmed_bases / 1e9
    trimmed_reads_m = trimmed_reads / 1e6

    # Effective rate = (bases after filtering / before) × 100
    effective = (trimmed_bases / raw_bases) * 100

    # Q20 and Q30 (percent)
    q20 = after['q20_rate'] * 100
    q30 = after['q30_rate'] * 100

    # Approximate error rate = 1 - Q30
    error_rate = (1 - after['q30_rate']) * 100

    # GC content percentage
    gc_content = after['gc_content']

    # Append the result for this sample
    rows.append({
        "Sample": sample_name,
        "Raw bases": raw_bases,
        "Raw bases (G)": f"{raw_bases_g:.2f}",
        "Raw reads": raw_reads,
        "Raw reads (M)": f"{raw_reads_m:.2f}",
        "Trimmed bases": trimmed_bases,
        "Trimmed bases (G)": f"{trimmed_bases_g:.2f}",
        "Trimmed reads": trimmed_reads,
        "Trimmed reads (M)": f"{trimmed_reads_m:.2f}",
        "Effective (%)": f"{effective:.2f}",
        "Error (%)": f"{error_rate:.2f}",
        "Q20 (%)": f"{q20:.2f}",
        "Q30 (%)": f"{q30:.2f}",
        "GC (%)": f"{gc_content:.2f}"
    })

# Convert to pandas DataFrame
df = pd.DataFrame(rows)

# Print the summary table
print(df.to_string(index=False))

# Save to CSV file
df.to_csv("fastp_summary_table.csv", index=False)
