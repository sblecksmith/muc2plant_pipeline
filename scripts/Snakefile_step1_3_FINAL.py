# Snakefile
# for FL119 processing metagenomes 
###############################
# USAGE INSTRUCTIONS
# ------------------
# 1. Navigate to project root directory:
#    cd /path/to/project_root
# 
# 2. Copy config/config.yaml and edit paths for your system
#
# 3. Create sample sheet (sample_sheet.txt) with columns: sample_name, r1_path, r2_path
#    See example at bottom of this file
#
# 4. Get required files from MetaPhlAn github (this may need to be updated if databases change)
#    wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/sgb_to_gtdb_profile.py -O scripts/sgb_to_gtdb_profile.py
#    wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/util_fun.py -O scripts/util_fun.py
#    wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv -O scripts/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv
#    wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/merge_metaphlan_tables.py -O scripts/merge_metaphlan_tables.py  
#
# 5. Load snakemake v9.11.4 into a conda environment (if necessary):
#    eval "$(mamba shell hook --shell bash)"
#    mamba create -n snakemake_env -c conda-forge -c bioconda snakemake=9.11.4
#    conda activate snakemake_env
#
# 6. Install slurm executor plugin for snakemake v8+ (only needs to be done once):
#    pip install snakemake-executor-plugin-slurm==1.9.0
#
# 7. Quick check to make sure there are no errors (dry run):
#    snakemake -s scripts/Snakefile --configfile config.yaml -n
#
# 8. Run the pipeline using one of these methods:
#
#    METHOD A - Submit via sbatch script (recommended):
#    sbatch scripts/submit_snakefile.sh
#
#    METHOD B - Run directly with snakemake:
#    snakemake -s scripts/Snakefile --configfile config.yaml --executor slurm --jobs 20 --use-conda \
#        --default-resources slurm_account=GROUPNAME mem_mb=4096 runtime=600
#
# 9. Monitor progress:
#    tail -f logs/snakemake_<jobid>.out
###############################
# REQUIRED DIRECTORY STRUCTURE:
# Make sure the .py scripts from the MetaPhlAn github have execute permissions
# ls -l scripts/sgb_to_gtdb_profile.py 
# chmod +x scripts/sgb_to_gtdb_profile.py 
# 
# project_root/
# │
# ├── config/
# │   └── config.yaml                [REQUIRED - edit with your paths]
# │
# ├── scripts/
# │   ├── Snakefile                  [this file]
# │   ├── sgb_to_gtdb_profile.py     [required script from MetaPhlAn github]
# │   ├── util_fun.py                [required script from MetaPhlAn github]
# │   ├── mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv [required for sgb to gtdb script]
# │   ├── merge_metaphlan_tables.py  [required script from MetaPhlAn github]
# │   └── submit_snakefile.sh        [submit this file with sbatch]
# │
# ├── envs/
# │   ├── fastqc_multiqc.yaml
# │   ├── bowtie.yaml
# │   ├── fastp.yaml              
# │   └── metaphlan.yaml  
# │
# ├── sample_sheet.txt              [REQUIRED - tab-separated file with sample info]
# │
# └── fastq_files/                  [your input files]
#     ├── sample1_R1.fastq.gz
#     ├── sample1_R2.fastq.gz
#     ├── sample2_R1.fastq.gz
#     └── sample2_R2.fastq.gz
#
# All output directories will be created automatically by Snakemake
#
###############################
# SAMPLE SHEET FORMAT (sample_sheet.txt):
# Tab-separated file with these columns:
# sample_name	long_sample	r1_path	r2_path
# 109_C1_E	109_C1_E_S18_L006	fastq_files/109_C1_E_S18_L006_R1_001.fastq.gz	fastq_files/109_C1_E_S18_L006_R2_001.fastq.gz
# 109_C1_N	109_C1_N_S66_L006	fastq_files/109_C1_N_S66_L006_R1_001.fastq.gz	fastq_files/109_C1_N_S66_L006_R2_001.fastq.gz
#
# To auto-generate from fastq_files directory, run:
# echo -e "sample_name\tlong_samples\tr1_path\tr2_path" > sample_sheet.txt
# for r1 in fastq_files/*_R1_001.fastq.gz; do
#     r2="${r1/_R1_/_R2_}"
#     sample=$(basename "$r1" | sed 's/_S[0-9]*_L[0-9]*_R1_001.fastq.gz//')
#     long_sample=$(basename "$r1" | sed 's/_R1_001.fastq.gz//')
#     echo -e "${sample}\t${long_sample}\t${r1}\t${r2}" >> sample_sheet.txt
# done
###############################
#
# CONFIGURATION
import pandas as pd

# Load configuration file
configfile: "config/config.yaml"

# Path to sample sheet (tab-separated: sample_name, r1_path, r2_path)
SAMPLE_SHEET = config["sample_sheet"]

# Read sample sheet and define samples
samples_df = pd.read_csv(SAMPLE_SHEET, sep="\t")
SAMPLES = samples_df['sample_name'].tolist()

# fastqc outputs long sample names, define that here
LONG_SAMPLES = samples_df['long_sample'].tolist()

# Create dictionaries for easy lookup
SAMPLE_R1 = dict(zip(samples_df['sample_name'], samples_df['r1_path']))
SAMPLE_R2 = dict(zip(samples_df['sample_name'], samples_df['r2_path']))

# Get configuration values
PROJECT_NAME = config["project_name"]
INPUT_DIR = config["input_dir"]

# Get MetaPhlAn version for output naming
METAPHLAN_VERSION = config["versions"]["metaphlan"]

# Debug: print detected samples
print(f"DEBUG: Detected {len(SAMPLES)} samples from {SAMPLE_SHEET}")
if len(SAMPLES) > 0:
    print(f"First 3 samples: {SAMPLES[:3]}")

# Define GB for downstream resource allocation. Do not change this
GB = 1024

# Snakemake works backwards. Specify the final files you want to check for
# if you don't specify all of the required output files, some rules will not run
rule all:
    input:
        # This will trigger the entire fastp pipeline
        expand("step2_fastp/merged/{sample}.merged.fastq.gz", sample=SAMPLES),
        # This will trigger fastqc + multiqc
        "fastqc_summary/multiqc_report.html",
        
        # This will trigger the summary reports
        "step2_fastp/pairedend/reports/pairedend_summary.txt",
        "step2_fastp/merged/reports/merged_summary.txt",
        "step1_bowtie2/human_read_loss_summary.txt",

        # This will trigger MetaPhlAn steps
        f"step3_metaphlan/{PROJECT_NAME}_sgb.metaphlan{METAPHLAN_VERSION}.txt",
        expand(f"step3_metaphlan/gtdb_profile/{{sample}}.gtdb.metaphlan{METAPHLAN_VERSION}.txt", sample=SAMPLES)

########################################
# Step 0: fastqc
########################################
rule fastqc:
    input:
        r1=f"{INPUT_DIR}/{{long_samples}}_R1_001.fastq.gz",
        r2=f"{INPUT_DIR}/{{long_samples}}_R2_001.fastq.gz"
    output:
        html1="fastqc_output/{long_samples}_R1_001_fastqc.html",
        zipped1="fastqc_output/{long_samples}_R1_001_fastqc.zip",
        html2="fastqc_output/{long_samples}_R2_001_fastqc.html",
        zipped2="fastqc_output/{long_samples}_R2_001_fastqc.zip"
    #threads: config["resources"]["fastqc"]["threads"]
    resources:
        mem_mb=config["resources"]["fastqc"]["mem_gb"]*GB,
        runtime=config["resources"]["fastqc"]["runtime_min"],
        threads=config["resources"]["fastqc"]["threads"]
    conda:
        config["conda_envs"]["fastqc_multiqc"]
    shell:
        """ 
        fastqc -o fastqc_output -f fastq -t {resources.threads} {input.r1} {input.r2}
        """

########################################
# Step 0.1: multiqc
########################################
rule multiqc:
    input:
        expand("fastqc_output/{long_samples}_R1_001_fastqc.html", long_samples=LONG_SAMPLES),
        expand("fastqc_output/{long_samples}_R2_001_fastqc.html", long_samples=LONG_SAMPLES)
    output:
        "fastqc_summary/multiqc_report.html"
    conda:
        config["conda_envs"]["fastqc_multiqc"]
    shell:
        """
        multiqc fastqc_output -o fastqc_summary
        """

#######################################
# Step 1: Remove human reads
# run bowtie2-build command if needed, only need to run once
# Example: bowtie2-build /path/to/human_genome.fna /path/to/human_genome
#######################################
rule remove_human:
    input:
        r1=lambda wildcards: SAMPLE_R1[wildcards.sample],
        r2=lambda wildcards: SAMPLE_R2[wildcards.sample]
    output:
        pe1="step1_bowtie2/{sample}_nohuman_pe.1.fastq",
        pe2="step1_bowtie2/{sample}_nohuman_pe.2.fastq",
        se="step1_bowtie2/{sample}_nohuman_se.fastq",
        stats="step1_bowtie2/{sample}_alignment_stats.txt"
    #threads: config["resources"]["bowtie2"]["threads"]
    params:
        index=config["databases"]["human_genome"],
        basename="step1_bowtie2/{sample}_nohuman_pe.fastq"
    resources:
        mem_mb=config["resources"]["bowtie2"]["mem_gb"]*GB,
        runtime=config["resources"]["bowtie2"]["runtime_min"],
        threads=config["resources"]["bowtie2"]["threads"]
    conda:
        config["conda_envs"]["bowtie"]
    shell:
        """
        # run bowtie2
        bowtie2 -q -p {resources.threads} -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            --un-conc {params.basename} \
            --un {output.se} \
            1>/dev/null 2>{output.stats}
        """

########################################
# Step 1.1: Zip reads to save space
########################################
rule zip_reads:
    input:
        pe1="step1_bowtie2/{sample}_nohuman_pe.1.fastq",
        pe2="step1_bowtie2/{sample}_nohuman_pe.2.fastq",
        se="step1_bowtie2/{sample}_nohuman_se.fastq"
    output:
        pe1="step1_bowtie2/{sample}_nohuman_pe.1.fastq.gz",
        pe2="step1_bowtie2/{sample}_nohuman_pe.2.fastq.gz",
        se="step1_bowtie2/{sample}_nohuman_se.fastq.gz"
    shell:
        """ 
        gzip {input.pe1} {input.pe2} {input.se}
        """

########################################
# Step 1.2: Summarize reads lost (per sample)
########################################
rule summarize_human_read_loss_per_sample:
    input:
        nh="step1_bowtie2/{sample}_nohuman_pe.1.fastq.gz",
        h=lambda wildcards: SAMPLE_R1[wildcards.sample]
    output:
        temp("step1_bowtie2/{sample}_read_loss.txt")
    resources:
        mem_mb=2*GB,
        runtime=30
    run:
        import gzip
        
        # Count total reads before filtering (from R1 only)
        total_reads = 0
        with gzip.open(input.h, 'rt') as f:
            total_reads = sum(1 for line in f) // 4
        
        # Count remaining reads after filtering (from R1 only)
        remaining_reads = 0
        with gzip.open(input.nh, 'rt') as f:
            remaining_reads = sum(1 for line in f) // 4
        
        # Calculate percentage
        percent_lost = 100 * (total_reads - remaining_reads) / total_reads if total_reads > 0 else 0.0
        
        # Write results for this sample
        with open(output[0], "w") as out:
            out.write(f"{wildcards.sample}\t{total_reads}\t{remaining_reads}\t{percent_lost:.2f}\n")

########################################
# Step 1.3: Aggregate all sample summaries
########################################
rule aggregate_human_read_loss:
    input:
        expand("step1_bowtie2/{sample}_read_loss.txt", sample=SAMPLES)
    output:
        "step1_bowtie2/human_read_loss_summary.txt"
    resources:
        mem_mb=1*GB,
        runtime=5
    shell:
        """
        echo -e "sample\ttotal_reads\tremaining_reads\tpercent_lost" > {output}
        cat {input} >> {output}
        """

########################################
# Step 2: Quality trimming with fastp
# This first pass saves paired end reads 
########################################
rule trim_quality:
    input:
        pe1="step1_bowtie2/{sample}_nohuman_pe.1.fastq.gz",
        pe2="step1_bowtie2/{sample}_nohuman_pe.2.fastq.gz"
    output:
        r1="step2_fastp/pairedend/{sample}_R1_paired.fastq.gz",
        r2="step2_fastp/pairedend/{sample}_R2_paired.fastq.gz",
        ur1="step2_fastp/pairedend/{sample}_R1_unpaired.fastq.gz",
        ur2="step2_fastp/pairedend/{sample}_R2_unpaired.fastq.gz",
        html="step2_fastp/pairedend/reports/{sample}.fastp.html",
        json="step2_fastp/pairedend/reports/{sample}.fastp.json"
    params:
        min_length=config["fastp_params"]["trim"]["min_length"],
        quality_phred=config["fastp_params"]["trim"]["quality_phred"]
    conda:
        config["conda_envs"]["fastp"]
    #threads: config["resources"]["fastp_trim"]["threads"]
    resources:
        mem_mb=config["resources"]["fastp_trim"]["mem_gb"]*GB,
        runtime=config["resources"]["fastp_trim"]["runtime_min"],
        threads=config["resources"]["fastp_trim"]["threads"]
    shell:
        """
        fastp --in1 {input.pe1} --in2 {input.pe2} \
            --thread {resources.threads} --detect_adapter_for_pe \
            --trim_poly_g --dedup \
            --length_required {params.min_length} \
            --qualified_quality_phred {params.quality_phred} \
            --unpaired1 {output.ur1} --unpaired2 {output.ur2} \
            --out1 {output.r1} --out2 {output.r2} \
            --json {output.json} \
            --html {output.html}
        """

########################################
# Step 2.1: Merge paired reads with fastp
########################################
rule merge_reads:
    input:
        r1="step2_fastp/pairedend/{sample}_R1_paired.fastq.gz",
        r2="step2_fastp/pairedend/{sample}_R2_paired.fastq.gz",
    output:
        merged="step2_fastp/merged/{sample}.merged.fastq.gz",
        um1="step2_fastp/merged/{sample}_R1_unmerged.fastq.gz",
        um2="step2_fastp/merged/{sample}_R2_unmerged.fastq.gz",
        html="step2_fastp/merged/reports/{sample}.fastp.html",
        json="step2_fastp/merged/reports/{sample}.fastp.json"
    params:
        overlap_len=config["fastp_params"]["merge"]["overlap_len_require"],
        overlap_diff=config["fastp_params"]["merge"]["overlap_diff_percent_limit"]
    conda:
        config["conda_envs"]["fastp"]
    resources:
        mem_mb=config["resources"]["fastp_merge"]["mem_gb"]*GB,
        runtime=config["resources"]["fastp_merge"]["runtime_min"],
        threads=config["resources"]["fastp_merge"]["threads"]
    #threads: config["resources"]["fastp_merge"]["threads"]
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
        --thread {resources.threads} \
        --merge \
        --merged_out {output.merged} \
        --out1 {output.um1} --out2 {output.um2} \
        --overlap_len_require {params.overlap_len} \
        --overlap_diff_percent_limit {params.overlap_diff} \
        --disable_trim_poly_g \
        --disable_adapter_trimming \
        --json {output.json} \
        --html {output.html}
        """

########################################
# Step 2.2: Summarize fastp reports
########################################
rule summarize_fastp_reports:
    input:
        pe_reports=expand("step2_fastp/pairedend/reports/{sample}.fastp.json", sample=SAMPLES),
        merged_reports=expand("step2_fastp/merged/reports/{sample}.fastp.json", sample=SAMPLES)
    output:
        pe_summary="step2_fastp/pairedend/reports/pairedend_summary.txt",
        merged_summary="step2_fastp/merged/reports/merged_summary.txt"
    run:
        import json
        # -----------------------------
        # Paired-end summary
        # -----------------------------
        with open(output.pe_summary, "w") as out_pe:
            # write header
            out_pe.write("sample_id\tbefore_filtering_total_reads\tafter_filtering_total_reads\tlow_quality_reads\ttoo_many_N_reads\ttoo_short_reads\tduplication_rate\tinsert_peak\tpercent_lost\n")
            
            for sample in SAMPLES:
                report_file = f"step2_fastp/pairedend/reports/{sample}.fastp.json"
                with open(report_file) as f:
                    data = json.load(f)

                # Extract fields from the JSON
                total_before = data["summary"]["before_filtering"]["total_reads"]
                total_after  = data["summary"]["after_filtering"]["total_reads"]
                low_qual     = data["filtering_result"]["low_quality_reads"]
                too_many_N   = data["filtering_result"]["too_many_N_reads"]
                too_short    = data["filtering_result"]["too_short_reads"]
                dup_rate     = data["duplication"]["rate"]
                insert_peak  = data["insert_size"]["peak"]
                percent_lost = 100 * (total_before - total_after) / total_before if total_before > 0 else 0

                # write line for this sample
                out_pe.write(f"{sample}\t{total_before}\t{total_after}\t{low_qual}\t{too_many_N}\t{too_short}\t{dup_rate}\t{insert_peak}\t{percent_lost}\n")

        # -----------------------------
        # Merged summary
        # -----------------------------
        with open(output.merged_summary, "w") as out_merged:
            out_merged.write("sample_id\tbefore_filtering_total_reads\tafter_filtering_total_reads\tlow_quality_reads\ttoo_many_N_reads\ttoo_short_reads\tduplication_rate\tinsert_peak\tpercent_lost\n")

            for sample in SAMPLES:
                report_file = f"step2_fastp/merged/reports/{sample}.fastp.json"
                with open(report_file) as f:
                    data = json.load(f)

                total_before = data["summary"]["before_filtering"]["total_reads"] / 2
                total_after  = data["summary"]["after_filtering"]["total_reads"]
                low_qual     = data["filtering_result"]["low_quality_reads"]
                too_many_N   = data["filtering_result"]["too_many_N_reads"]
                too_short    = data["filtering_result"]["too_short_reads"]
                dup_rate     = data["duplication"]["rate"]
                insert_peak  = data["insert_size"]["peak"]
                percent_lost = 100 * (total_before - total_after) / total_before if total_before > 0 else 0

                out_merged.write(f"{sample}\t{total_before}\t{total_after}\t{low_qual}\t{too_many_N}\t{too_short}\t{dup_rate}\t{insert_peak}\t{percent_lost}\n")

########################################
# Step 3: MetaPhlAn
# written for MetaPhlAn 4.2.2
# make sure to specify the index and check if any new databases are available 
########################################
rule metaphlan:
    input:
        r1="step2_fastp/pairedend/{sample}_R1_paired.fastq.gz",
        r2="step2_fastp/pairedend/{sample}_R2_paired.fastq.gz"
    output:
        bt="step3_metaphlan/bowtie_outs/{sample}.metaphlan.bz2",
        sams="step3_metaphlan/sams/{sample}.metaphlan.sam",
        profile=f"step3_metaphlan/sgb_profile/{{sample}}.sgb.metaphlan{METAPHLAN_VERSION}.txt"

    params:
        db_dir=config["databases"]["metaphlan_db"],
        index=config["databases"]["metaphlan_index"]
    conda:
        config["conda_envs"]["metaphlan"]
    #threads: config["resources"]["metaphlan"]["threads"]
    resources:
        mem_mb=config["resources"]["metaphlan"]["mem_gb"]*GB,
        runtime=config["resources"]["metaphlan"]["runtime_min"],
        threads=config["resources"]["metaphlan"]["threads"]
    shell:
        """
        metaphlan {input.r1},{input.r2} --input_type fastq \
        --db_dir {params.db_dir} \
        --nproc {resources.threads} \
        --index {params.index} \
        --mapout {output.bt} -s {output.sams} \
        -o {output.profile} --verbose --offline 
        """

########################################
# Step 3.1: GTDB taxonomy
# python script from MetaPhlAn github
# this script is specific to the index used in rule metaphlan,
# so there will probably be an updated script if a new index is released 
# wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/sgb_to_gtdb_profile.py -O scripts/sgb_to_gtdb_profile.py
# wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/util_fun.py -O scripts/util_fun.py
# wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv -O scripts/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv
########################################
rule gtdb:
    input: 
        profile=f"step3_metaphlan/sgb_profile/{{sample}}.sgb.metaphlan{METAPHLAN_VERSION}.txt"
    output: 
        profile2=f"step3_metaphlan/gtdb_profile/{{sample}}.gtdb.metaphlan{METAPHLAN_VERSION}.txt"
    shell:
        """
        python scripts/sgb_to_gtdb_profile.py -i {input.profile} -o {output.profile2}
        """

########################################
# Step 3.2: Combine metaphlan
# python script from MetaPhlAn github
# this script is specific to the index used in rule metaphlan,
# so there will probably be an updated script if a new index is released 
# wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/merge_metaphlan_tables.py -O scripts/merge_metaphlan_tables.py
########################################
rule combine_metaphlan_sgb:
    input:
        expand(f"step3_metaphlan/sgb_profile/{{sample}}.sgb.metaphlan{METAPHLAN_VERSION}.txt", sample=SAMPLES)
    output:
        f"step3_metaphlan/{PROJECT_NAME}_sgb.metaphlan{METAPHLAN_VERSION}.txt"
    shell:
        """
        python scripts/merge_metaphlan_tables.py {input} > {output}
        """