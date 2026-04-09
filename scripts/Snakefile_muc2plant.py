# Snakefile for muc2plant calculation using dbcan
# built from metagenomic pipeline by Nithya K Kumar, UC Davis
# https://github.com/nithyak2/lemaylab_metagenomics_pipeline
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
# 4. Load snakemake v9.11.4 into a conda environment (if necessary):
#    eval "$(mamba shell hook --shell bash)"
#    mamba create -n snakemake_env -c conda-forge -c bioconda snakemake=9.11.4
#    conda activate snakemake_env
#
# 5. Install slurm executor plugin for snakemake v8+ (only needs to be done once):
#    pip install snakemake-executor-plugin-slurm
#
# 6. Quick check to make sure there are no errors (dry run):
#    snakemake -s scripts/Snakefile --configfile config.yaml -n
#
# 7. Run the pipeline using one of these methods:
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
#
#
###############################
# REQUIRED DIRECTORY STRUCTURE:
# project_root/
# │
# ├── config/
# │   └── config.yaml                [REQUIRED - edit with your paths]
# │
# ├── scripts/
# │   ├── Snakefile_muc2plant.py     [this file]
# │   ├── submit_snakefile.sh        [submit this file with sbatch]
# │   ├── make_sample_sheet
# │   └── calculate_muc2plant.R       
# │
# ├── envs/
# │   ├── fastqc_multiqc.yaml
# │   ├── bowtie.yaml
# │   ├── fastp.yaml              
# │   ├── metaphlan.yaml
# │   ├── megahit.yaml
# │   ├── pyrodigal.yaml
# │   ├── bwa.yaml
# │   ├── samtools.yaml
# │   ├── rundbcan.yaml
# │   └── dbcan_utils.yaml 
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
configfile: "../config/config_muc2plant.yaml"

# Path to sample sheet (tab-separated: sample_name, r1_path, r2_path)
SAMPLE_SHEET = config["sample_sheet"]

# Read sample sheet and define samples
samples_df = pd.read_csv(SAMPLE_SHEET, sep=",")
SAMPLES = samples_df['sample_name'].tolist()

# fastqc outputs long sample names, define that here
LONG_SAMPLES = samples_df['long_sample'].tolist()

# Create dictionaries for easy lookup
SAMPLE_R1 = dict(zip(samples_df['sample_name'], samples_df['r1_path']))
SAMPLE_R2 = dict(zip(samples_df['sample_name'], samples_df['r2_path']))

# Get configuration values
PROJECT_NAME = config["project_name"]
INPUT_DIR = config["input_dir"]

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
        expand("fastp_output/merged/{sample}.merged.fastq.gz", sample=SAMPLES),
        
        # This will trigger fastqc + multiqc
        "fastqc_summary/multiqc_report.html",
        
        # This will trigger the summary reports
        "fastp_output/pairedend/reports/pairedend_summary.txt",
        "fastp_output/merged/reports/merged_summary.txt",
        "fastp_output/human_read_loss_summary.txt",
        
        # This will trigger megahit
        expand("megahit_output/{sample}.contigs.fa.gz", sample=SAMPLES),
        
        # This will trigger pyrodigal
        expand("pyrodigal_output/{sample}.faa.gz", sample=SAMPLES),
        
        # This will trigger dbcan
        expand("rundbcan_output/{sample}_dbCAN/PUL_blast.out", sample=SAMPLES),
        
        # This will trigger bwa for contigs
        expand("bwa_output/{sample}.bam.bai", sample=SAMPLES),
        
        # This will trigger cal_coverage
        expand("rundbcan_output/{sample}_abund/{sample}.depth.txt", sample=SAMPLES),
        
        # This will trigger fam_abund
        expand("rundbcan_output/{sample}_abund/fam_abund.out", sample=SAMPLES),
        
        # This will trigger calculate_muc2plant
        expand("muc2plant.tsv),
        
########################################
# fastqc
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
    threads: config["resources"]["fastqc"]["threads"]
    resources:
        mem_mb=config["resources"]["fastqc"]["mem_gb"]*GB,
        runtime=config["resources"]["fastqc"]["runtime_min"]
    conda:
        config["conda_envs"]["fastqc_multiqc"]
    shell:
        """ 
        fastqc -o fastqc_output -f fastq -t {threads} {input.r1} {input.r2}
        """

########################################
# multiqc
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
        pe1="bowtie2_output/{sample}_nohuman_pe.1.fastq",
        pe2="bowtie2_output/{sample}_nohuman_pe.2.fastq",
        se="bowtie2_output/{sample}_nohuman_se.fastq",
        stats="bowtie2_output/{sample}_alignment_stats.txt"
    threads: config["resources"]["bowtie2"]["threads"]
    params:
        index=config["databases"]["human_genome"],
        basename="bowtie2_output/{sample}_nohuman_pe.fastq"
    resources:
        mem_mb=config["resources"]["bowtie2"]["mem_gb"]*GB,
        runtime=config["resources"]["bowtie2"]["runtime_min"]
    conda:
        config["conda_envs"]["bowtie"]
    shell:
        """
        # run bowtie2
        bowtie2 -q -p {threads} -x {params.index} \
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
        pe1="bowtie2_output/{sample}_nohuman_pe.1.fastq",
        pe2="bowtie2_output/{sample}_nohuman_pe.2.fastq",
        se="bowtie2_output/{sample}_nohuman_se.fastq"
    output:
        pe1="bowtie2_output/{sample}_nohuman_pe.1.fastq.gz",
        pe2="bowtie2_output/{sample}_nohuman_pe.2.fastq.gz",
        se="bowtie2_output/{sample}_nohuman_se.fastq.gz"
    shell:
        """ 
        gzip {input.pe1} {input.pe2} {input.se}
        """

########################################
# Step 1.2: Summarize reads lost (per sample)
########################################
rule summarize_human_read_loss_per_sample:
    input:
        nh="bowtie2_output/{sample}_nohuman_pe.1.fastq.gz",
        h=lambda wildcards: SAMPLE_R1[wildcards.sample]
    output:
        temp("bowtie2_output/{sample}_read_loss.txt")
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
        expand("bowtie2_output/{sample}_read_loss.txt", sample=SAMPLES)
    output:
        "bowtie2_output/human_read_loss_summary.txt"
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
        pe1="bowtie2_output/{sample}_nohuman_pe.1.fastq.gz",
        pe2="bowtie2_output/{sample}_nohuman_pe.2.fastq.gz"
    output:
        r1="fastp_output/pairedend/{sample}_R1_paired.fastq.gz",
        r2="fastp_output/pairedend/{sample}_R2_paired.fastq.gz",
        ur1="fastp_output/pairedend/{sample}_R1_unpaired.fastq.gz",
        ur2="fastp_output/pairedend/{sample}_R2_unpaired.fastq.gz",
        html="fastp_output/pairedend/reports/{sample}.fastp.html",
        json="fastp_output/pairedend/reports/{sample}.fastp.json"
    params:
        min_length=config["fastp_params"]["trim"]["min_length"],
        quality_phred=config["fastp_params"]["trim"]["quality_phred"]
    conda:
        config["conda_envs"]["fastp"]
    threads: config["resources"]["fastp_trim"]["threads"]
    resources:
        mem_mb=config["resources"]["fastp_trim"]["mem_gb"]*GB,
        runtime=config["resources"]["fastp_trim"]["runtime_min"]
    shell:
        """
        fastp --in1 {input.pe1} --in2 {input.pe2} \
            --thread {threads} --detect_adapter_for_pe \
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
        r1="fastp_output/pairedend/{sample}_R1_paired.fastq.gz",
        r2="fastp_output/pairedend/{sample}_R2_paired.fastq.gz",
    output:
        merged="fastp_output/merged/{sample}.merged.fastq.gz",
        um1="fastp_output/merged/{sample}_R1_unmerged.fastq.gz",
        um2="fastp_output/merged/{sample}_R2_unmerged.fastq.gz",
        html="fastp_output/merged/reports/{sample}.fastp.html",
        json="fastp_output/merged/reports/{sample}.fastp.json"
    params:
        overlap_len=config["fastp_params"]["merge"]["overlap_len_require"],
        overlap_diff=config["fastp_params"]["merge"]["overlap_diff_percent_limit"]
    conda:
        config["conda_envs"]["fastp"]
    resources:
        mem_mb=config["resources"]["fastp_merge"]["mem_gb"]*GB,
        runtime=config["resources"]["fastp_merge"]["runtime_min"]
    threads: config["resources"]["fastp_merge"]["threads"]
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
        --thread {threads} \
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
        pe_reports=expand("fastp_output/pairedend/reports/{sample}.fastp.json", sample=SAMPLES),
        merged_reports=expand("fastp_output/merged/reports/{sample}.fastp.json", sample=SAMPLES)
    output:
        pe_summary="fastp_output/pairedend/reports/pairedend_summary.txt",
        merged_summary="fastp_output/merged/reports/merged_summary.txt"
    run:
        import json
        # -----------------------------
        # Paired-end summary
        # -----------------------------
        with open(output.pe_summary, "w") as out_pe:
            # write header
            out_pe.write("sample_id\tbefore_filtering_total_reads\tafter_filtering_total_reads\tlow_quality_reads\ttoo_many_N_reads\ttoo_short_reads\tduplication_rate\tinsert_peak\tpercent_lost\n")
            
            for sample in SAMPLES:
                report_file = f"fastp_output/pairedend/reports/{sample}.fastp.json"
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
# megahit
########################################
rule megahit:
    input:
        r1="/quobyte/dglemaygrp/sblecksmith/filterbytile_test/{sample}_R1_filtered.fastq.gz",
        r2="/quobyte/dglemaygrp/sblecksmith/filterbytile_test/{sample}_R2_filtered_synced.fastq.gz",
    output:
        contigs = "megahit_output/{sample}.contigs.fa.gz",
    conda: config["conda_envs"]["megahit"]
    params:
        prefix= "{sample}",
    threads: config["resources"]["megahit"]["threads"]
    shadow: "minimal"
    resources:
        mem_mb=config["resources"]["megahit"]["mem_gb"]*GB,
        runtime=config["resources"]["megahit"]["runtime_min"]
    shell:
        """ 
        megahit -m 0.5 \
        -t {threads} \
        -1 {input.r1} \
        -2 {input.r2} \
        --out-prefix {params.prefix} \
        --min-contig-len 1000
        
        pigz \
        --no-name \
        -p {threads} \
        megahit_out/*.fa 
        
        mv megahit_out/{params.prefix}.contigs.fa.gz {output.contigs}
        
        """

########################################
# pyrodigal
########################################
rule pyrodigal:
    input:
      contig="megahit_output/{sample}.contigs.fa.gz",
    output:
       score="pyrodigal_output/{sample}.score.gz",
       fna="pyrodigal_output/{sample}.fna.gz",
       faa="pyrodigal_output/{sample}.faa.gz",
       gff="pyrodigal_output/{sample}.gff.gz",
    params:
       prefix= "{sample}"
    shadow: "minimal"
    threads: config["resources"]["pyrodigal"]["threads"]
    conda: config["conda_envs"]["pyrodigal"]
    resources:
        mem_mb=config["resources"]["pyrodigal"]["mem_gb"]*GB,
        runtime=config["resources"]["pyrodigal"]["runtime_min"]
    shell:
        """
        pigz -cdf {input.contig} > pigz_fasta.fna
        
        pyrodigal -p meta \
        -i pigz_fasta.fna \
        -d {params.prefix}.fna \
        -a {params.prefix}.faa \
        -f gff \
        -o {params.prefix}.gff \
        -j {threads} \
        -s {params.prefix}.score
        
        pigz -nmf {params.prefix}*
        
        mv {params.prefix}.fna.gz {output.fna}
        mv {params.prefix}.faa.gz {output.faa}
        mv {params.prefix}.gff.gz {output.gff}
        mv {params.prefix}.score.gz {output.score}
        """

#######################################
# run_dbcan 
#######################################
rule run_dbcan:
    input:
        faa="pyrodigal_output/{sample}.faa.gz",
        gff="pyrodigal_output/{sample}.gff.gz",
    output:
        pul="rundbcan_output/{sample}_dbCAN/PUL_blast.out",
        cgc_faa="rundbcan_output/{sample}_dbCAN/CGC.faa",
        cgc_gff="rundbcan_output/{sample}_dbCAN/cgc.gff",
        cgc_std="rundbcan_output/{sample}_dbCAN/cgc_standard_out.tsv",
        cgc_std_sum="rundbcan_output/{sample}_dbCAN/cgc_standard_out_summary.tsv",
        diamond="rundbcan_output/{sample}_dbCAN/diamond.out",
        diamond_tc="rundbcan_output/{sample}_dbCAN/diamond.out.tc",
        diamond_tf="rundbcan_output/{sample}_dbCAN/diamond.out.tf",
        diamond_pep="rundbcan_output/{sample}_dbCAN/diamond.out.peptidase",
        diamond_sulf="rundbcan_output/{sample}_dbCAN/diamond.out.sulfatase",
        hmm="rundbcan_output/{sample}_dbCAN/dbCAN_hmm_results.tsv",
        sub_hmm="rundbcan_output/{sample}_dbCAN/dbCANsub_hmm_results.tsv",
        sub_hmm_raw="rundbcan_output/{sample}_dbCAN/dbCANsub_hmm_raw.tsv",
        overview="rundbcan_output/{sample}_dbCAN/overview.tsv",
        stp_hmm="rundbcan_output/{sample}_dbCAN/STP_hmm_results.tsv",
        substrate="rundbcan_output/{sample}_dbCAN/substrate_prediction.tsv",
        total_cgc="rundbcan_output/{sample}_dbCAN/total_cgc_info.tsv",
        uni="rundbcan_output/{sample}_dbCAN/uniInput.faa",
    threads: config["resources"]["rundbcan"]["threads"]
    conda: config["conda_envs"]["rundbcan"]
    shadow: "minimal"
    log:
        "logs/rundbcan/{sample}.log"
    params:
        outdir="{sample}_dbCAN",
        db_dir=config["databases"]["dbcan_db"],
    resources:
        mem_mb=config["resources"]["rundbcan"]["mem_gb"]*GB,
        runtime=config["resources"]["rundbcan"]["runtime_min"]
    shell:
        """
        gunzip -c {input.faa} > {wildcards.sample}.faa
        gunzip -c {input.gff} > {wildcards.sample}.gff
        
        run_dbcan easy_substrate \
        --input_raw_data {wildcards.sample}.faa \
        --db_dir {params.db_dir} \
        --mode protein  \
        --input_gff {wildcards.sample}.gff \
        --gff_type prodigal \
        --output_dir {params.outdir} \
        2>&1 | tee {log}
        
        sleep 10
        mv {params.outdir} rundbcan_output/
        """


########################################
# bwa_contigs - Read mapping to all contigs of each sample
########################################

rule bwa_sam_contigs:
    input:
        r1="/quobyte/dglemaygrp/sblecksmith/filterbytile_test/{sample}_R1_filtered.fastq.gz",
        r2="/quobyte/dglemaygrp/sblecksmith/filterbytile_test/{sample}_R2_filtered_synced.fastq.gz",
        contig="megahit/{sample}.contigs.fa.gz",
    output:
        bam="bwa_output/{sample}.bam",
        bai="bwa_output/{sample}.bam.bai",
        amb="bwa_output/bwa_{sample}/{sample}.amb",
        ann="bwa_output/bwa_{sample}/{sample}.ann",
        bwt="bwa_output/bwa_{sample}/{sample}.bwt",
        pac="bwa_output/bwa_{sample}/{sample}.pac",
        sa="bwa_output/bwa_{sample}/{sample}.sa",
    conda: config["conda_envs"]["bwa"],
    params:
        index_prefix="bwa/bwa_{sample}/{sample}",
    resources:
        mem_mb=config["resources"]["bwa"]["mem_gb"]*GB,
        runtime=config["resources"]["bwa"]["runtime_min"],
    threads: config["resources"]["bwa"]["threads"]
    shadow: "minimal"
    shell:
        """
        gunzip -c {input.contig} > {wildcards.sample}.fa
        mkdir -p bwa/bwa_{wildcards.sample} samtools
        
        bwa index -p {params.index_prefix} {wildcards.sample}.fa
        
        bwa mem \
            -t {threads} \
            {params.index_prefix} \
            {input.r1} {input.r2} \
        | samtools sort \
            -@ {threads} \
            -m 4G \
            -o {output.bam} \
            - 

        samtools index -@ {threads} {output.bam}
        """

#######################################
# cal_coverage - count reads per gene from a BAM and GFF
########################################
rule cal_coverage:
    input:
        gff="pyrodigal_output/{sample}.gff.gz",
        bam="bwa_output/{sample}.bam",
        bai="bwa_output/{sample}.bam.bai",
    output:
        depth="rundbcan_output/{sample}_abund/{sample}.depth.txt",
        #depth_dir=directory("samtools/{sample}_samtools_depth")
    conda: config["conda_envs"]["dbcan_utils"]
    threads: config["resources"]["dbcan_utils"]["threads"]
    shadow: "minimal"
    shell:
        """
        gunzip -c {input.gff} > {wildcards.sample}.gff
        
        dbcan_utils cal_coverage \
        -g {wildcards.sample}.gff \
        -i {input.bam} \
        -o {output.depth} \
        -t {threads} \
        --overlap_base_ratio 0.2 \
        --mapping_quality 30 \
        --identity 0.98\
        """


#######################################
# fam_abund
########################################
rule fam_abund:
    input:
        depth="rundbcan_output/{sample}_abund/{sample}.depth.txt",
        overview="rundbcan_output/{sample}_dbCAN/overview.tsv"
    output:
        fam="rundbcan_output/{sample}_abund/fam_abund.out",
        ec="rundbcan_output/{sample}_abund/EC_abund.out", 
        subfam="rundbcan_output/{sample}_abund/subfam_abund.out",
    params:
        outdir=directory("{sample}_dbCAN"),
        abunddir=directory("rundbcan_output/{sample}_abund"),
        depth="{sample}.depth.txt"
    conda: config["conda_envs"]["dbcan_utils"]
    shell:
        """
        cd {params.abunddir}
        dbcan_utils fam_abund -bt {params.depth} -i ../{params.outdir} -a TPM
        """        

#######################################
# aggregate_cazyme_families
########################################
import pandas as pd
from pathlib import Path

# Find all CAZyme files
files = Path("rundbcan_output").glob("*_abund/fam_abund.out")

dfs = []
for f in files:

    # Extract sample name
    sample_id = f.parent.name.replace("_abund", "")

    # Read the file
    df = pd.read_csv(f, sep="\t")

    # Keep only Family and Abundance columns
    df = df[['Family', 'Abundance']]

    # Rename Abundance to sample name
    df = df.rename(columns={'Abundance': sample_id})

    # Set Family as index for merging
    df = df.set_index('Family')

    dfs.append(df)

# Merge all dataframes on Family (outer join to keep all families)
combined = pd.concat(dfs, axis=1, join='outer')

# Fill missing values with 0
combined = combined.fillna(0)

# Save with Family as first column
combined.to_csv("aggregated_cazyme_family_TPM.tsv", sep="\t")


#######################################
# muc2plant
########################################
rule calculate_muc2plant:
    input:
        cazyme_abund="aggregated_cazyme_family_TPM.tsv",
        fams="smits_cazyme_substrates_sep.csv",
    output:
        muc2plant="muc2plant.tsv"
    params:
    shell:
        """
        Rscript calculate_muc2plant.R --input {input.cazyme_abund} --families {input.fams} --output {output.muc2plant}
        
        """        

