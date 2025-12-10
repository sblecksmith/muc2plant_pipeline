# Lemay Lab Metagenomics Pipeline (Snakemake)

A reproducible Snakemake workflow for processing 150 bp PE shotgun metagenomic sequencing data.  
Developed by Nithya K Kumar (Lemay Lab, UC Davis, USDA).

---

## Overview
This pipeline streamlines the analysis of shotgun metagenomic sequencing data.  
It performs quality control, removes host reads, profiles microbial taxa, and summarizes results in a reproducible and modular framework.

**Core features:**
- Automated workflow management using **Snakemake**
- Reproducible environments using **Conda**
- Modular structure for easy extension
- Example configuration for quick setup

---

##  Workflow Summary

| Step | Tool | Description |
|------|------|--------------|
| 1. Quality control | FastQC | Assess read quality |
| 2. Host read removal | Bowtie2 | Remove host reads |
| 3. Trimming | Fastp | Trim adapters and low-quality bases |
| 4. Merging | Fastp | Merge paired end reads 
| 5. Taxonomic profiling | MetaPhlAn | Read based classification from PE reads |

## snakemake simplified DAG (rulegraph)
![main](pipeline_DAG.png)


---


## Clone the repository
```bash
git clone https://github.com/nithyak2/lemaylab_metagenomics_pipeline
cd metagenomics-pipeline
```


 USAGE INSTRUCTIONS
 ------------------
1.   Navigate to project root directory:
```bash
cd /path/to/project_root
```
 
2. Copy config/config.yaml and edit paths for your system. **This should be the only file you need to edit, you should not need to edit the Snakefile**

3. Create sample sheet (sample_sheet.txt) with columns: sample_name, r1_path, r2_path
```bash
# To auto-generate from fastq_files directory, run:
# edit this based on your file names 
# Create the sample sheet from your existing files
 echo -e "sample_name\tlong_sample\tr1_path\tr2_path" > sample_sheet.txt
 for r1 in fastq_files/*_R1_001.fastq.gz; do
     r2="${r1/_R1_/_R2_}"
     sample=$(basename "$r1" | sed 's/_S[0-9]*_L[0-9]*_R1_001.fastq.gz//')
     long_sample=$(basename "$r1" | sed 's/_R1_001.fastq.gz//')
     echo -e "${sample}\t${long_sample}\t${r1}\t${r2}" >> sample_sheet.txt
 done
 
# Example output:
# Tab-separated file with these columns:
# sample_name	long_sample	r1_path	r2_path
# 109_C1_E	109_C1_E_S18_L006	fastq_files/109_C1_E_S18_L006_R1_001.fastq.gz	fastq_files/109_C1_E_S18_L006_R2_001.fastq.gz
# 109_C1_N	109_C1_N_S66_L006	fastq_files/109_C1_N_S66_L006_R1_001.fastq.gz	fastq_files/109_C1_N_S66_L006_R2_001.fastq.gz
```
Note: The "long_sample" column is required for the fastqc rule. If your raw fastq files do not have these added characters, you can simply make the long_sample column identical to the sample_name column. This is a quick fix so you don't need to modify the variable calls in downstream rules. 

4. Get required files from MetaPhlAn github (this may need to be updated if databases change)
```bash
wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/sgb_to_gtdb_profile.py -O scripts/sgb_to_gtdb_profile.py
wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/util_fun.py -O scripts/util_fun.py
wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv -O scripts/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv
wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/merge_metaphlan_tables.py -O scripts/merge_metaphlan_tables.py 

# Make sure these scripts have execute permissions
ls -l scripts/sgb_to_gtdb_profile.py 
chmod +x scripts/sgb_to_gtdb_profile.py 
```
5. Load snakemake v9.11.4 into a conda environment (if necessary):
```bash
eval "$(mamba shell hook --shell bash)"
mamba create -n snakemake_env -c conda-forge -c bioconda snakemake=9.11.4
conda activate snakemake_env
```
 6. Install slurm executor plugin for snakemake v8+ (only needs to be done once):
 ```bash
pip install snakemake-executor-plugin-slurm
```
 7. Quick check to make sure there are no errors (dry run):
```bash
snakemake -s scripts/Snakefile --configfile config.yaml -n
```

8. Run the pipeline using one of these methods (meant for using HPC with SLURM scheduler):
* METHOD A - Submit via sbatch script (recommended):
    ```
    sbatch scripts/submit_snakefile.sh
    ````

* METOD B - Run in terminal directly. 
    ```bash
    snakemake -s scripts/Snakefile --configfile config.yaml --executor slurm --jobs 20 --use-conda \
            --default-resources slurm_account=GROUPNAME mem_mb=4096 runtime=600
    ```

9. Monitor progress:
```bash
tail -f logs/snakemake_<jobid>.out
```
 REQUIRED DIRECTORY STRUCTURE:
-----------------------
```
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
# │     ├── envs/
# │        ├── fastqc_multiqc.yaml
# │        ├── bowtie.yaml
# │        ├── fastp.yaml              
# │        └── metaphlan.yaml  
# │
# ├── sample_sheet.txt              [REQUIRED - tab-separated file with sample info]
# │
# └── fastq_files/                  [your input files, demultiplexed]
#     ├── sample1_R1_001.fastq.gz
#     ├── sample1_R2_001.fastq.gz
#     ├── sample2_R1_001.fastq.gz
#     └── sample2_R2_001.fastq.gz
#
```






