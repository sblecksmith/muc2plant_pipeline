# run this in the terminal to make a sample sheet
# Create the sample sheet from your existing files
 echo -e "sample_name\tlong_sample\tr1_path\tr2_path" > sample_sheet.txt

# Add all your samples (this will extract the sample name from the filename)
# The pattern to match (_S[0-9]*_L[0-9]*_R1_001.fastq.gz) may need to change for your file names
 for r1 in fastq_files/*_R1_001.fastq.gz; do
     r2="${r1/_R1/_R2}"
     sample=$(basename "$r1" | sed 's/_S[0-9]*_L[0-9]*_R1_001.fastq.gz//')
     long_sample=$(basename "$r1" | 's/_R1_001.fastq.gz//')
     echo -e "${sample}\t${long_sample}\t${r1}\t${r2}" >> sample_sheet.txt
 done
 
# Check the first few lines
head sample_sheet.txt
