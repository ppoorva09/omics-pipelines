#!/bin/bash

# Pipeline in bash to initiate the pre-processing of Bulk RNA-seq data analysis
## Libraries Required:

### fastQC        - Quality control 
### trimmomatic   - Trimming 
### HISAT2        - Alignment 
### featureCounts - Quantification


SECONDS=0


# change Directory
cd /Users/bioinfo_poorva/Desktop/Github/Bulk_RNAseq/

# Step 1.1: fastqc - Quality control
# install fastqc using conda
fastqc data/demo.fastq -o data/

# Step 1.2: trimmomatic to trim reads with poor Quality
# install Trimmomatic by going to the web site and downloading the binary zipfile
# make sure you have java installed
java -jar /Users/bioinfo_poorva/Desktop/Github/Bulk_RNAseq/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 data/demo.fastq data/trimmed.fastq TRAILING:10 -phred33
echo "Trimmomatic finished!"

fastqc data/trimmed.fastq -o data/

# Step 2.1: Run HISAT2
# mkdir HISAT2
# genome index
# wget or curl https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# run Alignment
hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U data/trimmed.fastq | samtools sort -o HISAT2/trimmed.bam
echo "HISAT2 finished!"

# Step 3.1: featureCounts - Quntification
# wget or curl http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
featureCounts -S 2 -a /Users/bioinfo_poorva/Desktop/Github/Bulk_RNAseq/Homo_sapiens.GRCh38.115.gtf -o quants/demo_featurecounts.txt HISAT2/trimmed.bam
echo "featureCounts finished running!"


duration=$SECONDS
echo "Pipeline finished in $duration seconds"


## To run use commands on terminal: 
### chmod 755 Initial_BulkRNAseq_Pipeline.sh 
### ./Initial_BulkRNAseq_Pipeline.sh

