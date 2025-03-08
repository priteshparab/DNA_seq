#!/bin/bash

# STEP 1 : DOWNLOADING THE DATA
#SRA Database : https://www.ncbi.nlm.nih.gov/sra/ERX2405312%5baccn%5d

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR235/009/ERR2356709/ERR2356709_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR235/009/ERR2356709/ERR2356709_2.fastq.gz


#Downloading reference genome
#UCSC BROWSER : https://genome.ucsc.edu/

#COMMAND:
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz

#Uncompression of the data
gunzip ERR2356709_1.fastq.gz
gunzip ERR2356709_2.fastq.gz
gunzip chr17.fa.gz



#STEP 2 : QUALITY CHECK OF FORWARD AND REVERSE READS

#To install FastQC
#COMMAND : 
sudo apt-get install fastqc

#COMMAND:
fastqc ERR2356709_1.fastq
fastqc ERR2356709_2.fastq





#STEP 2: READ TRIMMING 
#To install fastp:
#COMMAND:
sudo apt-get install fastp

#COMMAND:
fastp -i ERR2356709_1.fastq -o fastp_ERR2356709_1.fastq -I ERR2356709_2.fastq -O fastp_ERR2356709_2.fastq --adapter_fasta adapter.fasta

#NOTE : adapter.fasta will contain all the adapter or overrepresented sequences to be trimmed in fasta format




#STEP 4 : AGAIN CHECKING THEIR QUALITY
#COMMAND:
fastqc fastp_ERR2356709_1.fastq
fastqc fastp_ERR2356709_2.fastq




#STEP 5: READ ALIGNMENT TO REFERENCE GENOME & INDEXING
# A) Indexing of reference genome using BWA
#Install BWA:
sudo apt-get install bwa

#COMMAND:
bwa index -a bwtsw chr17.fa

# B) Aligning reads on indexed genome
#Command:
bwa mem -t 2 chr17.fa fastp_ERR2356709_1.fastq fastp_ERR2356709_2.fastq > bwa_ERR2356709.bam




#STEP 6: CONVERT BAM TO SAM AND SORTING.
#BAM TO SAM
#COMMAND:
samtools view bwa_ERR2356709.bam > bwa_ERR2356709.sam

#SORTING
#COMMAND:
samtools sort bwa_ERR2356709.bam > sorted_ERR2356709.bam




#STEP 7: MARK DUPLICATES 
#SAMtools rmdup (deprecated)-Still can be used
#COMMAND:
samtools rmdup -sS sorted_ERR2356709.bam rmdup_ERR2356709.bam

#Re-index BAM File
samtools index rmdup_ERR2356709.bam




#STEP 8: VARIANT CALLING
#To Download GATK : 
wget -c https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip

#Unzip the file-
unzip gatk-4.3.0.0.zip

#To install picard tools-
sudo apt-get install picard-tools

#To create sequence dictionary of the reference genome:
picard-tools CreateSequenceDictionary R=chr17.fa O=chr17.dict

#Generating .fai file from the reference:
samtools faidx chr17.fa

#Preparing input BAM file for GATK in picard tool format by adding read group information and sorting the file.

#COMMAND:
picard-tools AddOrReplaceReadGroups I=rmdup_ERR2356709.bam O=picard_output.bam RGLB=lib1 RGPL=illumina RGPU=run RGSM=ERR2356709 SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

#Call variations using GATK
#COMMAND:
#Copy the path of the jar file and paste in the command:
java -jar /mnt/c/Users/pritesh/Desktop/dna_seq/tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar HaplotypeCaller -R chr17.fa -I picard_output.bam -O GATK_output.vcf




#STEP 9: VARIANT FILTERING 
#COMMAND
bcftools filter -i 'INFO/DP>10 && QUAL>30' -o filtered_GATK_output.vcf GATK_output.vcf


#STEP 10 : FUNCTIONAL ANNOTATION USING ANNOVAR (optional)
#DOWNLOADING ANNOVAR TOOL
wget -c http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz

#EXTRACT 
tar -xvzf annovar.latest.tar.gz

#MAKE A FOLDER WHERE ANNOVAR OUTPUT WILL BE STORED
mkdir output

##TRANSFER THE GATK_VCF FILE TO THE TOOL FOLDER
cp GATK_output.vcf tools/ 

#NOW TO ANNOTATE SO USE THIS COMMAND:
annovar/table_annovar.pl GATK_output.vcf annovar/humandb/ -buildver hg19 -out output/myannovar --thread 4 -remove -protocol refGene -operation g -nastring . -vcfinput -polish

#EXPLORE THE TXT FILE FOR FUNCTIONAL ANNOTATIONS
myannovar.hg19_multianno.txt




#STEP 11 : VARIATION VISUALIZATION USING IGV
#MAKE A FOLDER WHERE REQUIRED FILES FOR IGV WILL BE STORED
mkdir IGV

#TRANFER OF FILES
picard_output.bam
chr17.fa
GATK_output.vcf

#INDEXING FILES -.bam AND .fa
samtools index picard_output.bam
samtools faidx chr17.fa

echo "Pipeline execution completed. Load chr17.fa, picard_output.bam, and GATK_output.vcf in IGV for visualization."




