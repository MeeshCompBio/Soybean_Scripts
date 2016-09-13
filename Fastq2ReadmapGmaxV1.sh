#!/bin/bash

#This script takes in your forward reads, reverse reads, 6bp adapter, and output directory and returns
##all intermitant files with a sorted bam being the final output
#This script is to version one of the genome.

#How to use the script
#bash Unicon_Fastq2ReadmapGmaxV1.sh \
#/home/stuparr/michnoj0/Scratch/FNPipeline/FNseq**/*  \
#/home/stuparr/michnoj0/Scratch/FNPipeline/FNseq**/*  \
#NNNNNN \
#/home/stuparr/michnoj0/Scratch/FNPipeline/FNseq**/

set -euo pipefail
#use shebang to indicate the path to the interpretor used to execute this script
#set -e terminated the script if any command exits with a non-zero value
#set -u will prevent bash from running if there is an unset variable
#set -o pipefail This setting prevents errors in a pipeline from being masked. 
#If any command in a pipeline fails, that return code will be used as the return code of the whole pipeline.
module load riss_util
module load bwa
module load samtools
module load bedtools
module load fastqc
module load cutadapt
module load fastx_toolkit
module load perl

#first argument will be your first .fastq
FILE="$1"
#first argument will be your second .fastq
FILE2="$2"
#your 6 CAPITALIZED base pair adapter
ADAPTER="$3"
#where you want all of your files to output to
OUTPUTDIR="$4"


#just get the last part of the file name
filename=$(basename "$FILE")
#get what the extension is
extension="${filename##*.}"
#and now you just have the base name for renaming
basename="${filename%%.*}"

filename2=$(basename "$FILE2")
#get what the extension is
extension2="${filename2##*.}"
#and now you just have the base name for renaming
basename2="${filename2%%.*}"
samplename=$(echo ${basename} | cut -f1 -d"_")  


#Quality Checks
if [ -r $FILE ] && [ "$extension" = "fastq" ]; then
   echo "First .fastq file is valid"
else
   echo "This file is not valid, check to see that is a readable fastq file"
   exit
fi

if [ -r $FILE2 ] && [ "$extension2" = "fastq" ]; then
   echo "Second .fastq file is valid"
else
   echo "This file is not valid, check to see that is a readable fastq file"
   exit
fi

if [ 6 -eq ${#ADAPTER}  ]; then
   echo "Correct Adatper length"
else
   echo "Adatper length is the wrong size"
   exit
fi

if [ -d $4  ]; then
   echo "$4 is a valid directory"
else
   echo "$4 is not directory"
   exit
fi

#go to Output directory
 cd "${OUTPUTDIR}"

#load and run fastqc module for forward and reverse at the same time
module load fastqc
fastqc -f fastq ${FILE} &
fastqc -f fastq ${FILE2} 
wait
echo "Initial FastQC worked"
#load and run cutadapt

cutadapt \
-b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC\
${ADAPTER}\
ATCTCGTATGCCGTCTTCTGCTTG \
-b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
-f fastq \
-m 30 \
--quality-base=33 \
-o ${OUTPUTDIR}${basename}_cutadapt.fastq \
${FILE} \
> ${OUTPUTDIR}${basename}_cutadapt.log &

REVERSEADAPTER1=$(echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}ATCTCGTATGCCGTCTTCTGCTTG" | rev | tr ATGC TACG)
REVERSEADAPTER2=$(echo AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT | rev | tr ATGC TACG)

#Reverse reads
cutadapt \
-b ${REVERSEADAPTER1} \
-b ${REVERSEADAPTER2} \
-f fastq \
-m 30 \
--quality-base=33 \
-o ${OUTPUTDIR}${basename2}_cutadapt.fastq \
${FILE2} \
> ${OUTPUTDIR}${basename2}_cutadapt.log 
wait
fastqc -f fastq ${OUTPUTDIR}${basename}_cutadapt.fastq &
fastqc -f fastq ${OUTPUTDIR}${basename2}_cutadapt.fastq 
wait
echo "Cut adapt worked"


#load and run fastx for low complexity sequences
fastx_artifacts_filter -Q33 -v \
-i ${OUTPUTDIR}${basename}_cutadapt.fastq \
-o ${OUTPUTDIR}${basename}_cutadapt_complexity.fastq &

fastx_artifacts_filter -Q33 -v \
-i ${OUTPUTDIR}${basename2}_cutadapt.fastq \
-o ${OUTPUTDIR}${basename2}_cutadapt_complexity.fastq 
wait
fastqc -f fastq ${OUTPUTDIR}${basename}_cutadapt_complexity.fastq &
fastqc -f fastq ${OUTPUTDIR}${basename2}_cutadapt_complexity.fastq 
wait
echo "low complex trim worked"


#trim short sequnces
fastq_quality_trimmer -Q  33 -v -t 20 -l 30 \
-i ${OUTPUTDIR}${basename}_cutadapt_complexity.fastq \
-o ${OUTPUTDIR}${basename}_cutadapt_complexity_qual_trim.fastq &

fastq_quality_trimmer -Q  33 -v -t 20 -l 30 \
-i ${OUTPUTDIR}${basename2}_cutadapt_complexity.fastq \
-o ${OUTPUTDIR}${basename2}_cutadapt_complexity_qual_trim.fastq
wait

fastqc -f fastq ${OUTPUTDIR}${basename}_cutadapt_complexity_qual_trim.fastq &
fastqc -f fastq ${OUTPUTDIR}${basename2}_cutadapt_complexity_qual_trim.fastq
wait
echo "low quality trim worked"


#sync the reads
resync.pl \
-s \
${OUTPUTDIR}${basename}_cutadapt_complexity_qual_trim.fastq  \
${OUTPUTDIR}${basename2}_cutadapt_complexity_qual_trim.fastq


bwa mem -t 8 -w 100 -M -B 6 \
-R "@RG\tID:wgs_${samplename}\tLB:ES_${samplename}\tSM:WGS_${samplename}\tPL:ILLUMINA" \
/home/stuparr/shared/References/Gmax.a1.v1.1/assembly/Gmax_189.fa \
${OUTPUTDIR}${basename}_cutadapt_complexity_qual_trim.fastq.out \
${OUTPUTDIR}${basename2}_cutadapt_complexity_qual_trim.fastq.out  \
>${OUTPUTDIR}bwa${samplename}.sam
echo "BWA complete"


#convert the sam file to a bam file
samtools view -bSq 20 ${OUTPUTDIR}bwa${samplename}.sam > ${OUTPUTDIR}bwa${samplename}.bam

#sort the bam file
samtools sort -o ${OUTPUTDIR}bwa${samplename}.sorted -@ 8 ${OUTPUTDIR}bwa${samplename}.bam 
samtools index ${OUTPUTDIR}bwa${samplename}.sorted.bam



