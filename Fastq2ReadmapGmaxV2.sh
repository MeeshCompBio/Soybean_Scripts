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

###### FOR MSI USE
#Please modify the module loads commands to the appropritate software numbers
# module load fastqc
# module laod bowtie2
# module load bwa
# module load samtools
# module load bedtools
# module load fastqc
# module load cutadapt


#use getopts for command line arguments
while getopts hf:r::a::m::o: flag; do
    case $flag in
      #this is the help command
        h)
            echo "This is your help information:
            This script was designed for use on MSI by the Stupar Lab.
            Some of these file paths are hard coded for use on MSI but
            that can be easily changed.

            Usage: 
            Fastq2ReadmapGmaxV2.sh -f forward_read.fastq -r <reverse_read.fastq> -a <adapter> -m <use bowtie2 instead of bwa> -o output_directory \n"
            exit 2
            ;;
        f)
            echo "Forward read is: $OPTARG";
            FILE=$OPTARG
            I=1
            A="False"
            ;;
        r)
            echo "Reverse read is: $OPTARG";
            FILE2=$OPTARG
            I=2
            ;;
        a)
            echo "Adapter sequence is: $OPTARG";
            ADAPTER=$OPTARG
            A="True"
            ;;
        m)
            echo "You are opting to use bowtie2 instead of BWA";
            ALIGNER="bowtie2"
            ;;
        o)
            echo "everything will be output to this directory: ${OPTARG%/}"
            OUTPUTDIR="${OPTARG%/}"
            ;;
        ?)
            echo "not a valid argument, try using the -h option for usage information"
            exit 2
            ;;
    esac
done

#Letting the user know if the program will use SE or PE pipeline
if [ $I == 2 ]
then
    echo "setting up for paried-end read pipeline"
elif [ $I == 1 ]
then
    echo "setting up for single-end read pipeline"
else
    echo "error, you might be missing a variable try the -h command for help"
    exit 2
fi




#just get the last part of the file name
filename=$(basename "$FILE")
#get what the extension is
extension="${filename##*.}"
#and now you just have the base name for renaming
basename="${filename%%.*}"

#This command will split the base name and grab the first part of it
#It is currently splitting on "_"
#Ex: "Sample1_R1" ==> "Sample1"
samplename=$(echo ${basename} | cut -f1 -d"_")  


#Quality Checks before running the scripts


#check to see that the forward read exists and is readable
if [ -r $FILE ]
   then
   echo "Forward file exists and it readable"
   else
   echo "This file is not valid, check to see that is a readable fastq file"
   exit
fi

#check same for reverse read if option is flagged
if [ $I == 2 ]
   then
   if [ -r $FILE ]
      then
      echo "Forward file exists and it readable"
      filename2=$(basename "$FILE2")
      #get what the extension is
      extension2="${filename2##*.}"
      #and now you just have the base name for renaming
      basename2="${filename2%%.*}"
      #basename2="${filename2%%.anqrp.fastq.gz}"
   else
      echo "This file is not valid, check to see that is a readable fastq file"
      exit
   fi
fi

if [ $A == "True" ]
    then
    if [ 6 -eq ${#ADAPTER}  ]; then
        echo "Correct Illumina TRUSEQ barcode length"
    else
        echo "Adatper length is the wrong size"
        exit
    fi
fi

if [ -d $OUTPUTDIR  ]; then
   echo "$OUTPUTDIR is a valid directory"
else
   echo "$OUTPUTDIR is not directory"
   exit
fi


#go to Output directory
cd "${OUTPUTDIR}"

#load and run fastqc module for forward and reverse at the same time
if [ $I == 2 ]
   then
      #Run both fastqc programs at once
      fastqc -f fastq -o ${OUTPUTDIR} ${FILE} &
      fastqc -f fastq -o ${OUTPUTDIR} ${FILE2} 
      wait
      echo "Initial FastQC on forward and revese reads worked"
   else
      fastqc -f fastq -o ${OUTPUTDIR} ${FILE}
      echo "Initial FastQC on forward reads worked"
fi

#load and run cutadapt
echo "Starting read and or adapter trimming"


if [ $A == "True" ]
    then
    REVERSEADAPTER1=$(echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}ATCTCGTATGCCGTCTTCTGCTTG" | rev | tr ATGC TACG)
    REVERSEADAPTER2=$(echo AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT | rev | tr ATGC TACG)
    if [ $I == 2 ]
       then
          #load and run cutadapt
          cutadapt \
          -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}ATCTCGTATGCCGTCTTCTGCTTG \
          -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
          -B ${REVERSEADAPTER1} \
          -B ${REVERSEADAPTER2} \
          -f fastq \
          -m 40 \
          -q 30 \
          --quality-base=33 \
          -o ${OUTPUTDIR}${basename}_cutadapt.fastq \
          -p ${OUTPUTDIR}${basename2}_cutadapt.fastq \
          ${FILE} \
          ${FILE2} \
          > ${OUTPUTDIR}${basename}_cutadapt.log

          #run Fastqc to make sure reads were trimmed
          fastqc -f fastq ${OUTPUTDIR}${basename}_cutadapt.fastq &
          fastqc -f fastq ${OUTPUTDIR}${basename2}_cutadapt.fastq 

       else
          #load and run cutadapt for single end reads
          cutadapt \
          -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}ATCTCGTATGCCGTCTTCTGCTTG \
          -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
          -f fastq \
          -m 40 \
          -q 30 \
          --quality-base=33 \
          -o ${OUTPUTDIR}${basename}_cutadapt.fastq \
          ${FILE} \
          > ${OUTPUTDIR}${basename}_cutadapt.log

          fastqc -f fastq ${OUTPUTDIR}${basename}_cutadapt.fastq
    fi
fi
echo "Adapter trimming finished"
#load and run fastx for low complexity sequences


if [ $ALIGNER == "bowtie2" ]
   then
   echo "starting bowtie2 alignment"
   if [ $I == 2 ]
      then
      bowtie2 -N 6 -p 8 \
      --rg-id "@RG\tID:wgs_${samplename}" \
      -x /panfs/roc/groups/13/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0 \
      -1 ${OUTPUTDIR}${basename}_cutadapt.fastq \
      -2 ${OUTPUTDIR}${basename2}_cutadapt.fastq \
      -S ${OUTPUTDIR}bt2${samplename}.sam
      #convert the sam file to a bam file
      samtools view -bSq 20 ${OUTPUTDIR}bt2${samplename}.sam > ${OUTPUTDIR}bt2${samplename}.bam
      #sort and index the bam file
      samtools sort -o ${OUTPUTDIR}bt2${samplename}.sorted -@ 8 ${OUTPUTDIR}bt2${samplename}.bam 
      samtools index ${OUTPUTDIR}bt2${samplename}.sorted.bam
      echo "bowtie2 alignment complete"
      else
      bowtie2 -N 6 -p 8 \
      --rg-id "@RG\tID:wgs_${samplename}" \
      -x /panfs/roc/groups/13/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0 \
      -1 ${OUTPUTDIR}${basename}_cutadapt.fastq \
      -S ${OUTPUTDIR}bt2${samplename}.sam
      #convert the sam file to a bam file
      samtools view -bSq 20 ${OUTPUTDIR}bt2${samplename}.sam > ${OUTPUTDIR}bt2${samplename}.bam
      #sort and index the bam file
      samtools sort -o ${OUTPUTDIR}bt2${samplename}.sorted -@ 8 ${OUTPUTDIR}bt2${samplename}.bam 
      samtools index ${OUTPUTDIR}bt2${samplename}.sorted.bam
      echo "bowtie2 alignment complete"
   fi
   echo"starting bwa alignment"
   else
   if [ $I == 2 ]
      then
      bwa mem -t 8 -w 100 -M -B 6 \
         -R "@RG\tID:wgs_${samplename}\tLB:ES_${samplename}\tSM:WGS_${samplename}\tPL:ILLUMINA" \
         /panfs/roc/groups/13/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0.fa \
         ${OUTPUTDIR}${basename}_cutadapt.fastq \
         ${OUTPUTDIR}${basename2}_cutadapt.fastq  \
         > ${OUTPUTDIR}bwa${samplename}.sam
         echo "BWA alignmentcomplete"
         #convert the sam file to a bam file
         samtools view -bSq 20 ${OUTPUTDIR}bwa${samplename}.sam > ${OUTPUTDIR}bwa${samplename}.bam
         #sort and index the bam file
         samtools sort -o ${OUTPUTDIR}bwa${samplename}.sorted -@ 8 ${OUTPUTDIR}bwa${samplename}.bam 
         samtools index ${OUTPUTDIR}bwa${samplename}.sorted.bam
      else
         bwa mem -t 8 -w 100 -M -B 6 \
         -R "@RG\tID:wgs_${samplename}\tLB:ES_${samplename}\tSM:WGS_${samplename}\tPL:ILLUMINA" \
         /panfs/roc/groups/13/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0.fa \
         ${OUTPUTDIR}${basename}_cutadapt.fastq \
         > ${OUTPUTDIR}bwa${samplename}.sam
         #convert the sam file to a bam file
         samtools view -bSq 20 ${OUTPUTDIR}bwa${samplename}.sam > ${OUTPUTDIR}bwa${samplename}.bam
         #sort and index the bam file
         samtools sort -o ${OUTPUTDIR}bwa${samplename}.sorted -@ 8 ${OUTPUTDIR}bwa${samplename}.bam 
         samtools index ${OUTPUTDIR}bwa${samplename}.sorted.bam
         echo "BWA complete"
   fi
fi



