#!/bin/bash

#This script takes in your forward reads, reverse reads, 6bp adapter, and output directory and returns
##all intermitant files with a sorted bam being the final output
#This script is to version two of the genome.


set -euo pipefail
#use shebang to indicate the path to the interpretor used to execute this script
#set -e terminated the script if any command exits with a non-zero value
#set -u will prevent bash from running if there is an unset variable
#set -o pipefail This setting prevents errors in a pipeline from being masked. 
#If any command in a pipeline fails, that return code will be used as the return code of the whole pipeline.

###### FOR MSI USE
#Please modify the module loads commands to the appropritate software versions
# module load fastqc
# module load bowtie2
# module load bwa
# module load samtools/1.5
# module load bedtools
# module load fastqc
# module load cutadapt
# module load java
# SILENT_JAVA_OPTIONS="$_JAVA_OPTIONS"
# unset _JAVA_OPTIONS

#use getopts for command line arguments
while getopts hf:r::a::m::t::u::o: flag; do
    case $flag in
      #this is the help command
        h)
            echo "This is your help information:

            This script was designed for use on MSI by the Stupar Lab.
            Some of these file paths are hard coded for use on MSI but
            that can be easily changed.

            Usage:
            Fastq2ReadmapGmaxV2.sh -f forward_read.fastq -r <reverse_read.fastq> -o output_directory 
            
            Options:
                    -a <adapter> :6bp Illumina TruSeq barcode
                    -m True :If flag is used, then bowtie will be used (default bwa)
                    -t INT :number of threads (default 2)
                    -u Set the location to your .fa reference file
                    "
            exit 2
            ;;
        f)
            echo "Forward read is: $OPTARG";
            FILE=$OPTARG
            I=1
            A="False"
            ALIGNER="bwa"
            THREADS=2
            REFERENCE="/panfs/roc/groups/13/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0.fa"
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
        t)
            echo "$OPTARG threads will be used";
            THREADS=$OPTARG
            ;;
        u)
            echo "$OPTARG is your reference .fa location";
            REFERENCE=$OPTARG
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

if [ -r $REFERENCE ]
   then
   echo "Reference fasta file is readable"
   else
   echo "This file is not valid, check to see that is a readable fasta file"
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
        echo "Adatper length is the wrong size <6bp Illumina TruSeq barcode>"
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

#output everything into a log file for easy error identification
echo "Your log file will be located in ${OUTPUTDIR}/FQ2RM_${samplename}.log"
exec > >(tee "${OUTPUTDIR}/FQ2RM_${samplename}.log") 2>&1


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
echo "Starting adapter trimming"
REVERSEADAPTER2=$(echo AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT | rev | tr ATGC TACG)
if [ $A == "True" ]
    then
    REVERSEADAPTER1=$(echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}ATCTCGTATGCCGTCTTCTGCTTG" | rev | tr ATGC TACG)
    if [ $I == 2 ]
       then
          #load and run cutadapt for paired end with adaptor
          cutadapt \
          -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}ATCTCGTATGCCGTCTTCTGCTTG \
          -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
          -B ${REVERSEADAPTER1} \
          -B ${REVERSEADAPTER2} \
          -f fastq \
          -m 40 \
          -q 20 \
          --quality-base=33 \
          -o ${OUTPUTDIR}/${basename}_cutadapt.fastq \
          -p ${OUTPUTDIR}/${basename2}_cutadapt.fastq \
          ${FILE} \
          ${FILE2} \
          > ${OUTPUTDIR}/${basename}_cutadapt.log

          #run Fastqc to make sure reads were trimmed
          fastqc -f fastq ${OUTPUTDIR}/${basename}_cutadapt.fastq &
          fastqc -f fastq ${OUTPUTDIR}/${basename2}_cutadapt.fastq 

       else
          #load and run cutadapt for single end reads with 6bp adaptor
          cutadapt \
          -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}ATCTCGTATGCCGTCTTCTGCTTG \
          -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
          -f fastq \
          -m 40 \
          -q 20 \
          --quality-base=33 \
          -o ${OUTPUTDIR}/${basename}_cutadapt.fastq \
          ${FILE} \
          > ${OUTPUTDIR}/${basename}_cutadapt.log

          fastqc -f fastq ${OUTPUTDIR}/${basename}_cutadapt.fastq
    fi
fi

#if no adapter is applied then just run cutadadapt with standard illumina
#to trim based on quiality scores regardless
if [ $A == "False" ]
  then
  echo "Trimming based on quality scores"
      if [ $I == 2 ]
        then
          #load and run cutadapt just searching for standard illumina adapters
          cutadapt \
          -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
          -B ${REVERSEADAPTER2} \
          -f fastq \
          -m 40 \
          -q 20 \
          --quality-base=33 \
          -o ${OUTPUTDIR}/${basename}_cutadapt.fastq \
          -p ${OUTPUTDIR}/${basename2}_cutadapt.fastq \
          ${FILE} \
          ${FILE2} \
          > ${OUTPUTDIR}/${basename}_cutadapt.log

          #run Fastqc to make sure reads were trimmed
          fastqc -f fastq ${OUTPUTDIR}/${basename}_cutadapt.fastq &
          fastqc -f fastq ${OUTPUTDIR}/${basename2}_cutadapt.fastq 
          wait
        else
          #load and run cutadapt just searching for standard illumina adapters
          cutadapt \
          -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
          -f fastq \
          -m 40 \
          -q 20 \
          --quality-base=33 \
          -o ${OUTPUTDIR}/${basename}_cutadapt.fastq \
          ${FILE} \
          > ${OUTPUTDIR}/${basename}_cutadapt.log

          fastqc -f fastq ${OUTPUTDIR}/${basename}_cutadapt.fastq
      fi
fi

echo "Adapter trimming finished"


#this is the default aligner
if [ $ALIGNER == "bwa" ]
  then
  echo "starting bwa alignment"
   if [ $I == 2 ]
      then
      bwa mem -t ${THREADS} -k 8 -r 1.0 -M -T 85 \
         -R "@RG\tID:wgs_${samplename}\tLB:ES_${samplename}\tSM:WGS_${samplename}\tPL:ILLUMINA" \
         ${REFERENCE} \
         ${OUTPUTDIR}/${basename}_cutadapt.fastq \
         ${OUTPUTDIR}/${basename2}_cutadapt.fastq  \
         > ${OUTPUTDIR}/bwa${samplename}.sam
         echo "BWA alignmentcomplete"
         #convert the sam file to a bam file
         samtools view -@ ${THREADS} -bSq 20 ${OUTPUTDIR}/bwa${samplename}.sam > ${OUTPUTDIR}/bwa${samplename}.bam
         #we made the bam and no longer need the sam
         rm ${OUTPUTDIR}/bwa${samplename}.sam
         #sort and index the bam file
         samtools sort -@ ${THREADS} -m 800M -T ${samplename} -o ${OUTPUTDIR}/bwa${samplename}.sorted.bam ${OUTPUTDIR}/bwa${samplename}.bam 
         #don't need the original bam since he now have a sorted version of it
         rm ${OUTPUTDIR}/bwa${samplename}.bam
         samtools index -@ ${THREADS} ${OUTPUTDIR}/bwa${samplename}.sorted.bam
         #Mark duplicate reads
         #java -jar -Xmx12g \
         #/panfs/roc/groups/13/stuparr/mich0391/Software/picard/build/libs/picard.jar \
         #MarkDuplicates I=bwa${samplename}.sorted.bam \
         #O=bwa${samplename}_sorted_dedupped.bam \
         #CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
         #M=bwa${samplename}.output.metrics
      else
         bwa mem -t ${THREADS} -k 8 -r 1.0 -M -T 85 \
         -R "@RG\tID:wgs_${samplename}\tLB:ES_${samplename}\tSM:WGS_${samplename}\tPL:ILLUMINA" \
         ${REFERENCE} \
         ${OUTPUTDIR}/${basename}_cutadapt.fastq \
         > ${OUTPUTDIR}/bwa${samplename}.sam
         #convert the sam file to a bam file
         samtools view -@ ${THREADS} -bSq 20 ${OUTPUTDIR}/bwa${samplename}.sam > ${OUTPUTDIR}/bwa${samplename}.bam
         #we made the bam and no longer need the sam
         rm ${OUTPUTDIR}/bwa${samplename}.sam
         #sort and index the bam file
         samtools sort -@ ${THREADS} -m 800M -T ${samplename} -o ${OUTPUTDIR}/bwa${samplename}.sorted.bam ${OUTPUTDIR}/bwa${samplename}.bam
         rm ${OUTPUTDIR}/bwa${samplename}.bam
         samtools index -@ ${THREADS} ${OUTPUTDIR}/bwa${samplename}.sorted.bam
         #Mark duplicate reads
         #java -jar -Xmx12g \
         #/panfs/roc/groups/13/stuparr/mich0391/Software/picard/build/libs/picard.jar \
         #MarkDuplicates I=bwa${samplename}.sorted.bam \
         #O=bwa${samplename}_sorted_dedupped.bam \
         #CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
         #M=bwa${samplename}.output.metrics
         echo "BWA complete"
   fi
fi


#check if the user is using bowtie
if [ $ALIGNER == "bowtie2" ]
   then
   echo "starting bowtie2 alignment"
   if [ $I == 2 ]
      then
      bowtie2 -p ${THREADS} \
      --rg-id "@RG\tID:wgs_${samplename}" \
      -x /panfs/roc/groups/13/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0 \
      -1 ${OUTPUTDIR}/${basename}_cutadapt.fastq \
      -2 ${OUTPUTDIR}/${basename2}_cutadapt.fastq \
      -S ${OUTPUTDIR}/bt2${samplename}.sam
      #convert the sam file to a bam file
      samtools view -bSq 20 ${OUTPUTDIR}/bt2${samplename}.sam > ${OUTPUTDIR}/bt2${samplename}.bam
      #we made the bam and no longer need the sam
      rm ${OUTPUTDIR}/bwa${samplename}.sam
      #sort and index the bam file
      samtools sort -@ ${THREADS} -m 800M ${OUTPUTDIR}/bt2${samplename}.bam ${OUTPUTDIR}/bt2${samplename}.sorted
      samtools index ${OUTPUTDIR}/bt2${samplename}.sorted.bam
      echo "bowtie2 alignment complete"
      else
      bowtie2 -p ${THREADS} \
      --rg-id "@RG\tID:wgs_${samplename}" \
      -x /panfs/roc/groups/13/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0 \
      -1 ${OUTPUTDIR}/${basename}_cutadapt.fastq \
      -S ${OUTPUTDIR}/bt2${samplename}.sam
      #convert the sam file to a bam file
      samtools view -bSq 20 ${OUTPUTDIR}/bt2${samplename}.sam > ${OUTPUTDIR}/bt2${samplename}.bam
      #we made the bam and no longer need the sam
      rm ${OUTPUTDIR}/bwa${samplename}.sam
      #sort and index the bam file
      samtools sort -@ ${THREADS} -m 800M ${OUTPUTDIR}/bt2${samplename}.bam ${OUTPUTDIR}/bt2${samplename}.sorted
      samtools index ${OUTPUTDIR}/bt2${samplename}.sorted.bam
      echo "bowtie2 alignment complete"
   fi
fi




