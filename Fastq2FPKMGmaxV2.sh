#!/bin/bash

#This script takes in your forward reads, <optional> reverse reads, <optional> 6bp adapter, and output directory and returns
##all intermitant files with a FPKM list for your genes as a final output
#This script is to version two of the genome.
#currently set up to run on 8 threads

#How to use the script
#bash test.sh -f \
#forward_read.fastq \
#-r <reverse_read.fastq> \
#-a <adapter> \
#-o output_directory


#The parameters used in this script is inteded for soybean based parameters from:
#Hirsch, Cory D., Nathan M. Springer, and Candice N. Hirsch. 
#"Genomic limitations to RNA sequencing expression profiling." 
#The Plant Journal 84.3 (2015): 491-503.
#the only devation to these parameter is a minimum read length of 40bp instead of 50bp


##################IMPORTANT##################
##################IMPORTANT##################
##################IMPORTANT##################
#These are the following programs that need to be installed:
#FastQC, cutadapt, bowtie2, tophat2, cufflinks

#This is what needs to be hard coded in this script
#lines 287, 299, 313 Change this to the path of your .gff or .gtf file
#lines 292, 304 Change this to the prefix path to your genome/index file
#For instance if your genome was ~/home/user/reference.fa it woule be ~/home/user/reference
##################IMPORTANT##################
##################IMPORTANT##################
##################IMPORTANT##################






set -euo pipefail
#use shebang to indicate the path to the interpretor used to execute this script
#set -e terminated the script if any command exits with a non-zero value
#set -u will prevent bash from running if there is an unset variable
#set -o pipefail This setting prevents errors in a pipeline from being masked. 
#If any command in a pipeline fails, that return code will be used as the return code of the whole pipeline.

#use getopts for command line arguments
while getopts hf:r::a::o: flag; do
    case $flag in
      #this is the help command
        h)
            echo "This is your help information \n
Usage: test.sh -f forward_read.fastq -r <reverse_read.fastq> -a <adapter> -o output_directory \n"
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
        o)
            echo "everything will be output to this directory: $OPTARG"
            OUTPUTDIR=$OPTARG
            ;;
        ?)
            echo "not a valid argument, try using the -h option for usage information"
            exit 2
            ;;
    esac
done

#Letting the user know if the program will use SE or PE pipeline
#echo "$FORWARD"
if [ $I == 2 ]
then
    echo "setting up for PE pipeline"
elif [ $I == 1 ]
then
    echo "setting up for SE pipeline"
else
    echo "error, you might be missing a variable"
    exit 2
fi





#just get the last part of the file name
filename=$(basename "$FILE")
#get what the extension is
extension="${filename##*.}"
#and now you just have the base name for renaming
basename="${filename%%.*}"
#this command is for Stacey data
#basename="${filename%%.anqrp.fastq.gz}"


#check to make sure that there is a reverse read
if [ $I == 2 ]
   then
      filename2=$(basename "$FILE2")
      #get what the extension is
      extension2="${filename2##*.}"
      #and now you just have the base name for renaming
      basename2="${filename2%%.*}"
      #basename2="${filename2%%.anqrp.fastq.gz}"
fi

##################IMPORTANT##################
##################IMPORTANT##################
##################IMPORTANT##################
#This line is to modify your sample name so it cuts everything after an "_"
# So File_AAAAA.fastq would then become File
#Comment this out if you dont want to shorted the continueing files names
samplename=$(echo ${basename} | cut -f1 -d"_")
#samplename=${basename} 
##################IMPORTANT##################
##################IMPORTANT##################
##################IMPORTANT##################

#echo ${extension}
#Quality Checks
# if [ -r $FILE ] && [ "$extension" == "fastq" ]; then
#    echo "First .fastq file is valid"
# else
#    echo "First file does not have 'fastq' as an extension, check to see that is a readable fastq file"
#    exit 2
# fi

# if [ $I == 2 ]
#    then
#       if [ -r $FILE2 ] && [ "$extension2" = "fastq" ]; then
#          echo "Second .fastq file is valid"
#       else
#          echo "This file is not valid, check to see that is a readable fastq file"
#          exit 2
#       fi
# fi
if [ $A == "True" ]; 
    then
    if [ 6 -eq ${#ADAPTER}  ]; then
        echo "Correct Illumina adatper length"
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
if [ $I == 2 ];
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

# #these adapters are assuming that you are using the Illumina standard library prep
# echo "Starting read and or adapter trimming"
# REVERSEADAPTER2=$(echo AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT | rev | tr ATGC TACG)

if [ $A == "True" ]
    then
    REVERSEADAPTER1=$(echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}ATCTCGTATGCCGTCTTCTGCTTG" | rev | tr ATGC TACG)
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
          -q 10 \
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
          #load and run cutadapt
          cutadapt \
          -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${ADAPTER}ATCTCGTATGCCGTCTTCTGCTTG \
          -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
          -f fastq \
          -m 40 \
          -q 10 \
          --quality-base=33 \
          -o ${OUTPUTDIR}${basename}_cutadapt.fastq \
          ${FILE} \
          > ${OUTPUTDIR}${basename}_cutadapt.log

          fastqc -f fastq ${OUTPUTDIR}${basename}_cutadapt.fastq
    fi
    else
      if [ $I == 2 ]
        then
          #load and run cutadapt just searching for standard illumina adapters
          cutadapt \
          -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
          -B ${REVERSEADAPTER2} \
          -f fastq \
          -m 40 \
          -q 10 \
          --quality-base=33 \
          -o ${OUTPUTDIR}${basename}_cutadapt.fastq \
          -p ${OUTPUTDIR}${basename2}_cutadapt.fastq \
          ${FILE} \
          ${FILE2} \
          > ${OUTPUTDIR}${basename}_cutadapt.log

          #run Fastqc to make sure reads were trimmed
          fastqc -f fastq ${OUTPUTDIR}${basename}_cutadapt.fastq &
          fastqc -f fastq ${OUTPUTDIR}${basename2}_cutadapt.fastq 
          wait
          echo "Cut adapt worked"

        else
          #load and run cutadapt just searching for standard illumina adapters
          cutadapt \
          -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
          -f fastq \
          -m 40 \
          -q 10 \
          --quality-base=33 \
          -o ${OUTPUTDIR}${basename}_cutadapt.fastq \
          ${FILE} \
          > ${OUTPUTDIR}${basename}_cutadapt.log

          fastqc -f fastq ${OUTPUTDIR}${basename}_cutadapt.fastq
      fi
fi

# ###To build the transcriptome index
# tophat \
# #GFF file
# -G /home/misc00/jmichno/MeeshCM/Data/References/Wm82.a2.v1/Gmax_275_Wm82.a2.v1.gene_exons.gff3 \
# #what you want your transcriptome index directory to be called with the prefix
# --transcriptome-index=transcriptome_data/Gmax_275_v2.0 \
# #number of threads
# -p 8 \
# #where your bowtie indexes are kept with prefix
# /home/misc00/jmichno/MeeshCM/Data/References/Wm82.a2.v1/Gmax_275_v2.0



if [ $I == 2 ]
  then
      #run tophat2 to generate the bam file
      tophat2 \
      -G /home/misc00/jmichno/MeeshCM/Data/References/Wm82.a2.v1/Gmax_275_Wm82.a2.v1.gene_exons.gff3 \
      -p 7 \
      -i 5 \
      -I 20000 \
      -o ${OUTPUTDIR}${samplename} \
      /home/misc00/jmichno/MeeshCM/Data/References/Wm82.a2.v1/Gmax_275_v2.0 \
      ${OUTPUTDIR}${basename}_cutadapt.fastq \
      ${OUTPUTDIR}${basename2}_cutadapt.fastq

    else
      #run tophat2 to generate the bam file
      tophat2 \
      -G /home/misc00/jmichno/MeeshCM/Data/References/Wm82.a2.v1/Gmax_275_Wm82.a2.v1.gene_exons.gff3 \
      -p 7 \
      -i 5 \
      -I 20000 \
      -o ${OUTPUTDIR}${samplename} \
      /home/misc00/jmichno/MeeshCM/Data/References/Wm82.a2.v1/Gmax_275_v2.0 \
      ${OUTPUTDIR}${basename}_cutadapt.fastq

fi

cd ${OUTPUTDIR}${samplename}

#use cufflinks to generate FPKM values
cufflinks \
-G /home/misc00/jmichno/MeeshCM/Data/References/Wm82.a2.v1/Gmax_275_Wm82.a2.v1.gene_exons.gff3 \
-I 20000 \
--min-intron-length 5 \
${OUTPUTDIR}${samplename}/accepted_hits.bam

#rename the files with the sample name in front
mv genes.fpkm_tracking "${samplename}_genes.fpkm_tracking"
mv isoforms.fpkm_tracking "${samplename}_isoforms.fpkm_tracking"
mv skipped.gtf "${samplename}_skipped.gtf"












