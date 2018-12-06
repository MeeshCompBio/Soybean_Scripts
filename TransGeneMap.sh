#!/bin/bash

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
# module load samtools/1.6
# module load bedtools
# module load fastqc
# module load trimmomatic
# module load java
# SILENT_JAVA_OPTIONS="$_JAVA_OPTIONS"
# unset _JAVA_OPTIONS

#use getopts for command line arguments
while getopts hf:r:a::t::u::s:v:i::g::o: flag; do
    case $flag in
      # this is the help command
        h)
            echo "This is your help information:

            This script was designed for use on MSI by the Stupar Lab to 
            detect soybean trangenic insertion events. You must use paired 
            end reads.

            Usage:
            TransGeneMap.sh. -f forward_read.fastq -r reverse_read.fastq -s Path_to_trimmomatic_DIR -o output_directory 
            
            Version: 1.0

            Options:
                    -a Trimmomatic adapter file name:
                        must be located in Trimmomatic/adapter directory.
                    -t INT :number of threads (default 2).
                    -u Set the location to your reference genome (not vector)
                    -s Path to trimmomatic DIR.
                    -v Path to vector fasta file
                    -i True :If flag is set to 'True', then bwa and bowtie2
                             will index your vector.fa file
                    -g True :If flag is set to 'True', then bwa will also 
                             run a regular genome alignment.
                    -o output DIR.
                    "
            exit 2
            ;;
        f)
            echo "Forward read is: $OPTARG";
            FILE=$OPTARG
            I=1
            A="False"
            THREADS=2
            # hard coded paths are for me, will be overwritten
                # if other users use the corresponding flag.
            REFERENCE="/panfs/roc/groups/13/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0.fa"
            TRIMMOMATIC="/panfs/roc/groups/13/stuparr/mich0391/Software/Trimmomatic-0.36"
            INDEX="False"
            GENOMEALIGN="False"
            ;;
        r)
            echo "Reverse read is: $OPTARG";
            FILE2=$OPTARG
            I=2
            ;;
        a)
            echo "Adapter filename is: $OPTARG";
            ADAPTER=$OPTARG
            A="True"
            ;;
        t)
            echo "$OPTARG threads will be used";
            THREADS=$OPTARG
            ;;
        u)
            echo "$OPTARG is your species reference";
            REFERENCE=$OPTARG
            ;;
        s)
            echo "$OPTARG is the DIR with you trimmomatic.jar file";
            TRIMMOMATIC="${OPTARG%/}"
            ;;
        v)
            echo "$OPTARG is your vector file (must be .fa)";
            VECTOR=$OPTARG
            ;;
        i)
            echo "Your vector file will be indexed";
            INDEX="True"
            ;;
        g)
            echo "Your vector file will be indexed";
            GENOMEALIGN="True"
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




#just get the last part of the file name
filename=$(basename "$FILE")
#get what the extension is
extension="${filename##*.}"
#and now you just have the base name for renaming
basename="${filename%%.*}"

#get this the filepath of the reference for indexing
REFDIR=$(dirname "${REFERENCE}")
REFNAME=$(basename $REFERENCE)
#Remove the .fa to get the prefix for indexing
REFPREFIX=$(echo $REFNAME | awk -F'.fa' '{print $1}')


#get this the filepath of the reference for indexing
TGDIR=$(dirname "${VECTOR}")
TRANGENEMAP=$(basename $VECTOR)
#Remove the .fa to get the prfix for indexing
TGPREFIX=$(echo $TRANGENEMAP | awk -F'.fa' '{print $1}')

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
   echo "$FILE is not valid, check to see that is a readable fastq file"
   exit
fi

if [ -r $FILE2 ]
    then
    echo "Reverse file exists and it readable"
    filename2=$(basename "$FILE2")
    #get what the extension is
    extension2="${filename2##*.}"
    #and now you just have the base name for renaming
    basename2="${filename2%%.*}"
    #basename2="${filename2%%.anqrp.fastq.gz}"
else
    echo "$FILE2 is not valid, check to see that is a readable fastq file"
    exit
fi

if [ -r $REFERENCE ]
    then
    echo "Reference .fa file is readable"
else
    echo "Reference .fa file is not readable"
    exit
fi


# check to see that trimmomatic DIR has a .jar file
if [ -r ${TRIMMOMATIC}/trimmomatic-0.36.jar ]
   then
   echo "Trimmomatic jar file is present"
else
   echo "$TRIMMOMATIC is not valid, check to see that the specified directory has the .jar file"
   exit
fi

# check to see that the vector file is readable
if [ -r $VECTOR ]
   then
   echo "Vector file exists"
else
   echo "$VECTOR is not valid, check to see that is a readable GFF"
   exit
fi



#Chech trimmomatic adapter file
if [ $A == "True" ]
    then
    if [ -r ${TRIMMOMATIC}/adapters/${ADAPTER} ]; then
        echo "Trimmomatic adaptor file is readable"
else
        echo "Error: make sure adaptor file in in trimomatic/adapter DIR"
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
echo "Your log file will be located in ${OUTPUTDIR}/TGRM_${samplename}.log"
exec > >(tee "${OUTPUTDIR}/TGRM_${samplename}.log") 2>&1


#load and run fastqc module for forward and reverse at the same time

#Run both fastqc programs at once
parallel --gnu -j 2 \
fastqc -f fastq -o ${OUTPUTDIR} \
{1} \
::: ${FILE} ${FILE2}

# some servers have preset options with don't play nice when using this 
# script in parallel
unset _JAVA_OPTIONS

#load and run trimmomatic
echo "Starting adapter trimming"
if [ $INDEX == "True" ]
    then
    echo "Starting index"
    bwa index -p ${TGDIR}/${TGPREFIX} ${VECTOR}
    bowtie2-build ${VECTOR} ${TGDIR}/${TGPREFIX}
fi

#if no adapter is applied then just run cutadadapt with standard illumina
#to trim based on quiality scores regardless
if [ $A == "False" ]
    then
    echo "Trimming based on quality scores"
    #load and run trimmomatic just searching for standard illumina adapters
    java -Xmx${THREADS}G -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar \
    PE \
    -threads ${THREADS} \
    ${FILE} \
    ${FILE2} \
    ${OUTPUTDIR}/${basename}_trimmomatic.fastq \
    ${OUTPUTDIR}/${basename}_trimmomatic_unpaired.fastq \
    ${OUTPUTDIR}/${basename2}_trimmomatic.fastq \
    ${OUTPUTDIR}/${basename2}_trimmomatic_unpaired.fastq \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:40

    #run Fastqc to make sure reads were trimmed
    parallel --gnu -j 2 \
    fastqc -f fastq -o ${OUTPUTDIR} \
    {1} \
    ::: ${OUTPUTDIR}/${basename}_trimmomatic.fastq ${OUTPUTDIR}/${basename2}_trimmomatic.fastq
fi

echo "Adapter trimming finished"


if [ $A == "True" ]
    then
    #load and run trimmomatic for paired end with adaptor
    java -Xmx${THREADS}G -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar \
    PE \
    -threads ${THREADS} \
    ${FILE} \
    ${FILE2} \
    ${OUTPUTDIR}/${basename}_trimmomatic.fastq \
    ${OUTPUTDIR}/${basename}_trimmomatic_unpaired.fastq \
    ${OUTPUTDIR}/${basename2}_trimmomatic.fastq \
    ${OUTPUTDIR}/${basename2}_trimmomatic_unpaired.fastq \
    ILLUMINACLIP:${TRIMMOMATIC}/adapters/${ADAPTER}:2:30:10:4:TRUE \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:40
    #run Fastqc to make sure reads were trimmed
    parallel --gnu -j 2 \
    fastqc -f fastq -o ${OUTPUTDIR} \
    {1} \
    ::: ${OUTPUTDIR}/${basename}_trimmomatic.fastq ${OUTPUTDIR}/${basename2}_trimmomatic.fastq
fi


#this is the default aligner
if [ $GENOMEALIGN == "True" ]
    then
    echo "starting bwa alignment with trimmed reads to genome"
    bwa mem -t ${THREADS} -k 8 -r 1.0 -M -T 85 \
    -R "@RG\tID:wgs_${samplename}\tLB:ES_${samplename}\tSM:WGS_${samplename}\tPL:ILLUMINA" \
    ${REFERENCE} \
    ${OUTPUTDIR}/${basename}_trimmomatic.fastq \
    ${OUTPUTDIR}/${basename2}_trimmomatic.fastq \
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
fi


echo "Mapping reads to transgene"
echo ${VECTOR}
bwa mem -t ${THREADS} \
-R "@RG\tID:wgs_TG${samplename}\tLB:ES_TG${samplename}\tSM:WGS_TG${samplename}\tPL:ILLUMINA" \
${TGDIR}/${TGPREFIX} \
${OUTPUTDIR}/${basename}_trimmomatic.fastq \
${OUTPUTDIR}/${basename2}_trimmomatic.fastq \
> ${OUTPUTDIR}/bwa_TGmap_${samplename}.sam

# Convert sam to bam
samtools view -bhS -@ ${THREADS} ${OUTPUTDIR}/bwa_TGmap_${samplename}.sam > ${OUTPUTDIR}/bwa_TGmap_${samplename}.bam
rm ${OUTPUTDIR}/bwa_TGmap_${samplename}.sam
#sort and compress the sam to a bam
samtools sort -@ ${THREADS} -m 800M -T TGmap_${samplename} -o ${OUTPUTDIR}/bwa_TGmap_${samplename}.sorted.bam ${OUTPUTDIR}/bwa_TGmap_${samplename}.bam
#don't need the original bam since he now have a sorted version of it
rm ${OUTPUTDIR}/bwa_TGmap_${samplename}.bam
#index bam
samtools index -@ ${THREADS} ${OUTPUTDIR}/bwa_TGmap_${samplename}.sorted.bam


echo "Pulling out orphaned reads"
samtools view ${OUTPUTDIR}/bwa_TGmap_${samplename}.sorted.bam \
 | Orphan.pl \
 -gs * \
 -o TGmap_${samplename}_orphan

echo "Mapping orphaned reads to genome ref"
bowtie2 \
--local \
--very-sensitive-local \
-x ${REFDIR}/${REFPREFIX} \
-U TGmap_${samplename}_orphan.fastq \
-S TGmap_BacktoRef_${samplename}_orphan.sam

echo "Pipeline complete"
