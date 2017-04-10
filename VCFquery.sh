#ted by Jean-Michel Michno (Meesh)
#email: mich039@umn.edu
# module load cutadapt
set -euo pipefail

#use getopts for command line arguments
while getopts hc:s:t:f::o: flag; do
    case $flag in
      #this is the help command
        h)
            echo "This is your help information:
            This script was designed for use on MSI by the Stupar Lab.
            The VCF file paths are hardcoded to a public directory
            Usage:
            VCFquery.sh -c Chromosome -s start_position -t stop position  -o output_prefix
            Ex:
            VCFquery.sh -c 01 -s 1234 -t 5678  -o myfile
            
            Options:
                    -c Chromosome to use <only need a TWO digit number 01-20>
                    -s start position (number)
                    -t stop position (number)
                    -f which lines to search, use wither 106 or NAM (106 is default)
                    -o output file prefix, what you want to file to be called 
                        without an extension
                    "
            exit 2
            ;;
        c)
            echo " Looking for: Chr$OPTARG";
            CHROM=$OPTARG
            FILE="/panfs/roc/groups/13/stuparr/shared/VariantCalls/106_Genomes_Platypus.vcf"
            HEADER="/panfs/roc/groups/13/stuparr/shared/VariantCalls/106_Genomes_Headers.vcf"
            TYPE="106"
            ;;
        s)
            echo "Start postion is: $OPTARG";
            START=$OPTARG
            ;;
        t)
            echo "Stop postion is: $OPTARG";
            STOP=$OPTARG
            I=2
            ;;
        f)
            echo "Using: $OPTARG samples";
            TYPE=$OPTARG
            ;;
        o)
            echo "The file will be output to: ${OPTARG%/}.vcf"
            OUTPUT="${OPTARG%/}"
            ;;
        ?)
            echo "not a valid argument, try using the -h option for usage information"
            exit 2
            ;;
    esac
done

#just madke a list of available chromosomes
CHR="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20"

#see if the option is not in the list, if so kill script
if echo $CHR | grep -vw $CHROM > /dev/null; then
    echo "error: you must use a number between 01-20 for chromosome (-t)" >&2; exit 1
fi

#make sure that the start and stop options are numbers
NUM='^[0-9]+$'
if ! [[ $START =~ $NUM ]] ; then
   echo "error: you must use a number for start postion (-s)" >&2; exit 1
fi

NUM='^[0-9]+$'
if ! [[ $STOP =~ $NUM ]] ; then
   echo "error: you must use a number for stop postion (-t)" >&2; exit 1
fi

if [ $START -gt $STOP ] ; then
   echo "error: your start position is larger than your stop" >&2; exit 1
fi

if [ $TYPE == "NAM" ]
    then
    FILE="/panfs/roc/groups/13/stuparr/shared/VariantCalls/NAM_Parent_Platypus.vcf"
    HEADER="/panfs/roc/groups/13/stuparr/shared/VariantCalls/NAM_Parent_Headers.vcf"
fi

#grab the header from a pre-parsed file to save time
$(grep "#" $HEADER > $OUTPUT.vcf)
#give a warning so they don't think it hangs
echo "Searching through VCF file, this can take up to:
        5 minutes for the 106
        30 seconds for the NAM"

#grep by chromosome since it is faster then filter based on pos using awk
$(grep -e "Chr$CHROM"  $FILE | awk -v start="$START" -v stop="$STOP" '$2 >= start && $2 <= stop' >> $OUTPUT.vcf)
