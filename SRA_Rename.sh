#!/bin/bash
set -euo pipefail

#use getopts for command line arguments
while getopts hf:d::t:: flag; do
    case $flag in
      #this is the help command
        h)
            echo "This is your help information:

            This script will rename your file.sra according to the header of
            your .fastq reads after running fastq-dump with the -F option

            Usage:
            SRA_Rename.sh -f Filename.sra
            "
            exit 2
            ;;
        f)
            echo "File to be modified is $OPTARG";
            FILE=$OPTARG
            TEST="False"
            DIR="False"
            ;;
        d)
            echo "running a test instance";
            DIR="True"
            ;;
        t)
            echo "running a test instance";
            TEST="True"
            ;;
        ?)
            echo "not a valid argument, try using the -h option for usage information"
            exit 2
            ;;
    esac
done

#Check to see that it is a readable file 
if [ -r $FILE ]
   then
   echo "Forward file exists and it readable"
   else
   echo "This file is not valid, check to see that is a readable fastq file"
   exit
fi

FILENAME=$(basename "$FILE1")
POSTFIX=$(echo ${filename} | cut -f2- -d"_")
NEWNAME=$(head -n1 $FILE | cut -f2 -d"-" | cut -f1 -d":")

if [ $TEST == "True" ]
    then
    echo "$NEWNAME"
    else
    mv ${FILE} ${NEWNAME}${POSTFIX} 
fi