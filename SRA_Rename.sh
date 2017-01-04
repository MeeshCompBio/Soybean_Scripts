#!/bin/bash
#created by Meesh mich0391@umn.edu

set -euo pipefail

#use getopts for command line arguments
while getopts hf::d::t:: flag; do
    case $flag in
      #this is the help command
        h)
            echo "This is your help information:

            This script will rename your SRRXXX.fastq according to the header of
            your .fastq reads after running fastq-dump with the -F option. You
            can with run this on a single file or on an entire directory of SRR
            files. 

            Usage:
            SRA_Rename.sh -f SRR******.fastq
            or
            SRA_Rename.sh -d <directory>

            Options:
                    -d <directory>: give a directory instead of a file name. NOTE:
                                    This will only look for files starting with SRR
                    -t True       : will print what the new file names would be 
                                    without renaming them
                    -d <directory>: give a directory instead of a file name. NOTE:
                                    This will only loo for files starting with SRR
                    -h            : help command
            "
            exit 2
            ;;
        f)
            echo "File to be modified is $OPTARG";
            FILE=$OPTARG
            TEST="False"
            F="True"
            D="False"
            ;;
        d)
            echo "Changing all file in your specified directory";
            DIR="${OPTARG%/}"
            TEST="False"
            D="True"
            F="False"
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


if [ $F == "True" ]
    then
    if [ -r $FILE ]
        then
        echo "$FILE is readable"
        FILENAME=$(basename "$FILE")
        POSTFIX=$(echo ${FILENAME} | cut -f2- -d"_")
        NEWNAME=$(head -n1 $FILE | cut -f2 -d"-" | cut -f1 -d":")
        if [ $TEST == "True" ]
            then
            echo "Your file would be renamed: ${NEWNAME}_${POSTFIX} "
            else
            mv ${FILE} ${NEWNAME}_${POSTFIX} 
        fi
    else
        echo "$DIR is not directory"
        exit
    fi
fi


#if you used the directory flag
if [ $D == "True" ]
    then
    #check if its actually a directory
    if [ -d $DIR  ]
        then
        echo "$DIR is a valid directory"
        #If you are using the testing parameter
        if [ $TEST == "True" ]
            then
            #loop through every file in the directory
            for i in $( ls ${DIR} | grep ^SRR);
                do
                    echo $i
                    #for each dir split on "_" and take from 2:end
                    POSTFIX=$(echo $i | cut -f2- -d"_")
                    #pull the sequence name from the header
                    NEWNAME=$(head -n1 ${DIR}/$i | cut -f2 -d"-" | cut -f1 -d":")
                    #print what the file name would be to terminal
                    echo "Your file would be renamed: ${NEWNAME}_${POSTFIX} "
                done
            else
            for i in $( ls ${DIR} | grep ^SRR);
                do
                    POSTFIX=$(echo $i | cut -f2- -d"_")
                    NEWNAME=$(head -n1 ${DIR}/$i | cut -f2 -d"-" | cut -f1 -d":")
                    #rewrite file names
                    mv ${DIR}/${i} ${DIR}/${NEWNAME}_${POSTFIX}
                done
        fi
        else
        echo "$DIR is not directory"
        exit
    fi
fi


