#!/bin/bash
set -euo pipefail

#use getopts for command line arguments
while getopts hf::d::t:: flag; do
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
            F="True"
            ;;
        d)
            echo "Changing all file in your specified directory";
            DIR="${OPTARG%/}"
            TEST="False"
            D="True"
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

# if [ $F == "True" ]
#     then
#     #Check to see that it is a readable file 
#     if [ -r $FILE ]
#        then
#        echo "Forward file exists and it readable"
#        else
#        echo "This file is not valid, check to see that is a readable fastq file"
#        exit
#     fi

#     FILENAME=$(basename "$FILE")
#     POSTFIX=$(echo ${FILENAME} | cut -f2- -d"_")
#     NEWNAME=$(head -n1 $FILE | cut -f2 -d"-" | cut -f1 -d":")

#     if [ $TEST == "True" ]
#         then
#         echo "Your file would be renamed: ${NEWNAME}_${POSTFIX} "
#         else
#         mv ${FILE} ${NEWNAME}_${POSTFIX} 
#     fi
# fi


if [ $D == "True" ]
    then
    if [ -d $DIR  ]; 
        then
        echo "$DIR is a valid directory"
        if [ $TEST == "True" ]
            then
            for i in $( ls $DIR);
                do
                    echo $i
                    POSTFIX=$(echo ${DIR}/$i | cut -f2- -d"_")
                    NEWNAME=$(head -n1 ${DIR}/$i | cut -f2 -d"-" | cut -f1 -d":")
                    echo "Your file would be renamed: ${NEWNAME}_${POSTFIX} "
                done
            else
            for i in $( ls $DIR);
                do
                    echo $i
                    POSTFIX=$(echo ${DIR}/$i | cut -f2- -d"_")
                    NEWNAME=$(head -n1 ${DIR}/$i | cut -f2 -d"-" | cut -f1 -d":")
                    mv ${DIR}/${i} ${DIR}/${NEWNAME}_${POSTFIX}
                done
        fi
    else
        echo "$DIR is not directory"
        exit
    fi
fi


#     if [ $TEST == "True" ]
#         then
#         for i in $( ls $DIR);
#             do
#                 echo $i
#                 POSTFIX=$(echo ${DIR}/$i | cut -f2- -d"_")
#                 NEWNAME=$(head -n1 ${DIR}/$i | cut -f2 -d"-" | cut -f1 -d":")
#                 echo "Your file would be renamed: ${NEWNAME}_${POSTFIX} "
#             done
#         else
#         for i in $( ls $DIR);
#             do
#                 echo $i
#                 POSTFIX=$(echo ${DIR}/$i | cut -f2- -d"_")
#                 NEWNAME=$(head -n1 ${DIR}/$i | cut -f2 -d"-" | cut -f1 -d":")
#                 mv ${DIR}/${i} ${DIR}/${NEWNAME}_${POSTFIX}
#             done
#     fi
# fi