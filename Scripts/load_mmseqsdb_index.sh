#!/bin/bash

helpFunction()
{
   echo "Usage: load_mmseqsdb_index.sh -n EMUDB"
   echo -e "\t-n <str> Database Name."
   exit 1 # Exit script after printing help
}

while getopts "n:" opt
do
    case "$opt" in
    n )
        db="$OPTARG"
        ;;
    ? ) helpFunction 
        ;;
    esac
done

if [ -z "$db" ]
    then
    echo "Please provide the database name.";
    helpFunction
fi

if [ "$db" == "REFSEQ" ]; then

    MMSEQS_IDX=$(grep MMSEQS_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export MMSEQS_REFSEQ="//;s/"//g;s/$/\/REFSEQ_MMSEQS.idx/')

elif [ "$db" == "GTDB" ]; then
    
    MMSEQS_IDX=$(grep MMSEQS_GTDB ~/.bashrc | tail -n 1 | sed 's/export MMSEQS_GTDB="//;s/"//g;s/$/\/GTDB_MMSEQS.idx/')

elif [ "$db" == "MIMT" ]; then

    MMSEQS_IDX=$(grep MMSEQS_MIMT ~/.bashrc | tail -n 1 | sed 's/export MMSEQS_MIMT="//;s/"//g;s/$/\/MIMT_MMSEQS.idx/')

elif [ "$db" == "GSR" ]; then

    MMSEQS_IDX=$(grep MMSEQS_GSR ~/.bashrc | tail -n 1 | sed 's/export MMSEQS_GSR="//;s/"//g;s/$/\/GSR_MMSEQS.idx/')

elif [ "$db" == "EMUDB" ]; then

    MMSEQS_IDX=$(grep MMSEQS_EMU ~/.bashrc | tail -n 1 | sed 's/export MMSEQS_EMU="//;s/"//g;s/$/\/EMU_MMSEQS.idx/')

fi

cat $MMSEQS_IDX > /dev/null