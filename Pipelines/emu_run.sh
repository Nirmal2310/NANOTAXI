#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: emu_run.sh -p /path/to/the/directory -k kit-name -t 16 -m 1400 -M 1800 -q 10 -n EMUDB"
   echo -e "\t-p <path> Path to directory containing passed raw data."
   echo -e "\t-k <str> Kit-Name."
   echo -e "\t-t <int> Number of threads to be used for the analysis. [default: 16]"
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-q <int> Minimum Q-Score. [default: 10]"
   echo -e "\t-n <str> Database Name. [default: EMUDB]"
   exit 1 # Exit script after printing help
}

# Default values for EMU

threads=16
min=1400
max=1800
q_score=10
db="EMUDB"

while getopts "p:k:t:m:M:q:n:" opt
do
    case "$opt" in
    p )
    	path="$OPTARG"
    	;;
    k )
        kit_name="$OPTARG"
        ;;
    t )
        threads="$OPTARG"
        ;;
    m )
    	min="$OPTARG"
    	;;
    M )
    	max="$OPTARG"
    	;;
    q )
        q_score="$OPTARG"
        ;;
    n )
        db="$OPTARG"
        ;;
    ? ) helpFunction ;;
    esac
done

if [ -z "$path" ]
    then
    echo "Please provide the path to the directory containing raw data";
    helpFunction
fi


if [ ! -f $path/barcode_list ]
	then
	for i in `ls -d $path/barcode*/`
    do 
	    [ "$(find $i -type f)" ] && echo $i

    done | sed "s|$path||g;s/\///g" > $path/barcode_list
fi

TAXONKIT_DB=$(grep TAXONKIT_DB ~/.bashrc | tail -n 1 | sed 's/export TAXONKIT_DB="//;s/"//g')

if [ "$db" == "REFSEQ" ]; then

    EMU_DB=$(grep EMU_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export EMU_REFSEQ="//;s/"//g')

elif [ "$db" == "GTDB" ]; then
    
    EMU_DB=$(grep EMU_GTDB ~/.bashrc | tail -n 1 | sed 's/export EMU_GTDB="//;s/"//g')

elif [ "$db" == "MIMT" ]; then

    EMU_DB=$(grep EMU_MIMT ~/.bashrc | tail -n 1 | sed 's/export EMU_MIMT="//;s/"//g')

elif [ "$db" == "GSR" ]; then

    EMU_DB=$(grep EMU_GSR ~/.bashrc | tail -n 1 | sed 's/export EMU_GSR="//;s/"//g')

elif [ "$db" == "EMUDB" ]; then

    EMU_DB=$(grep EMU_DB ~/.bashrc | tail -n 1 | sed 's/export EMU_DB="//;s/"//g')

fi

while read barcode

do

    conda activate nanofilt

    zcat $path/$barcode/*fastq.gz | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | NanoFilt -q $q_score -l $min --maxlength $max | sed -n '1~4s/^@/>/p;2~4p' > $path/$barcode/${barcode}_16s.fasta

    conda activate emu

    emu abundance --type map-ont --db $EMU_DB --output-dir $path/$barcode/ --output-basename $barcode --keep-counts --threads $threads $path/$barcode/${barcode}_16s.fasta

    conda activate taxonkit

    line_count=$(wc -l $path/$barcode/${barcode}_rel-abundance.tsv | awk '{print $1}')

    awk -v l=$line_count 'BEGIN{FS="\t";OFS="\t"}{if(NR>1 && NR<l) print $1,$10,$9,$8,$7,$6,$5,$4,$3; else if(NR>1 && NR==l) print "Unclassified",$10,"Unclassified","Unclassified","Unclassified","Unclassified","Unclassified","Unclassified","Unclassified"}' $path/$barcode/${barcode}_rel-abundance.tsv > $path/$barcode/${barcode}_final_emu_result.txt

done < "$path/barcode_list"

if [ ! -d $path/EMU_Results/$db ]; then

    mkdir -p $path/EMU_Results/$db

fi

echo "Success"

mv $path/barcode*/*_final_emu_result.txt $path/EMU_Results/$db/
