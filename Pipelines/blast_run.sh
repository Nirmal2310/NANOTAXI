#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: blast_run.sh -p /path/to/the/directory -k kit-name -t 16 -m 1400 -M 1800 -i 85 -c 85 -q 10 -n REFSEQ"
   echo -e "\t-p <path> Path to directory containing passed raw data."
   echo -e "\t-k <str> Kit-Name."
   echo -e "\t-t <int> Number of threads to be used for the analysis. [default: 16]"
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-i <int> Minimum BLAST Identity(%). [default: 85]"
   echo -e "\t-c <int> Minimum BLAST Coverage(%). [default: 85]"
   echo -e "\t-q <int> Minimum Q-Score. [default: 10]"
   echo -e "\t-n <str> Database Name. [default: REFSEQ]"
   exit 1 # Exit script after printing help
}

# Default values for BLASTN

threads=16
min=1400
max=1800
identity=85
coverage=85
q_score=10
db="REFSEQ"

while getopts "p:k:t:m:M:i:c:q:n:" opt
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
    i )
    	identity="$OPTARG"
    	;;
    c )
        coverage="$OPTARG"
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

if [ "$db" == "REFSEQ" ]; then

    BLAST_DB=$(grep BLAST_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export BLAST_REFSEQ="//;s/"//g;s/$/\/REFSEQ_BLAST/')

    TAXA_DATA=$(grep BLAST_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export BLAST_REFSEQ="//;s/"//g;s/$/\/RefSeq_taxa.txt/')

elif [ "$db" == "GTDB" ]; then
    
    BLAST_DB=$(grep BLAST_GTDB ~/.bashrc | tail -n 1 | sed 's/export BLAST_GTDB="//;s/"//g;s/$/\/GTBD_BLAST/')

    TAXA_DATA=$(grep BLAST_GTDB ~/.bashrc | tail -n 1 | sed 's/export BLAST_GTDB="//;s/"//g;s/$/\/GTDB_taxa.txt/')

elif [ "$db" == "MIMT" ]; then

    BLAST_DB=$(grep BLAST_MIMT ~/.bashrc | tail -n 1 | sed 's/export BLAST_MIMT="//;s/"//g;s/$/\/MIMT_BLAST/')

    TAXA_DATA=$(grep BLAST_MIMT ~/.bashrc | tail -n 1 | sed 's/export BLAST_MIMT="//;s/"//g;s/$/\/MIMT_taxa.txt/')

elif [ "$db" == "GSR" ]; then

    BLAST_DB=$(grep BLAST_GSR ~/.bashrc | tail -n 1 | sed 's/export BLAST_GSR="//;s/"//g;s/$/\/GSR_BLAST/')

    TAXA_DATA=$(grep GSR ~/.bashrc | tail -n 1 | sed 's/export GSR="//;s/"//g;s/$/\/GSR_taxa.txt/')

elif [ "$db" == "EMUDB" ]; then

    BLAST_DB=$(grep BLAST_EMU ~/.bashrc | tail -n 1 | sed 's/export BLAST_EMU="//;s/"//g;s/$/\/EMU_BLAST/')

    TAXA_DATA=$(grep BLAST_EMU ~/.bashrc | tail -n 1 | sed 's/export BLAST_EMU="//;s/"//g;s/$/\/EMU_taxa.txt/')

fi

TAXONKIT_DB=$(grep TAXONKIT_DB ~/.bashrc | tail -n 1 | sed 's/export TAXONKIT_DB="//;s/"//g')


if [ ! -f $path/barcode_list ]
	then
	for i in `ls -d $path/barcode*/`
    do 
	    [ "$(find $i -type f)" ] && echo $i

    done | sed "s|$path||g;s/\///g" > $path/barcode_list
fi

script_dir=$(dirname "$(readlink -f "$0")")

while read barcode

do

    conda activate nanofilt

    zcat $path/$barcode/*fastq.gz | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | NanoFilt -q $q_score -l $min --maxlength $max | sed -n '1~4s/^@/>/p;2~4p' > $path/$barcode/${barcode}_16s.fasta

    conda activate blast

    blastn -db $BLAST_DB -query $path/$barcode/${barcode}_16s.fasta -out $path/${barcode}/${barcode}_blast.txt -num_threads $threads -max_target_seqs 1 -max_hsps 1 -perc_identity $identity -qcov_hsp_perc $coverage -outfmt "6"

    conda activate minimap2
    
    python $script_dir/add_taxon_info.py -c <(awk 'BEGIN{FS="\t";OFS="\t"}{print $2}' $path/${barcode}/${barcode}_blast.txt | sort | uniq -c | awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1}') -t $TAXA_DATA | \
    awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) print $1, $2, $4, $5, $6, $7, $8, $9, $1}' | sort -k1 -n -r | uniq > $path/${barcode}/${barcode}_final_blast_result.txt

done < "$path/barcode_list"

if [ ! -d $path/Blast_Results/$db ]; then

    mkdir -p $path/Blast_Results/$db

fi

mv $path/barcode*/*_final_blast_result.txt $path/Blast_Results/$db/