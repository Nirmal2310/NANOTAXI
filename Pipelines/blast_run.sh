#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: blast_run.sh -p /path/to/the/directory -k kit-name -t 16 -m 1400 -M 1800 -i 75"
   echo -e "\t-p <path> Path to directory containing passed raw data."
   echo -e "\t-k <str> Kit-Name."
   echo -e "\t-t <int> Number of threads to be used for the analysis. [default: 16]"
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-i <int> Minimum BLAST Identity(%). [default: 75]"
   exit 1 # Exit script after printing help
}

# Default values for BLASTN

threads=16
min=1400
max=1800
identity=75

BLAST_DB=$(grep BLAST_DB ~/.bashrc | tail -n 1 | sed 's/export BLAST_DB="//;s/"//g')

TAXONKIT_DB=$(grep TAXONKIT_DB ~/.bashrc | tail -n 1 | sed 's/export TAXONKIT_DB="//;s/"//g')

echo $BLAST_DB

while getopts "p:k:t:m:M:i:" opt
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

while read barcode

do

    conda activate nanofilt

    zcat $path/$barcode/*fastq.gz | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | NanoFilt -q 10 -l $min --maxlength $max | sed -n '1~4s/^@/>/p;2~4p' > $path/$barcode/${barcode}_16s.fasta

    conda activate blast

    blastn -db $BLAST_DB/16S_ribosomal_RNA -query $path/$barcode/${barcode}_16s.fasta -out $path/${barcode}/${barcode}_blast.txt -num_threads $threads -max_target_seqs 1 -max_hsps 1 -outfmt "6 std staxids stitle"

    conda activate taxonkit

    awk -v i="$identity" 'BEGIN{FS="\t";OFS="\t"}{if($3>=i) print $13}' $path/${barcode}/${barcode}_blast.txt | sort | uniq -c | awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1}' | taxonkit reformat --threads $threads --data-dir $TAXONKIT_DB --taxid-field 1 - | sed 's/;/\t/g' > $path/${barcode}/${barcode}_final_blast_result.txt

done < "$path/barcode_list"

if [ ! -d $path/Blast_Results ]; then

    mkdir $path/Blast_Results

fi

mv $path/barcode*/*_final_blast_result.txt $path/Blast_Results/
