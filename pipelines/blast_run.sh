#!/bin/bash

helpFunction()
{
   echo "Usage: blast_run.sh -p '/path/to/the/directory'-t 16 -m 1400 -M 1800 -i 75"
   echo -e "\t-p <path> Path to directory containing passed raw data"
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

while getopts "p:t:m:M:i:" opt
do
    case "$opt" in
    p )
    	path="$OPTARG"
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

conda_path=$(which conda | sed "s/\b\/conda\b//g")


if [ ! -f $path/barcode_list ]
	then
	for x in {01..24..01}
	do
	echo barcode$x
	done > $path/barcode_list
fi

while read barcode

do

    source $conda_path/activate nanofilt

    zcat $path/$barcode/*fastq.gz | Nanofilt -q 10 -l 1400 --maxlength 1800 | sed -n '1~4s/^@/>/p;2~4p' > $path/$barcode/${barcode}_16s.fasta

    source $conda_path/activate blast

    blastn -db $BLAST_DB/16S_ribosomal_RNA -query $path/$barcode/${barcode}_16s.fasta -out $path/${barcode}/${barcode}_blast.txt -num_threads $threads -max_target_seqs 1 -max_hsps 1 -outfmt "6 std staxids stitle"

    source $conda_path/activate taxonkit

    awk -v i="$identity" 'BEGIN{FS="\t";OFS="\t"}{if($3>=i) print $13}' $path/${barcode}/${barcode}_blast.txt | sort | uniq -c | awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1}' | taxonkit reformat --data-dir $TAXONKIT_DB --taxid-field 1 - | sed 's/;/\t/g' > $path/${barcode}/${barcode}_final_blast_result.txt

done < "$path/barcode_list"

if [ ! -d $path/Blast_Results ]; then

    mkdir $path/Blast_Results

fi

mv $path/barcode*/*_final_blast_result.txt $path/Blast_Results/
