#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: kraken_run.sh -p '/path/to/the/directory'-t 16 -m 1400 -M 1800 -r Species"
   echo -e "\t-p <path> Path to directory containing passed raw data"
   echo -e "\t-t <int> Number of threads to be used for the analysis. [default: 16]"
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-r <str> Minimum Taxonomy Rank. [default: Species]"
   exit 1 # Exit script after printing help
}

# Default values for BLASTN

threads=16
min=1400
max=1800
rank=Species

while getopts "p:t:m:M:r:" opt
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
    r )
    	rank="$OPTARG"
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
	for x in {01..24..01}
	do
	echo barcode$x
	done > $path/barcode_list
fi

r=$(echo $rank | cut -c1-1)

while read barcode

do

    conda activate nanofilt

    zcat $path/$barcode/*fastq.gz | Nanofilt -q 10 -l $min --maxlength $max | sed -n '1~4s/^@/>/p;2~4p' > $path/$barcode/${barcode}_16s.fasta

    conda activate kraken2

    kraken2 --db $KRAKEN_DB --threads $threads $path/$barcode/${barcode}_16s.fasta --output $path/$barcode/${barcode}_kraken2_output.txt --report $path/$barcode/${barcode}_kraken2_report.txt

    kraken-biom --max $r --min $r -o $path/$barcode/${barcode}_kraken_biom.txt --fmt tsv $path/$barcode/${barcode}_kraken2_report.txt

    conda activate taxonkit

    sed '1,2d' $path/$barcode/${barcode}_kraken_biom.txt | taxonkit reformat --data-dir $TAXONKIT_DB --taxid-field 1 - | sed 's/;/\t/g' > $path/$barcode/${barcode}_final_kraken2_result.txt

    rm -r $path/$barcode/${barcode}_kraken_biom.txt

done < "barcode_list"

if [ ! -d $path/Kraken2_Results ]; then
    
    mkdir $path/Kraken2_Results

fi

mv $path/barcode*/*_final_kraken2_result.txt $path/Kraken2_Results/