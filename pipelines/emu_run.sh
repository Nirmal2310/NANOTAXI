#!/bin/bash

helpFunction()
{
   echo "Usage: emu_run.sh -p '/path/to/the/directory'-t 16 -m 1400 -M 1800"
   echo -e "\t-p <path> Path to directory containing passed raw data"
   echo -e "\t-t <int> Number of threads to be used for the analysis. [default: 16]"
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   exit 1 # Exit script after printing help
}

# Default values for BLASTN

threads=16
min=1400
max=1800

while getopts "p:t:m:M:" opt
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

    source $conda_path/activate emu

    emu abundance --type map-ont --db $EMU_DB --output-dir $path/$barcode/ --output-basename $barcode --keep-counts --threads $threads $path/$barcode/${barcode}_16s.fasta

    source $conda_path/activate taxonkit

    awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) print $1,$14}' $path/$barcode/${barcode}_rel-abundance.tsv | sed '$ d' | taxonkit reformat --data-dir $TAXONKIT_DB --taxid-field 1 - | sed 's/;/\t/g' > $path/$barcode/temp

    awk 'BEGIN {for (i = 1; i <= 7; i++) printf "%s\t", "Unclassified"}' | paste -d "\t" <(echo "Unclassified") <(tail -n 1 $path/$barcode/${barcode}_rel-abundance.tsv | awk -F "\t" '{print $14}') - | cat $path/$barcode/temp - > $path/$barcode/${barcode}_final_emu_result.txt

    rm -r $path/$barcode/temp

done < "$path/barcode_list"

if [ ! -d $path/EMU_Results ]; then

    mkdir $path/EMU_Results

fi

mv $path/barcode*/*_final_emu_result.txt $path/EMU_Results/