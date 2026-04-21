#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: create_consensus_sequences.sh -i input.fasta -t taxa_file.txt -c consensus.fasta -n 16 -o reference_map.txt"
   echo -e "\t-i <str> Input Reference FASTA File."
   echo -e "\t-t <str> Input Reference taxa File."
   echo -e "\t-c <str> Output consensus FASTA File."
   echo -e "\t-n <int> Number of Threads."
   echo -e "\t-o <str> Output Reference consensus map File."
   exit 1 # Exit script after printing help
}

while getopts "i:t:c:n:o:" opt
do
    case "$opt" in
    i )
        input_fasta="$OPTARG"
        ;;
    t )
        taxa_file="$OPTARG"
        ;;
    c )
        consensus_fasta="$OPTARG"
        ;;
    n )
        threads="$OPTARG"
        ;;
    o )
        output_file="$OPTARG"
        ;;
    ? )
        helpFunction
        ;;
    esac
done

if [ -z "$input_fasta" ] || [ -z "$taxa_file" ] || [ -z "$output_file" ] || [ -z "$threads" ] || [ -z $consensus_fasta ]
    then
    echo "Please provide the Input Reference FASTA, Input Reference Taxa file, Output FASTA file, Number of threads and Output map file."
    helpFunction
fi

tmp_counts=$(mktemp)

tmp_out=$(mktemp)

tmp_singleton=$(mktemp)

tmp_multiton=$(mktemp)

script_dir=$(dirname "$(readlink -f "$0")")

awk -F "\t" '{if(NR>1) print $8}' $taxa_file | sort | uniq -c | sed 's/^[^0-9]*//' | sed 's/ /\t/' > $tmp_counts

conda activate bbtools

if [ -f $input_fasta.fai ]; then
    rm -r $input_fasta.fai && seqkit faidx $input_fasta
else 
    seqkit faidx $input_fasta
fi

while IFS=$'\t' read -r count species; do
    sp_name=$(echo $species | sed "s/ /_/g;s/(/_/g;s/)/_/g;s/\:/_/g;s/\//_/g;s/$/_consensus/;s/'//g")
    if [ $count -eq 2 ]; then
        echo "grep -F \"$species\" $taxa_file | awk -F \"\t\" '{print \$1}' | grep -Ff - $input_fasta.fai | awk 'BEGIN{FS=OFS=\"\t\"}{print \$1,\$2}' | \
        sort -k2 -n -r | head -n 1 | sed 's/\t.*$//g' | seqkit grep -f - --skip-file-check --quiet -j 1 $input_fasta | sed 's/>.*$/>$sp_name/'"
        echo -e "$species\t$sp_name" >> $output_file
    elif [ $count -gt 2 ]; then
        echo "grep -F \"$species\" $taxa_file | awk -F \"\t\" '{print \$1}' | seqkit grep -f - -j 1 --skip-file-check --quiet $input_fasta > ${sp_name}_tmp.fasta && \
              python $script_dir/get_medoid_sequence.py ${sp_name}_tmp.fasta ${sp_name}_tmp_rep.fasta && \
              sed -i 's/>.*$/>$sp_name/g' ${sp_name}_tmp_rep.fasta && cat ${sp_name}_tmp_rep.fasta && rm -r ${sp_name}_tmp.fasta ${sp_name}_tmp_rep.fasta"
        echo -e "$species\t$sp_name" >> $output_file
    fi
done < "$tmp_counts" | parallel -j $threads {} > $tmp_out

awk 'BEGIN{FS=OFS="\t"} NR==FNR {if($1==1) a[$2]=$0; next} $8 in a {print a[$8], $0}' $tmp_counts $taxa_file | awk 'BEGIN{FS="\t";OFS="\t"}{sp_name=$2; gsub(" ","_", sp_name); gsub("$","_consensus", sp_name); print $3,$2,sp_name}' > $tmp_singleton

awk -F "\t" '{print $1}' $tmp_singleton | seqkit grep -f - -j $threads $input_fasta | \
seqkit replace -j $threads -p "^(\S+)" -r "{kv}" -k <(awk -F "\t" '{print $1"\t"$3}' $tmp_singleton) >> $tmp_out

mv $tmp_out $consensus_fasta

awk 'BEGIN{FS=OFS="\t"}{print $2,$3}' $tmp_singleton >> $output_file

sort $output_file | uniq > temp && mv temp $output_file

rm -r "$tmp_counts" "$tmp_singleton" "$input_fasta.fai"

echo "Finished creating consensus sequence per species."