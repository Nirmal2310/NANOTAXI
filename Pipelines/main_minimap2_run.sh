#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: main_minimap2_run.sh -d /path/to/the/data -k kit-name -s /path/to/the/script -m 1400 -M 1800 -i 85 -c 85 -q 10 -n REFSEQ -e metadata.csv"
   echo -e "\t-d <str> Path Containing the Raw Data."
   echo -e "\t-k <str> Kit-Name."
   echo -e "\t-s <str> Path Containing the Scripts."
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-i <int> Minimum Percent Identity. [default: 85]"
   echo -e "\t-c <int> Minimum Percent Coverage. [default: 85]"
   echo -e "\t-q <int> Minimum Q-Score. [default: 10]"
   echo -e "\t-n <str> Database Name. [default: REFSEQ]"
   echo -e "\t-e <str> Metadata File."
   exit 1 # Exit script after printing help
}

min=1400
max=1800
identity=85
coverage=85
q_score=10
db="REFSEQ"

while getopts "d:k:s:m:M:i:c:q:n:e:" opt
do
    case "$opt" in
    d )
        data_path="$OPTARG"
        ;;
    k )
        kit_name="$OPTARG"
        ;;
    s )
	    script_path="$OPTARG"
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
    e )
        metadata="$OPTARG"
        ;;
    ? ) helpFunction ;;
    esac
done

if [ -z "$data_path" ] && [ -z "$script_path" ];then
    
	echo "Please provide the path to the directory containing the data and scripts";
    
	helpFunction
fi

for i in `ls -d $data_path/barcode*/`
do 
	[ "$(find $i -type f)" ] && echo $i

done | sed "s|$data_path||g;s/\///g" | grep -f <(awk -F "," '{if(NR>1) print $1}' $metadata) - > $data_path/barcode_list

if [ "$(cat $data_path/barcode_list | wc -l)" -eq 0 ]; then

    echo "No Data is written yet."

else

    conda activate seqkit

    n=$(cat "$data_path/barcode_list" | wc -l)

    if [ "$n" -gt 8 ]; then parallel_jobs=8; else parallel_jobs="$n"; fi

    while read barcode

    do

    	echo "bash $script_path/real_time_minimap2.sh -d $data_path -k $kit_name -b $barcode -m $min -M $max -i $identity -c $coverage -q $q_score -n $db"

    done < "$data_path/barcode_list" | parallel -j "$parallel_jobs" {}
fi