#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: main.sh -d /path/to/the/data -k kit-name -s /path/to/the/script -m 1400 -M 1800 -t Species"
   echo -e "\t-d <str> Path Containing the Raw Data."
   echo -e "\t-k <str> Kit-Name."
   echo -e "\t-s <str> Path Containing the Scripts."
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-t <str> Minimum Taxonomy Rank. [default: Species]"
   exit 1 # Exit script after printing help
}

min=1400
max=1800
rank=Species

while getopts "d:k:s:m:M:t:" opt
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
    t )
        rank="$OPTARG"
        ;;
    ? ) helpFunction ;;
    esac
done

r=$(echo $rank | cut -c1-1)

if [ -z "$data_path" ] && [ -z "$script_path" ];then
    
    echo "Please provide the path to the directory containing the data and scripts";
    
	helpFunction
fi

for i in `ls -d $data_path/barcode*/`
do 
	[ "$(find $i -type f)" ] && echo $i

done | sed "s|$data_path||g;s/\///g" > $data_path/barcode_list

if [ "$(cat $data_path/barcode_list | wc -l)" -eq 0 ]; then

    echo "No Data is written yet."

else

    conda activate parallel


    n=$(cat "$data_path/barcode_list" | wc -l)

    echo $n

    if [ "$n" -gt 8 ]; then parallel_jobs=8; else parallel_jobs="$n"; fi

    while read barcode

    do

    echo "bash $script_path/real_time_kraken.sh -d $data_path -k $kit_name -b $barcode -m $min -M $max -t $r"

    done < "$data_path/barcode_list" | parallel -j "$parallel_jobs" {}
fi
