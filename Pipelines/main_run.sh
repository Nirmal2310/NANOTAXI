#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: main_run.sh -d /path/to/the/data -k kit-name -s /path/to/the/script -m 1400 -M 1800 -i 85 -c 85 -q 10 -n REFSEQ -t 24 -p Minimap2 -r Species -v 0.0 -b 500"
   echo -e "\t-d <str> Path Containing the Raw Data."
   echo -e "\t-k <str> Kit-Name."
   echo -e "\t-s <str> Path Containing the Scripts."
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-i <int> Minimum Percent Identity. [default: 85]"
   echo -e "\t-c <int> Minimum Percent Coverage. [default: 85]"
   echo -e "\t-q <int> Minimum Q-Score. [default: 10]"
   echo -e "\t-n <str> Database Name. [default: REFSEQ]"
   echo -e "\t-t <int> Total Number of Threads. [default: 24]"
   echo -e "\t-p <str> Analysis Pipeline. [default: Minimap2]"
   echo -e "\t-r <str> Minimum Taxonomy Rank. [default: Species]"
   echo -e "\t-v <int> Confidence Value. [default: 0.0]"
   echo -e "\t-b <int> Number of reads to process per barcode. [default: 500]"
   exit 1 # Exit script after printing help
}

min=1400
max=1800
identity=85
coverage=85
q_score=10
db="REFSEQ"
threads=24
pipeline="Minimap2"
rank=Species
conf=0.0
batch_size=500

while getopts "d:k:s:m:M:i:c:q:n:t:p:r:v:b:" opt
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
    t )
        threads="$OPTARG"
        ;;
    p )
        pipeline="$OPTARG"
        ;;
    r )
        rank="$OPTARG"
        ;;
    v )
        conf="$OPTARG"
        ;;
    b )
        batch_size="$OPTARG"
        ;;
    ? ) helpFunction ;;
    esac
done

r=$(echo $rank | cut -c1-1)

if [ -z "$data_path" ] && [ -z "$script_path" ];then
    
	echo "Please provide the path to the directory containing the data and scripts";
    
	helpFunction
fi


if ls -d "$data_path"/barcode*/ >/dev/null 2>&1; then
    for i in "$data_path"/barcode*/; do
        if [ "$(find "$i" -type f -name '*.gz' | wc -l)" -gt 0 ]; then
            if [ -f "$i"/"$db"/processed_reads.txt ]; then
                barcode_name=$(echo "$i" | sed "s|$data_path||g;s|/||g")
                if [ "$(zcat $data_path/$barcode_name/*fastq.gz | paste - - - - | cut -f1 | sed 's/ .*$//g;s/@//g' | grep -vf $data_path/$barcode_name/$db/processed_reads.txt | wc -l)" -ge $batch_size ]; then
                    echo "$i"
                fi
            else
                barcode_name=$(echo "$i" | sed "s|$data_path||g;s|/||g")
                if [ "$(zcat $data_path/$barcode_name/*fastq.gz | wc -l | awk '{print $1/4}')" -ge $batch_size ]; then
                    echo "$i"
                else
                    continue
                fi
            fi    
        fi

    done | sed "s|$data_path||g;s|/||g" > "$data_path"/barcode_list
fi

if [ ! -s "$data_path/barcode_list" ]; then

    echo "No new files found."

else

    echo "New files found, now processing."

    conda activate seqkit

    n=$(cat "$data_path/barcode_list" | wc -l)

    max_process=$(("$threads"/4))

    if [ "$n" -ge "$max_process" ]; then
        per_job_thread=4
        parallel_jobs=$max_process
    else
        per_job_thread=$(("$threads"/"$n"))
        parallel_jobs=$n
    fi

    while read -r barcode

    do

    	if [ $pipeline == "Minimap2" ]; then

            echo "bash $script_path/real_time_minimap2.sh -d $data_path -k $kit_name -b $barcode -m $min -M $max -i $identity -c $coverage -q $q_score -n $db -t $per_job_thread -s $batch_size > /dev/null"
        
        elif [ $pipeline == "Kraken2" ]; then
            
            echo "bash $script_path/real_time_kraken.sh -d $data_path -k $kit_name -b $barcode -m $min -M $max -t $r -c $conf -q $q_score -n $db -p $per_job_thread -s $batch_size > /dev/null"
        
        elif [ $pipeline == "BLASTn" ]; then

            echo "bash $script_path/real_time_blash.sh -d $data_path -k $kit_name -b $barcode -m $min -M $max -i $identity -c $coverage -q $q_score -n $db -t $per_job_thread -s $batch_size > /dev/null"
        
        elif [ $pipeline == "EMU" ]; then

            echo "bash $script_path/real_time_emu.sh -d $data_path -k $kit_name -b $barcode -m $min -M $max -q $q_score -n $db -t $per_job_thread -s $batch_size > /dev/null"

        fi


    done < "$data_path/barcode_list" | parallel -j "$parallel_jobs" {}
fi