#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
echo "Usage: blast_run.sh -p /path/to/the/directory -k kit-name -s /path/to/scripts -t 16 -m 1400 -M 1800 -i 85 -c 85 -q 10 -n REFSEQ -e metadata.csv"
echo -e "\t-p <path> Path to directory containing passed raw data."
echo -e "\t-k <str> Kit-Name."
echo -e "\t-s <str> Path Containing the Scripts."
echo -e "\t-t <int> Number of threads to be used for the analysis. [default: 16]"
echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
echo -e "\t-i <int> Minimum BLAST Identity(%). [default: 85]"
echo -e "\t-c <int> Minimum BLAST Coverage(%). [default: 85]"
echo -e "\t-q <int> Minimum Q-Score. [default: 10]"
echo -e "\t-n <str> Database Name. [default: REFSEQ]"
echo -e "\t-e <str> Metadata File."
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

while getopts "p:k:s:t:m:M:i:c:q:n:e:" opt
do
    case "$opt" in
    p )
        path="$OPTARG"
        ;;
    k )
        kit_name="$OPTARG"
        ;;
    s )
        script_path="$OPTARG"
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
    e )
        metadata="$OPTARG"
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

    done | sed "s|$path||g;s/\///g" | grep -f <(awk -F "," '{if(NR>1) print $1}' $metadata) - > $path/barcode_list
fi

if [ "$db" == "REFSEQ" ]; then

    MINIMAP_DB=$(grep MINIMAP2_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_REFSEQ="//;s/"//g;s/$/\/refseq_final_seqs.fasta/')

    TAXA_DATA=$(grep MINIMAP2_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_REFSEQ="//;s/"//g;s/$/\/RefSeq_taxa.txt/')

elif [ "$db" == "GTDB" ]; then
    
    MINIMAP_DB=$(grep MINIMAP2_GTDB ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_GTDB="//;s/"//g;s/$/\/GTBD_final_seqs.fasta/')

    TAXA_DATA=$(grep MINIMAP2_GTDB ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_GTDB="//;s/"//g;s/$/\/GTDB_taxa.txt/')

elif [ "$db" == "MIMT" ]; then

    MINIMAP_DB=$(grep MINIMAP2_MIMT ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_MIMT="//;s/"//g;s/$/\/MIMT_final_seqs.fasta/')

    TAXA_DATA=$(grep MINIMAP2_MIMT ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_MIMT="//;s/"//g;s/$/\/MIMT_taxa.txt/')

elif [ "$db" == "GSR" ]; then

    MINIMAP_DB=$(grep MINIMAP2_GSR ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_GSR="//;s/"//g;s/$/\/GSR-DB_full-16S_filt_seqs.fasta/')

    TAXA_DATA=$(grep MINIMAP2_GSR ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_GSR="//;s/"//g;s/$/\/GSR-DB_full-16S_filt_taxa.txt/')

elif [ "$db" == "EMUDB" ]; then

    MINIMAP_DB=$(grep MINIMAP2_EMUDB ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_EMUDB="//;s/"//g;s/$/\/EMU.fasta/')

    TAXA_DATA=$(grep MINIMAP2_EMUDB ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_EMUDB="//;s/"//g;s/$/\/EMU_taxa.txt/')

fi

while read barcode

do

    conda activate nanofilt

    zcat $path/$barcode/*fastq.gz | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | NanoFilt -q $q_score -l $min --maxlength $max | sed -n '1~4s/^@/>/p;2~4p' > $path/$barcode/${barcode}_16s.fasta

    conda activate minimap2

    minimap2 -ax map-ont -t $threads -I 16G --eqx $MINIMAP_DB $path/$barcode/${barcode}_16s.fasta | \
        samtools view -@ $threads -F 3844 -bS | samtools sort -@ $threads -o $path/$barcode/${barcode}_16S.bam; samtools index -@ $threads $path/$barcode/${barcode}_16S.bam
    
    python $script_path/alignment_filter.py -b $path/$barcode/${barcode}_16S.bam -t $TAXA_DATA -i $identity -c $coverage | \
        awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) print $1, $2, $4, $5, $6, $7, $8, $9, $1}' | sort -k1 -n -r | uniq > $path/$barcode/${barcode}_final_minimap2_result.txt

done < "$path/barcode_list"

if [ ! -d $path/Minimap2_Results/$db ]; then

    mkdir -p $path/Minimap2_Results/$db

fi

mv $path/barcode*/*_final_minimap2_result.txt $path/Minimap2_Results/$db/