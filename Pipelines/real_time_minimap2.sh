#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: real_time_analysis_blast.sh -d /path/to/data/directory -k kit-name -b barcode01 -m 1400 -M 1800 -i 85 -q 10 -n REFSEQ"
   echo -e "\t-d <str> Path Containing Sequencing Data."
   echo -e "\t-k <str> Kit-name."
   echo -e "\t-b <str> Barcode Name."
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-i <int> Minimum Percent Identity. [default: 85]"
   echo -e "\t-c <int> Minimum Percent Coverage. [default: 85]"
   echo -e "\t-q <int> Minimum Q-Score. [default: 10]"
   echo -e "\t-n <str> Database Name. [default: REFSEQ]"
   exit 1 # Exit script after printing help
}

min=1400
max=1800
identity=85
coverage=85
q_score=10
db="REFSEQ"

while getopts "d:k:b:m:M:i:c:q:n:" opt
do
    case "$opt" in
    d )
        data_path="$OPTARG"
        ;;
    k )
        kit_name="$OPTARG"
        ;;
    b )
	    barcode="$OPTARG"
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

script_path=$(dirname "$(readlink -f "$0")")

if [ -z "$data_path" ]
    then
    echo "Please provide the path to the Data Directory.";
    helpFunction
fi

if [ "$db" == "REFSEQ" ]; then

    MINIMAP_DB=$(grep MINIMAP2_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_REFSEQ="//;s/"//g;s/$/\/refseq_final_seqs.fasta/')

    TAXA_DATA=$(grep MINIMAP2_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_REFSEQ="//;s/"//g;s/$/\/RefSeq_taxa.txt/')

elif [ "$db" == "GTDB" ]; then
    
    MINIMAP_DB=$(grep MINIMAP2_GTDB ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_GTDB="//;s/"//g;s/$/\/GTDB_final_seqs.fasta/')

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

if [ ! -d $data_path/$barcode/$db ]; then

    mkdir -p $data_path/$barcode/$db

fi


if [ ! -f $data_path/$barcode/processed_files.txt ]; then
    
    conda activate nanofilt
    
    ls $data_path/$barcode/*fastq.gz | xargs zcat | dorado trim --threads 1 --sequencing-kit $kit_name --emit-fastq | NanoFilt -q $q_score --length $min --maxlength $max | sed -n '1~4s/^@/>/p;2~4p' > ${data_path}/$barcode/${barcode}_16S.fasta

    conda activate bbtools

    ls $data_path/$barcode/*fastq.gz | xargs zcat | dorado trim --threads 1 --sequencing-kit $kit_name --emit-fastq | readlength.sh in=stdin.fq out=$data_path/$barcode/${barcode}_hist_temp.txt bin=10 round=t -ignorebadquality

    paste -d "\t" <(echo $barcode) <(grep "#Avg" $data_path/$barcode/${barcode}_hist_temp.txt | cut -f2) > $data_path/$barcode/${barcode}_average_length.txt

    grep -v "#" $data_path/$barcode/${barcode}_hist_temp.txt | awk '{print $1"\t"$2}' > $data_path/$barcode/${barcode}_hist.txt

    ls $data_path/$barcode/*fastq.gz | xargs zcat | dorado trim --threads 1 --sequencing-kit $kit_name --emit-fastq | bioawk -c fastx '{print meanqual($qual)}' > $data_path/$barcode/${barcode}_quality.txt

    if [ "$(grep ">" $data_path/$barcode/${barcode}_16S.fasta | wc -l)" -gt 0 ]; then

        conda activate minimap2

        minimap2 -ax map-ont -t 2 --eqx $MINIMAP_DB $data_path/$barcode/${barcode}_16S.fasta | \
	    samtools view -@ 1 -F 3844 -bS | samtools sort -@ 1 -o $data_path/$barcode/${barcode}_16S.bam; samtools index -@ 1 $data_path/$barcode/${barcode}_16S.bam

        python $script_path/alignment_filter.py -b $data_path/$barcode/${barcode}_16S.bam -t $TAXA_DATA -i $identity -c $coverage | \
        awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) print $1, $4, $5, $6, $7, $8, $9, $1}' | sort -k1 -n -r | uniq > $data_path/$barcode/$db/${barcode}_final_minimap2_result.txt

        rm -r $data_path/$barcode/${barcode}_hist_temp.txt $data_path/$barcode/${barcode}_16S.bam
        
        ls $data_path/$barcode/*fastq.gz > $data_path/$barcode/processed_files.txt
    
    else
        
        echo "No Reads to Classify."
    
    fi

else

    if [ "$(echo $data_path/$barcode/*fastq.gz | tr ' ' '\n' | grep -vf $data_path/$barcode/processed_files.txt | wc -l)" -eq 0 ]; then
        
        echo "No New Files Found."

    else

        conda activate nanofilt

        ls $data_path/$barcode/*fastq.gz | grep -vf $data_path/$barcode/processed_files.txt | xargs zcat | dorado trim --threads 1 --sequencing-kit $kit_name --emit-fastq | NanoFilt -q $q_score --length $min --maxlength $max | sed -n '1~4s/^@/>/p;2~4p' > $data_path/$barcode/${barcode}_16S.fasta

        conda activate bbtools

        ls $data_path/$barcode/*fastq.gz | xargs zcat | dorado trim --threads 1 --sequencing-kit $kit_name --emit-fastq | readlength.sh in=stdin.fq out=$data_path/$barcode/${barcode}_hist_temp.txt bin=10 round=t -ignorebadquality

        paste -d "\t" <(echo $barcode) <(grep "#Avg" $data_path/$barcode/${barcode}_hist_temp.txt | cut -f2) > $data_path/$barcode/${barcode}_average_length.txt

        grep -v "#" $data_path/$barcode/${barcode}_hist_temp.txt | awk '{print $1"\t"$2}' > $data_path/$barcode/${barcode}_hist.txt

        ls $data_path/$barcode/*fastq.gz | grep -vf $data_path/$barcode/processed_files.txt | xargs zcat | dorado trim --threads 1 --sequencing-kit $kit_name --emit-fastq | bioawk -c fastx '{print meanqual($qual)}' >> $data_path/$barcode/${barcode}_quality.txt

        if [ "$(grep ">" $data_path/$barcode/${barcode}_16S.fasta | wc -l)" -gt 0 ]; then

            conda activate minimap2

            minimap2 -ax map-ont -t 2 --eqx $MINIMAP_DB $data_path/$barcode/${barcode}_16S.fasta | \
	        samtools view -@ 1 -F 3844 -bS | samtools sort -@ 1 -o $data_path/$barcode/${barcode}_16S.bam

            python $script_path/alignment_filter.py -b $data_path/$barcode/${barcode}_16S.bam -t $TAXA_DATA -i $identity -c $coverage | \
            awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) print $1, $4, $5, $6, $7, $8, $9, $1}' | sort -k1 -n -r | uniq >> $data_path/$barcode/$db/${barcode}_final_minimap2_result.txt
            
            rm -r $data_path/$barcode/${barcode}_hist_temp.txt $data_path/$barcode/${barcode}_16S.bam

            ls $data_path/$barcode/*fastq.gz | grep -vf $data_path/$barcode/processed_files.txt >> $data_path/$barcode/processed_files.txt
        
        else
            
            echo "No Reads to Classify."
        
        fi
    fi
fi