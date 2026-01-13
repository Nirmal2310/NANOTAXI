#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: real_time_minimap2.sh -d /path/to/data/directory -k kit-name -b barcode01 -m 1400 -M 1800 -i 85 -q 10 -n REFSEQ -t 4 -s 500"
   echo -e "\t-d <str> Path Containing Sequencing Data."
   echo -e "\t-k <str> Kit-name."
   echo -e "\t-b <str> Barcode Name."
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-i <int> Minimum Percent Identity. [default: 85]"
   echo -e "\t-c <int> Minimum Percent Coverage. [default: 85]"
   echo -e "\t-q <int> Minimum Q-Score. [default: 10]"
   echo -e "\t-n <str> Database Name. [default: REFSEQ]"
   echo -e "\t-t <int> Number of Threads."
   echo -e "\t-s <int> Number of reads to process per barcode. [default: 500]"
   exit 1 # Exit script after printing help
}

min=1400
max=1800
identity=85
coverage=85
q_score=10
db="REFSEQ"
batch_size=500

tmp_file=$(mktemp)

while getopts "d:k:b:m:M:i:c:q:n:t:s:" opt
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
    t )
        threads="$OPTARG"
        ;;
    s )
        batch_size="$OPTARG"
        ;;
    ? ) 
        helpFunction
        ;;
    esac
done

script_path=$(dirname "$(readlink -f "$0")")

if [ -z "$data_path" ]
    then
    echo "Please provide the path to the Data Directory.";
    helpFunction
fi

if [ "$db" == "REFSEQ" ]; then

    MINIMAP_DB=$(grep MINIMAP2_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_REFSEQ="//;s/"//g;s/$/\/refseq_final_seqs.mmi/')

    TAXA_DATA=$(grep MINIMAP2_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_REFSEQ="//;s/"//g;s/$/\/RefSeq_taxa.txt/')

elif [ "$db" == "GTDB" ]; then
    
    MINIMAP_DB=$(grep MINIMAP2_GTDB ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_GTDB="//;s/"//g;s/$/\/GTDB_final_seqs.mmi/')

    TAXA_DATA=$(grep MINIMAP2_GTDB ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_GTDB="//;s/"//g;s/$/\/GTDB_taxa.txt/')

elif [ "$db" == "MIMT" ]; then

    MINIMAP_DB=$(grep MINIMAP2_MIMT ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_MIMT="//;s/"//g;s/$/\/MIMT_final_seqs.mmi/')

    TAXA_DATA=$(grep MINIMAP2_MIMT ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_MIMT="//;s/"//g;s/$/\/MIMT_taxa.txt/')

elif [ "$db" == "GSR" ]; then

    MINIMAP_DB=$(grep MINIMAP2_GSR ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_GSR="//;s/"//g;s/$/\/GSR-DB_full-16S_filt_seqs.mmi/')

    TAXA_DATA=$(grep MINIMAP2_GSR ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_GSR="//;s/"//g;s/$/\/GSR-DB_full-16S_filt_taxa.txt/')

elif [ "$db" == "EMUDB" ]; then

    MINIMAP_DB=$(grep MINIMAP2_EMUDB ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_EMUDB="//;s/"//g;s/$/\/EMU.mmi/')

    TAXA_DATA=$(grep MINIMAP2_EMUDB ~/.bashrc | tail -n 1 | sed 's/export MINIMAP2_EMUDB="//;s/"//g;s/$/\/EMU_taxa.txt/')

fi

if [ ! -d $data_path/$barcode/$db ]; then

    mkdir -p $data_path/$barcode/$db

fi


if [ ! -f $data_path/$barcode/$db/processed_reads.txt ]; then

    if [ "$(zcat $data_path/$barcode/*fastq.gz | wc -l | awk '{print $1/4}')" -ge $batch_size ]; then

        zcat $data_path/$barcode/*fastq.gz | paste - - - - | cut -f1 | sed 's/ .*$//g;s/@//g' | head -n $batch_size > $data_path/$barcode/$db/processing_reads.txt
        
        conda activate chopper

        zcat $data_path/$barcode/*fastq.gz | seqkit --threads $threads grep -f $data_path/$barcode/$db/processing_reads.txt - | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | chopper -q $q_score --minlength $min --maxlength $max --threads $threads | sed -n '1~4s/^@/>/p;2~4p' > ${data_path}/$barcode/$db/${barcode}_16S.fasta

        conda activate bbtools
        
        zcat $data_path/$barcode/*fastq.gz | seqkit --threads $threads grep -f $data_path/$barcode/$db/processing_reads.txt - | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | readlength.sh in=stdin.fq out=$data_path/$barcode/$db/${barcode}_hist_temp.txt bin=10 round=t -ignorebadquality

        paste -d "\t" <(echo $barcode) <(grep "#Avg" $data_path/$barcode/$db/${barcode}_hist_temp.txt | cut -f2) > $data_path/$barcode/$db/${barcode}_average_length.txt

        grep -v "#" $data_path/$barcode/$db/${barcode}_hist_temp.txt | awk '{print $1"\t"$2}' > $data_path/$barcode/$db/${barcode}_hist.txt

        zcat $data_path/$barcode/*fastq.gz | seqkit --threads $threads grep -f $data_path/$barcode/$db/processing_reads.txt - | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | bioawk -c fastx '{print meanqual($qual)}' > $data_path/$barcode/$db/${barcode}_quality.txt

        cp $data_path/$barcode/$db/processing_reads.txt $data_path/$barcode/$db/processed_reads.txt && rm -r $data_path/$barcode/$db/processing_reads.txt

        paste -d "\t" <(echo $barcode) <(cat $data_path/$barcode/$db/processed_reads.txt | wc -l) > $data_path/$barcode/$db/${barcode}_processed_reads.txt
    
    else
    
        conda activate chopper
    
        ls $data_path/$barcode/*fastq.gz | xargs zcat | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | chopper -q $q_score --minlength $min --maxlength $max --threads $threads | sed -n '1~4s/^@/>/p;2~4p' > ${data_path}/$barcode/$db/${barcode}_16S.fasta

        conda activate bbtools

        ls $data_path/$barcode/*fastq.gz | xargs zcat | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | readlength.sh in=stdin.fq out=$data_path/$barcode/$db/${barcode}_hist_temp.txt bin=10 round=t -ignorebadquality

        paste -d "\t" <(echo $barcode) <(grep "#Avg" $data_path/$barcode/$db/${barcode}_hist_temp.txt | cut -f2) > $data_path/$barcode/$db/${barcode}_average_length.txt

        grep -v "#" $data_path/$barcode/$db/${barcode}_hist_temp.txt | awk '{print $1"\t"$2}' > $data_path/$barcode/$db/${barcode}_hist.txt

        ls $data_path/$barcode/*fastq.gz | xargs zcat | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | bioawk -c fastx '{print meanqual($qual)}' > $data_path/$barcode/$db/${barcode}_quality.txt

        zcat $data_path/$barcode/*fastq.gz | paste - - - - | cut -f1 | sed 's/ .*$//g;s/@//g' > $data_path/$barcode/$db/processed_reads.txt

        paste -d "\t" <(echo $barcode) <(cat $data_path/$barcode/$db/processed_reads.txt | wc -l) > $data_path/$barcode/$db/${barcode}_processed_reads.txt
    fi

    if [ "$(grep ">" $data_path/$barcode/$db/${barcode}_16S.fasta | wc -l)" -gt 0 ]; then

        conda activate minimap2

        minimap2 -ax map-ont -t $threads -N 50 -p .9 -K 500000000 --eqx $MINIMAP_DB $data_path/$barcode/$db/${barcode}_16S.fasta | \
	    samtools view -@ $threads -F 3844 -bS | samtools sort -@ $threads -o $data_path/$barcode/$db/${barcode}_16S.bam; samtools index -@ $threads $data_path/$barcode/$db/${barcode}_16S.bam

        python $script_path/alignment_filter.py -b $data_path/$barcode/$db/${barcode}_16S.bam -t $TAXA_DATA -i $identity -c $coverage | \
        awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) print $1, $2, $4, $5, $6, $7, $8, $9, $1}' | sort -k1 -n -r | uniq > $data_path/$barcode/$db/${barcode}_final_minimap2_result.txt

        rm -r $data_path/$barcode/$db/${barcode}_hist_temp.txt $data_path/$barcode/$db/${barcode}_16S.bam*
    
    else
        
        echo "No Reads to Classify."
    
    fi

else

    if [ "$(zcat $data_path/$barcode/*fastq.gz | paste - - - - | cut -f1 | sed 's/ .*$//g;s/@//g' | grep -vf $data_path/$barcode/$db/processed_reads.txt | wc -l)" -lt $batch_size ]; then
        
        echo "No New Reads Found."

    else

        zcat $data_path/$barcode/*fastq.gz | paste - - - - | cut -f1 | sed 's/ .*$//g;s/@//g' | grep -vf $data_path/$barcode/$db/processed_reads.txt | head -n $batch_size > $data_path/$barcode/$db/processing_reads.txt

        cat $data_path/$barcode/$db/processing_reads.txt $data_path/$barcode/$db/processed_reads.txt > $tmp_file && cp $tmp_file $data_path/$barcode/$db/processed_reads.txt && rm -r $tmp_file
        
        conda activate chopper

        zcat $data_path/$barcode/*fastq.gz | seqkit --threads $threads grep -f $data_path/$barcode/$db/processing_reads.txt - | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | chopper -q $q_score --minlength $min --maxlength $max --threads $threads | sed -n '1~4s/^@/>/p;2~4p' > $data_path/$barcode/$db/${barcode}_16S.fasta

        conda activate bbtools

        zcat $data_path/$barcode/*fastq.gz | seqkit --threads $threads grep -f $data_path/$barcode/$db/processed_reads.txt - | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | readlength.sh in=stdin.fq out=$data_path/$barcode/$db/${barcode}_hist_temp.txt bin=10 round=t -ignorebadquality

        paste -d "\t" <(echo $barcode) <(grep "#Avg" $data_path/$barcode/$db/${barcode}_hist_temp.txt | cut -f2) > $data_path/$barcode/$db/${barcode}_average_length.txt

        grep -v "#" $data_path/$barcode/$db/${barcode}_hist_temp.txt | awk '{print $1"\t"$2}' > $data_path/$barcode/$db/${barcode}_hist.txt

        zcat $data_path/$barcode/*fastq.gz | seqkit --threads $threads grep -f $data_path/$barcode/$db/processed_reads.txt - | dorado trim --threads $threads --sequencing-kit $kit_name --emit-fastq | bioawk -c fastx '{print meanqual($qual)}' > $data_path/$barcode/$db/${barcode}_quality.txt

        rm -r $data_path/$barcode/$db/processing_reads.txt

        paste -d "\t" <(echo $barcode) <(cat $data_path/$barcode/$db/processed_reads.txt | wc -l) > $data_path/$barcode/$db/${barcode}_processed_reads.txt

        if [ "$(grep ">" $data_path/$barcode/$db/${barcode}_16S.fasta | wc -l)" -gt 0 ]; then

            conda activate minimap2

            minimap2 -ax map-ont -t $threads -N 50 -p .9 -K 500000000 --eqx $MINIMAP_DB $data_path/$barcode/$db/${barcode}_16S.fasta | \
	        samtools view -@ $threads -F 3844 -bS | samtools sort -@ $threads -o $data_path/$barcode/$db/${barcode}_16S.bam; samtools index -@ $threads $data_path/$barcode/$db/${barcode}_16S.bam

            python $script_path/alignment_filter.py -b $data_path/$barcode/$db/${barcode}_16S.bam -t $TAXA_DATA -i $identity -c $coverage | \
            awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) print $1, $2, $4, $5, $6, $7, $8, $9, $1}' | sort -k1 -n -r | uniq >> $data_path/$barcode/$db/${barcode}_final_minimap2_result.txt
            
            rm -r $data_path/$barcode/$db/${barcode}_hist_temp.txt $data_path/$barcode/$db/${barcode}_16S.bam*
        
        else
            
            echo "No Reads to Classify."
        
        fi
    fi
fi