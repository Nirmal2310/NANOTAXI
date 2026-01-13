#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: real_time_kraken.sh -d /path/to/data/directory -k kit-name -b barcode01 -m 1400 -M 1800 -t Species -c 0.0 -q 10 -n REFSEQ -p 4 -s 500"
   echo -e "\t-d <str> Path Containing Sequencing Data."
   echo -e "\t-k <str> Kit-name."
   echo -e "\t-b <str> Barcode Name."
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-r <str> Minimum Taxonomy Rank. [default: Species]"
   echo -e "\t-c <int> Confidence Score. [default: 0.0]"
   echo -e "\t-q <int> Minimum Q-Score. [default: 10]"
   echo -e "\t-n <str> Database Name. [default: REFSEQ]"
   echo -e "\t-p <int> Number of Threads."
   echo -e "\t-s <int> Number of reads to process per barcode. [default: 500]"
   exit 1 # Exit script after printing help
}

min=1400
max=1800
rank=S
conf=0.0
q_score=10
db="REFSEQ"
batch_size=500

tmp_file=$(mktemp)

while getopts "d:k:b:m:M:t:c:q:n:p:s:" opt
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
    t )
    	rank="$OPTARG"
    	;;
    c)
        conf="$OPTARG"
        ;;
    q )
        q_score="$OPTARG"
        ;;
    n )
        db="$OPTARG"
        ;;
    p )
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

if [ -z "$data_path" ]
    then
    echo "Please provide the path to the Data Directory.";
    helpFunction
fi

TAXONKIT_DB=$(grep TAXONKIT_DB ~/.bashrc | tail -n 1 | sed 's/export TAXONKIT_DB="//;s/"//g')

if [ "$db" == "REFSEQ" ]; then

    KRAKEN_DB=$(grep KRAKEN_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export KRAKEN_REFSEQ="//;s/"//g')

elif [ "$db" == "GTDB" ]; then
    
    KRAKEN_DB=$(grep KRAKEN_GTDB ~/.bashrc | tail -n 1 | sed 's/export KRAKEN_GTDB="//;s/"//g')

elif [ "$db" == "MIMT" ]; then

    KRAKEN_DB=$(grep KRAKEN_MIMT ~/.bashrc | tail -n 1 | sed 's/export KRAKEN_MIMT="//;s/"//g')

elif [ "$db" == "GSR" ]; then

    KRAKEN_DB=$(grep KRAKEN_GSR ~/.bashrc | tail -n 1 | sed 's/export KRAKEN_GSR="//;s/"//g')

    TAXONKIT_DB=$(grep KRAKEN_GSR ~/.bashrc | tail -n 1 | sed 's/export KRAKEN_GSR="//;s/"//g;s/$/\/taxonomy\//g')

elif [ "$db" == "EMUDB" ]; then

    KRAKEN_DB=$(grep KRAKEN_EMUDB ~/.bashrc | tail -n 1 | sed 's/export KRAKEN_EMUDB="//;s/"//g')

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

        conda activate kraken2

        kraken2 --db $KRAKEN_DB --confidence $conf --threads $threads --output $data_path/$barcode/$db/${barcode}_kraken2_output.txt --report $data_path/$barcode/$db/${barcode}_kraken2_report.txt $data_path/$barcode/$db/${barcode}_16S.fasta

        kraken-biom --max P --min $rank -o $data_path/$barcode/$db/${barcode}_kraken_biom.txt --fmt tsv $data_path/$barcode/$db/${barcode}_kraken2_report.txt

        conda activate taxonkit

        sed '1,2d' $data_path/$barcode/$db/${barcode}_kraken_biom.txt | taxonkit reformat --threads $threads --data-dir $TAXONKIT_DB --taxid-field 1 - | \
        sed 's/;/\t/g;s/[k,p,c,o,f,g,s]__//g' | awk 'BEGIN{FS="\t";OFS="\t"}{for (i=2;i<=NF;i++) gsub(/_/, " ", $i)} 1' > $data_path/$barcode/$db/${barcode}_final_kraken2_result.txt

        rm -r $data_path/$barcode/$db/${barcode}_kraken_biom.txt $data_path/$barcode/$db/${barcode}_hist_temp.txt
        
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

            conda activate kraken2

            kraken2 --db $KRAKEN_DB --confidence $conf --threads $threads --output $data_path/$barcode/$db/${barcode}_kraken2_output.txt --report $data_path/$barcode/$db/${barcode}_kraken2_report.txt $data_path/$barcode/$db/${barcode}_16S.fasta

            kraken-biom --max P --min $rank -o $data_path/$barcode/$db/${barcode}_kraken_biom.txt --fmt tsv $data_path/$barcode/$db/${barcode}_kraken2_report.txt

            conda activate taxonkit

            sed '1,2d' $data_path/$barcode/$db/${barcode}_kraken_biom.txt | taxonkit reformat --data-dir $TAXONKIT_DB --threads $threads --taxid-field 1 - | \
            sed 's/;/\t/g;s/[k,p,c,o,f,g,s]__//g' | awk 'BEGIN{FS="\t";OFS="\t"}{for (i=2;i<=NF;i++) gsub(/_/, " ", $i)} 1' >> $data_path/$barcode/$db/${barcode}_final_kraken2_result.txt

            rm -r $data_path/$barcode/$db/${barcode}_kraken_biom.txt $data_path/$barcode/$db/${barcode}_hist_temp.txt

        else
            echo "No Reads to Classify."
        fi
    fi
fi