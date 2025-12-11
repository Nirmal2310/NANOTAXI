#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: real_time_analysis.sh -d /path/to/data/directory -k kit-name -b barcode01 -m 1400 -M 1800 -t Species -c 0.0 -q 10 -n REFSEQ"
   echo -e "\t-d <str> Path Containing Sequencing Data."
   echo -e "\t-k <str> Kit-name."
   echo -e "\t-b <str> Barcode Name."
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-r <str> Minimum Taxonomy Rank. [default: Species]"
   echo -e "\t-c <int> Confidence Score. [default: 0.0]"
   echo -e "\t-q <int> Minimum Q-Score. [default: 10]"
   echo -e "\t-n <str> Database Name. [default: REFSEQ]"
   exit 1 # Exit script after printing help
}

min=1400
max=1800
rank=S
conf=0.0
q_score=10
db="REFSEQ"

while getopts "d:k:b:m:M:t:c:q:n" opt
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
    ? ) helpFunction ;;
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

elif [ "$db" == "EMUDB" ]; then

    KRAKEN_DB=$(grep KRAKEN_EMUDB ~/.bashrc | tail -n 1 | sed 's/export KRAKEN_EMUDB="//;s/"//g')    

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

        conda activate kraken2

        kraken2 --db $KRAKEN_DB --confidence $conf --output $data_path/$barcode/${barcode}_kraken2_output.txt --report $data_path/$barcode/${barcode}_kraken2_report.txt $data_path/$barcode/${barcode}_16S.fasta

        kraken-biom --max D --min $rank -o $data_path/$barcode/${barcode}_kraken_biom.txt --fmt tsv $data_path/$barcode/${barcode}_kraken2_report.txt

        conda activate taxonkit

        sed '1,2d' $data_path/$barcode/${barcode}_kraken_biom.txt | taxonkit reformat --data-dir $TAXONKIT_DB --taxid-field 1 - | sed 's/;/\t/g' > $data_path/$barcode/$db/${barcode}_final_kraken2_result.txt

        rm -r $data_path/$barcode/${barcode}_kraken_biom.txt $data_path/$barcode/${barcode}_hist_temp.txt
        
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

        ls $data_path/$barcode/*fastq.gz | grep -vf $data_path/$barcode/processed_files.txt | xargs zcat | bioawk -c fastx '{print meanqual($qual)}' >> $data_path/$barcode/${barcode}_quality.txt

        if [ "$(grep ">" $data_path/$barcode/${barcode}_16S.fasta | wc -l)" -gt 0 ]; then

            conda activate kraken2

            kraken2 --db $KRAKEN_DB --output $data_path/$barcode/${barcode}_kraken2_output.txt --report $data_path/$barcode/${barcode}_kraken2_report.txt $data_path/$barcode/${barcode}_16S.fasta

            kraken-biom --max D --min $rank -o $data_path/$barcode/${barcode}_kraken_biom.txt --fmt tsv $data_path/$barcode/${barcode}_kraken2_report.txt

            conda activate taxonkit

            sed '1,2d' $data_path/$barcode/${barcode}_kraken_biom.txt | taxonkit reformat --data-dir $TAXONKIT_DB --taxid-field 1 - | sed 's/;/\t/g' >> $data_path/$barcode/$db/${barcode}_final_kraken2_result.txt

            rm -r $data_path/$barcode/${barcode}_kraken_biom.txt $data_path/$barcode/${barcode}_hist_temp.txt

            ls $data_path/$barcode/*fastq.gz | grep -vf $data_path/$barcode/processed_files.txt >> $data_path/$barcode/processed_files.txt
        else
            echo "No Reads to Classify."
        fi
    fi
fi