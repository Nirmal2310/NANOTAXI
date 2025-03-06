#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: real_time_analysis_blast.sh -d /path/to/data/directory -b barcode01 -m 1400 -M 1800 -i 85"
   echo -e "\t-d <str> Path Containing Sequencing Data."
   echo -e "\t-b <str> Barcode Name."
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-i <int> Minimum Percent Identity. [default: 85]"
   echo -e "\t-c <int> Minimum Percent Coverage. [default: 85]"
   exit 1 # Exit script after printing help
}

while getopts "d:b:m:M:i:c:" opt
do
    case "$opt" in
    d )
        data_path="$OPTARG"
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
    ? ) helpFunction ;;
    esac
done

min=1400
max=1800
identity=85
coverage=85

script_path=$(dirname "$(readlink -f "$0")")

if [ -z "$data_path" ]
    then
    echo "Please provide the path to the Data Directory.";
    helpFunction
fi

GSR_DB=$(grep GSR_DB ~/.bashrc | tail -n 1 | sed 's/export GSR_DB="//;s/"//g')

TAXONKIT_DB=$(grep TAXONKIT_DB ~/.bashrc | tail -n 1 | sed 's/export TAXONKIT_DB="//;s/"//g')


if [ ! -f $data_path/$barcode/processed_files.txt ]; then
    
    conda activate nanofilt
    
    ls $data_path/$barcode/*fastq.gz | xargs zcat | NanoFilt -q 10 --length $min --maxlength $max | sed -n '1~4s/^@/>/p;2~4p' > ${data_path}/$barcode/${barcode}_16s.fasta

    conda activate bbtools

    ls $data_path/$barcode/*fastq.gz | xargs zcat | readlength.sh in=stdin.fq out=$data_path/$barcode/${barcode}_hist_temp.txt bin=10 round=t -ignorebadquality

    paste -d "\t" <(echo $barcode) <(grep "#Avg" $data_path/$barcode/${barcode}_hist_temp.txt | cut -f2) > $data_path/$barcode/${barcode}_average_length.txt

    grep -v "#" $data_path/$barcode/${barcode}_hist_temp.txt | awk '{print $1"\t"$2}' > $data_path/$barcode/${barcode}_hist.txt

    ls $data_path/$barcode/*fastq.gz | xargs zcat | bioawk -c fastx '{print meanqual($qual)}' > $data_path/$barcode/${barcode}_quality.txt

    if [ "$(grep ">" $data_path/$barcode/${barcode}_16s.fasta | wc -l)" -gt 0 ]; then

        conda activate minimap2

        minimap2 -ax map-ont -t 2 --eqx $GSR_DB/GSR-DB_full-16s_filt_seqs.fasta $data_path/$barcode/${barcode}_16s.fasta | \
	 samtools view -@ 1 -F 2034 -bS | samtools sort -@ 1 -o $data_path/$barcode/${barcode}_16s.bam; samtools index -@ 1 $data_path/$barcode/${barcode}_16s.bam

        python $script_path/alignment_filter.py -b $data_path/$barcode/${barcode}_16s.bam -t $GSR_DB/GSR-DB_full-16s_filt_taxa.txt -i $identity -c $coverage | \
        awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) print $0}' > $data_path/$barcode/${barcode}_16s.temp

        conda activate taxonkit

        taxonkit name2taxid -i 1 -j 2 --data-dir $TAXONKIT_DB $data_path/$barcode/${barcode}_16s.temp | awk 'BEGIN{FS="\t";OFS="\t"}{print $3, $2}' | \
            taxonkit reformat --threads 2 --data-dir $TAXONKIT_DB --taxid-field 1 - | sed 's/;/\t/g' > $data_path/$barcode/${barcode}_final_minimap2_result.txt

        rm -r $data_path/$barcode/${barcode}_hist_temp.txt $data_path/$barcode/${barcode}_16s.temp $data_path/$barcode/${barcode}_16s.bam
        
        ls $data_path/$barcode/*fastq.gz > $data_path/$barcode/processed_files.txt
    
    else
        
        echo "No Reads to Classify."
    
    fi

else

    if [ "$(echo $data_path/$barcode/*fastq.gz | tr ' ' '\n' | grep -vf $data_path/$barcode/processed_files.txt | wc -l)" -eq 0 ]; then
        
        echo "No New Files Found."

    else

        conda activate nanofilt

        ls $data_path/$barcode/*fastq.gz | grep -vf $data_path/$barcode/processed_files.txt | xargs zcat | NanoFilt -q 10 --length $min --maxlength $max | sed -n '1~4s/^@/>/p;2~4p' > $data_path/$barcode/${barcode}_16s.fasta

        conda activate bbtools

        ls $data_path/$barcode/*fastq.gz | grep -vf $data_path/$barcode/processed_files.txt | xargs zcat | readlength.sh in=stdin.fq out=$data_path/$barcode/${barcode}_hist_temp.txt bin=10 round=t -ignorebadquality

        paste -d "\t" <(echo $barcode) <(grep "#Avg" $data_path/$barcode/${barcode}_hist_temp.txt | cut -f2) > $data_path/$barcode/${barcode}_average_length.txt

        grep -v "#" $data_path/$barcode/${barcode}_hist_temp.txt | awk '{print $1"\t"$2}' >> $data_path/$barcode/${barcode}_hist.txt

        ls $data_path/$barcode/*fastq.gz | grep -vf $data_path/$barcode/processed_files.txt | xargs zcat | bioawk -c fastx '{print meanqual($qual)}' >> $data_path/$barcode/${barcode}_quality.txt

        if [ "$(grep ">" $data_path/$barcode/${barcode}_16s.fasta | wc -l)" -gt 0 ]; then

            conda activate minimap2

            minimap2 -ax map-ont -t 2 --eqx $GSR_DB/GSR-DB_full-16s_filt_seqs.fasta $data_path/$barcode/${barcode}_16s.fasta | samtools view -@ 1 -F 2034 -bS | samtools sort -@ 1 -o $data_path/$barcode/${barcode}_16s.bam

            python $script_path/alignment_filter.py -b $data_path/$barcode/${barcode}_16s.bam -t $GSR_DB/GSR-DB_full-16s_filt_taxa.txt -i $identity -c $coverage | \
                awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) print $0}' > $data_path/$barcode/${barcode}_16s.temp

            conda activate taxonkit

            taxonkit name2taxid -i 1 -j 2 --data-dir $TAXONKIT_DB $data_path/$barcode/${barcode}_16s.temp | awk 'BEGIN{FS="\t";OFS="\t"}{print $3, $2}' | \
                taxonkit reformat --threads 2 --data-dir $TAXONKIT_DB --taxid-field 1 - | sed 's/;/\t/g' > $data_path/$barcode/${barcode}_final_minimap2_result.txt
            
            rm -r $data_path/$barcode/${barcode}_hist_temp.txt $data_path/$barcode/${barcode}_16s.temp $data_path/$barcode/${barcode}_16s.bam

            ls $data_path/$barcode/*fastq.gz | grep -vf $data_path/$barcode/processed_files.txt >> $data_path/$barcode/processed_files.txt
        
        else
            
            echo "No Reads to Classify."
        
        fi
    fi
fi
