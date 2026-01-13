#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: real_time_emu.sh -d /path/to/data/directory -k kit-name -b barcode01 -m 1400 -M 1800 -q 10 -n REFSEQ -t 4 -s 500"
   echo -e "\t-d <str> Path Containing Sequencing Data."
   echo -e "\t-k <str> Kit-name."
   echo -e "\t-b <str> Barcode Name."
   echo -e "\t-m <int> Minimum Read Length. [default: 1400]"
   echo -e "\t-M <int> Maximum Read Length. [default: 1800]"
   echo -e "\t-q <int> Minimum Q-Score. [default: 10]"
   echo -e "\t-n <str> Database Name. [default: REFSEQ]"
   echo -e "\t-t <int> Number of Threads."
   echo -e "\t-s <int> Number of reads to process per barcode. [default: 500]"
   exit 1 # Exit script after printing help
}

min=1400
max=1800
q_score=10
db="REFSEQ"
batch_size=500

while getopts "d:k:b:m:M:q:n:t:s:" opt
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

if [ -z "$data_path" ]
    then
    echo "Please provide the path to the Data Directory.";
    helpFunction
fi

if [ "$db" == "REFSEQ" ]; then

    EMU_DB=$(grep EMU_REFSEQ ~/.bashrc | tail -n 1 | sed 's/export EMU_REFSEQ="//;s/"//g')

elif [ "$db" == "GTDB" ]; then
    
    EMU_DB=$(grep EMU_GTDB ~/.bashrc | tail -n 1 | sed 's/export EMU_GTDB="//;s/"//g')

elif [ "$db" == "MIMT" ]; then

    EMU_DB=$(grep EMU_MIMT ~/.bashrc | tail -n 1 | sed 's/export EMU_MIMT="//;s/"//g')

elif [ "$db" == "GSR" ]; then

    EMU_DB=$(grep EMU_GSR ~/.bashrc | tail -n 1 | sed 's/export EMU_GSR="//;s/"//g')

elif [ "$db" == "EMUDB" ]; then

    EMU_DB=$(grep EMU_DB ~/.bashrc | tail -n 1 | sed 's/export EMU_DB="//;s/"//g')

fi

if [ ! -d $data_path/$barcode/$db ]; then

    mkdir -p $data_path/$barcode/$db

fi

tmp_file=$(mktemp)


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

        conda activate emu

        emu abundance --type map-ont --db $EMU_DB --output-dir $data_path/$barcode/$db/ --output-basename $barcode --keep-counts --threads $threads ${data_path}/$barcode/$db/${barcode}_16S.fasta

        conda activate taxonkit

        line_count=$(wc -l ${data_path}/$barcode/$db/${barcode}_rel-abundance.tsv | awk '{print $1}')

        awk -v l=$line_count 'BEGIN{FS="\t";OFS="\t"}{if(NR>1 && NR<l) print $1,$10,$9,$8,$7,$6,$5,$4,$3; else if(NR>1 && NR==l) print "Unclassified",$10,"Unclassified","Unclassified","Unclassified","Unclassified","Unclassified","Unclassified","Unclassified"}' $data_path/$barcode/$db/${barcode}_rel-abundance.tsv > $data_path/$barcode/$db/${barcode}_final_emu_result.txt

        rm -r $data_path/$barcode/$db/${barcode}_hist_temp.txt $data_path/$barcode/$db/${barcode}_rel-abundance.tsv
    
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

            conda activate emu

            emu abundance --type map-ont --db $EMU_DB --output-dir $data_path/$barcode/$db/ --output-basename $barcode --keep-counts --threads $threads ${data_path}/$barcode/$db/${barcode}_16S.fasta

            conda activate taxonkit

            line_count=$(wc -l ${data_path}/$barcode/$db/${barcode}_rel-abundance.tsv | awk '{print $1}')

            awk -v l=$line_count 'BEGIN{FS="\t";OFS="\t"}{if(NR>1 && NR<l) print $1,$10,$9,$8,$7,$6,$5,$4,$3; else if(NR>1 && NR==l) print "Unclassified",$10,"Unclassified","Unclassified","Unclassified","Unclassified","Unclassified","Unclassified","Unclassified"}' $data_path/$barcode/$db/${barcode}_rel-abundance.tsv >> $data_path/$barcode/$db/${barcode}_final_emu_result.txt

            rm -r $data_path/$barcode/$db/${barcode}_hist_temp.txt $data_path/$barcode/$db/${barcode}_rel-abundance.tsv
        
        else
            
            echo "No Reads to Classify."
        
        fi
    fi
fi