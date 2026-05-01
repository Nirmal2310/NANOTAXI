#!/bin/bash

eval "$(conda shell.bash hook)"

helpFunction()
{
   echo "Usage: run_picrust2.sh -i /path/to/input/counts/file -f /path/to/consensus/fasta/file -t 16"
   echo -e "\t-i <str> Input table of sequence abundance."
   echo -e "\t-f <str> Consensus Fasta file."
   echo -e "\t-t <int> Number of threads."
   exit 1 # Exit script after printing help
}

threads=16

while getopts "i:f:t:" opt
do
    case "$opt" in
    i )
        input_counts="$OPTARG"
        ;;
    f )
        consensus_fasta="$OPTARG"
        ;;
    t )
	    threads="$OPTARG"
	    ;;
    ? ) 
        helpFunction
        ;;
    esac
done

if [ $(which conda | grep "condabin") ]; then

    path=$(which conda | sed 's/\/condabin.*$//g')

else

    path=$(which conda | sed 's/\/bin.*$//g')

fi

if [ -z "$consensus_fasta" ] || [ -z "$input_counts" ]
    then
    echo "Please provide the ASV Fasta file and Abundance file.";
    helpFunction
fi

input_path=$(dirname "$(readlink -f "$input_counts")")

output_dir=$(echo "$input_path/picrust2_out")

tmp_file=$(mktemp)

if [ ! -d $output_dir ]; then

    awk -F "\t" '{if(NR>1) print $0}' $input_counts | sort | uniq | cat <(head -n 1 $input_counts) - > $tmp_file && mv $tmp_file $input_counts
    
    conda activate bbtools

    awk -F "\t" '{if(NR>1) print $1}' $input_counts | seqkit grep -f - --threads $threads $consensus_fasta > $input_path/otu_representatives.fasta

    grep ">" $input_path/otu_representatives.fasta | sed 's/>//g' > $input_path/otu_ref_ids

    conda activate picrust2

    python $path/envs/picrust2/bin/place_seqs.py -s $input_path/otu_representatives.fasta -o $input_path/temp.tre -p $threads --intermediate $input_path/placement_workdir --min_align 0.8

    grep ">" $input_path/placement_workdir/study_seqs_hmmalign.fasta | sed 's/>//g' > $input_path/mapped_ref_ids

    rm -r $input_path/temp.tre $input_path/placement_workdir

    conda activate bbtools

    grep -vFf $input_path/mapped_ref_ids $input_path/otu_ref_ids | seqkit grep -f - --threads $threads $input_path/otu_representatives.fasta | \
    seqkit seq - -t DNA --threads $threads -r -p > $input_path/rc_representatives.fasta

    seqkit grep -f $input_path/mapped_ref_ids --threads $threads $input_path/otu_representatives.fasta | \
    cat - $input_path/rc_representatives.fasta > $input_path/temp && mv $input_path/temp $input_path/otu_representatives.fasta

    rm -r $input_path/mapped_ref_ids $input_path/otu_ref_ids $input_path/rc_representatives.fasta

    conda activate picrust2

    python $path/envs/picrust2/bin/picrust2_pipeline.py -s $input_path/otu_representatives.fasta -i $input_counts -o $output_dir -p $threads --remove_intermediate --min_align 0.8

    python $path/envs/picrust2/bin/add_descriptions.py -i $output_dir/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -o $output_dir/KO_metagenome_out/pred_metagenome_unstrat_annotated.tsv -m KO

    python $path/envs/picrust2/bin/add_descriptions.py -i $output_dir/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -o $output_dir/EC_metagenome_out/pred_metagenome_unstrat_annotated.tsv -m EC

    python $path/envs/picrust2/bin/add_descriptions.py -i $output_dir/pathways_out/path_abun_unstrat.tsv.gz -o $output_dir/pathways_out/path_abun_unstrat_annotated.tsv -m METACYC
else
    
    awk -F "\t" '{if(NR>1) print $0}' $input_counts | sort | uniq | cat <(head -n 1 $input_counts) - > $tmp_file && cp $tmp_file $input_counts && rm -r $tmp_file

    echo "Stopping since picrust2 output directory already exists"

fi