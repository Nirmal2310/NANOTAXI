#!/bin/bash
if which conda >/dev/null; then
        
        echo "Conda Exist"

else
        source ~/.bashrc
        
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
        && chmod +x miniconda.sh && bash miniconda.sh -b -p miniconda
        
        base_dir=$(echo $PWD)
        
        export PATH=$base_dir/miniconda/bin:$PATH
        
        source ~/.bashrc
        
        echo -e "$base_dir/miniconda/etc/profile.d/conda.sh" >> ~/.profile
        
        conda init bash

fi

if [ $(which conda | grep "condabin") ]; then

    path=$(which conda | sed 's/\/condabin.*$//g')

else

    path=$(which conda | sed 's/\/bin.*$//g')

fi

base_dir=$PWD

source ~/.bashrc

script_dir=$(dirname "$(readlink -f "$0")")

if { conda env list |  grep "mmseqs"; } > /dev/null 2>&1; then

        conda list -n mmseqs --explicit > _current_env.txt

        if diff -q $script_dir/mmseqs.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name mmseqs --file $script_dir/mmseqs.txt -y && rm -r _current_env.txt
		conda list -n mmseqs --explicit > $script_dir/mmseqs.txt
        fi

else

        conda create --name mmseqs --file $script_dir/mmseqs.txt
	conda list -n mmseqs --explicit > $script_dir/mmseqs.txt

fi

if { conda env list | grep "chopper";} > /dev/null 2>&1; then

        conda list -n chopper --explicit > _current_env.txt

        if diff -q $script_dir/chopper.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name chopper --file $script_dir/chopper.txt -y && rm -r _current_env.txt
		conda list -n chopper --explicit > $script_dir/chopper.txt
        fi

else
        
        conda create -n chopper --file $script_dir/chopper.txt
	conda list -n chopper --explicit > $script_dir/chopper.txt
        
fi

if { conda env list | grep "taxonkit";} > /dev/null 2>&1; then

        conda list -n taxonkit --explicit > _current_env.txt

        if diff -q $script_dir/taxonkit.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name taxonkit --file $script_dir/taxonkit.txt -y && rm -r _current_env.txt
		conda list -n taxonkit --explicit > $script_dir/taxonkit.txt
        fi

else
        
        conda create -n taxonkit --file $script_dir/taxonkit.txt
	conda list -n taxonkit --explicit > $script_dir/taxonkit.txt
        
fi

if [ ! -d DATA ]; then
        
	mkdir DATA
fi

cd DATA

if [ ! -d TAXONKIT_DATA ]; then
                
        mkdir TAXONKIT_DATA

        cd TAXONKIT_DATA

        wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 

        tar -zxvf taxdump.tar.gz

        rm -r taxdump.tar.gz

        grep -qF "export TAXONKIT_DB=\"$PWD\"" ~/.bashrc || echo "export TAXONKIT_DB=\"$PWD\"" >> ~/.bashrc

        wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

        wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz

        zcat bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz > refseq_16S.fasta && rm -r bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz

        zcat nucl_gb.accession2taxid.gz | grep -w -f <(grep -oP '^>[^\s]+' refseq_16S.fasta | sed 's/^>//') | \
        awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3}' > refseq_taxid.txt

        rm -r nucl_gb.accession2taxid.gz refseq_16S.fasta

        source ~/.bashrc

        cd $base_dir

fi

cd $base_dir/DATA/TAXONKIT_DATA

grep -qF "export TAXONKIT_DB=\"$PWD\"" ~/.bashrc || echo "export TAXONKIT_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA

if [ ! -d MMSEQS ]; then
        
        mkdir MMSEQS
fi

cd MMSEQS

if [ ! -d REFSEQ ]; then
                
	mkdir REFSEQ

        cd REFSEQ
        
        zcat $base_dir/DATA/TMP_DIR/bacteria.16SrRNA.fna.gz $base_dir/DATA/TMP_DIR/archaea.16SrRNA.fna.gz > refseq_16S.fasta

        threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)

        source $path/bin/activate taxonkit

        taxonkit reformat2 --data-dir $TAXONKIT_DB --threads $threads -f "{domain};{phylum};{class};{order};{family};{genus};{species}" -I 2 $TAXONKIT_DB/refseq_taxid.txt | \
        awk 'BEGIN{FS=OFS="\t"}{print $1,$3}' | sort | uniq | sed 's/;/\t/g' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tfamily\tGenus\tSpecies") - > RefSeq_taxa.txt

        source $path/bin/activate bbtools

        seqkit faidx refseq_16S.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=1200 && $2<=1800) print $1}' refseq_16S.fasta.fai > refseq_filtered_ids

        seqkit faidx -X refseq_filtered_ids refseq_16S.fasta > refseq_final_seqs.fasta

        source $path/bin/activate mmseqs

        mmseqs createdb refseq_final_seqs.fasta REFSEQ_MMSEQS -v 0 --threads $threads

        mmseqs createindex -v 0 --threads $threads --remove-tmp-files 1 REFSEQ_MMSEQS tmp --search-type 3

        rm -r refseq_filtered_ids refseq_16S.fasta* refseq_final_seqs.fasta seqid_taxid.txt

        grep -qF "export MMSEQS_REFSEQ=\"$PWD\"" ~/.bashrc || echo "export MMSEQS_REFSEQ=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc
fi
        
cd $base_dir/DATA/MMSEQS/REFSEQ

grep -qF "export MMSEQS_REFSEQ=\"$PWD\"" ~/.bashrc || echo "export MMSEQS_REFSEQ=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

source $path/bin/activate base
        
cd $base_dir/DATA/MMSEQS

if [ ! -d MIMT ]; then

        mkdir MIMT

        cd MIMT

        zcat $base_dir/DATA/TMP_DIR/MIMt.fasta.gz >  MIMt.fasta

        sed -i 's/rrna_//g' MIMt.fasta

        zcat $base_dir/DATA/TMP_DIR/MIMT_taxa.txt.gz > MIMT_taxa.txt

        sed -i 's/rrna_//' MIMT_taxa.txt

        sed 's/;/\t/g;s/[K,P,C,O,F,G,S]__//g' MIMT_taxa.txt | awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) for (i=2;i<=NF;i++) gsub(/_/, " ", $i) split($8, a, " "); $8=a[1]" "a[2]} 1' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > temp && mv temp MIMT_taxa.txt

        source $path/bin/activate bbtools

        seqkit faidx MIMt.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=1200 && $2<=1800) print $1}' MIMt.fasta.fai > mimt_filtered_ids

        seqkit faidx -X mimt_filtered_ids MIMt.fasta > MIMT_final_seqs.fasta

        source $path/bin/activate mmseqs

        mmseqs createdb MIMT_final_seqs.fasta MIMT_MMSEQS -v 0 --threads $threads

        mmseqs createindex -v 0 --threads $threads --remove-tmp-files 1 MIMT_MMSEQS tmp --search-type 3

        rm -r mimt_filtered_ids MIMt.fasta* MIMT_final_seqs.fasta*

        grep -qF "export MMSEQS_MIMT=\"$PWD\"" ~/.bashrc || echo "export MMSEQS_MIMT=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/MMSEQS/MIMT

grep -qF "export MMSEQS_MIMT=\"$PWD\"" ~/.bashrc || echo "export MMSEQS_MIMT=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

source $path/bin/activate base
        
cd $base_dir/DATA/MMSEQS

if [ ! -d GTDB ]; then

        mkdir GTDB && cd GTDB

        zcat $base_dir/DATA/TMP_DIR/bac120_ssu_reps.fna.gz $base_dir/DATA/TMP_DIR/ar53_ssu_reps.fna.gz > GTDB_16S_reps.fasta

        zcat $base_dir/DATA/TMP_DIR/bac120_metadata.tsv.gz $base_dir/DATA/TMP_DIR/ar53_metadata.tsv.gz | grep -v "ncbi" | awk -F "\t" '{print $1"\t"$82}' | sed 's/;/\t/g;s/[d,p,c,o,f,g,s]__//g' | \
        awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8}' - | sort -k1 -n -r | uniq | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > GTDB_taxa.txt

        source $path/bin/activate bbtools

        seqkit faidx GTDB_16S_reps.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=1200 && $2<=1800) print $1}' GTDB_16S_reps.fasta.fai > gtdb_filtered_ids

        seqkit faidx -X gtdb_filtered_ids GTDB_16S_reps.fasta > GTBD_final_seqs.fasta

        source $path/bin/activate mmseqs

        mmseqs createdb GTBD_final_seqs.fasta GTDB_MMSEQS -v 0 --threads $threads

        mmseqs createindex -v 0 --threads $threads --remove-tmp-files 1 GTDB_MMSEQS tmp --search-type 3

        rm -r gtdb_filtered_ids GTDB_16S_reps.fasta* GTBD_final_seqs.fasta*

        grep -qF "export MMSEQS_GTDB=\"$PWD\"" ~/.bashrc || echo "export MMSEQS_GTDB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/MMSEQS/GTDB

grep -qF "export MMSEQS_GTDB=\"$PWD\"" ~/.bashrc || echo "export MMSEQS_GTDB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

source $path/bin/activate base
        
cd $base_dir/DATA/MMSEQS

if [ ! -d GSR ]; then

        mkdir GSR && cd GSR

        cp $base_dir/DATA/TMP_DIR/GSR-DB_full-16S_filt_seqs.fasta GSR-DB_full-16S_filt_seqs.fasta 

        awk '{if(NR>1) print $0}' $base_dir/DATA/TMP_DIR/GSR-DB_full-16S_filt_taxa.txt | sed 's/ //g;s/;/\t/g;s/[k,p,c,o,f,g,s]__//g' | \
        awk 'BEGIN{FS="\t";OFS="\t"}{for (i=2;i<=NF;i++) gsub(/_/, " ", $i)} 1' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > GSR_taxa.txt

        source $path/bin/activate bbtools

        seqkit faidx GSR-DB_full-16S_filt_seqs.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=1200 && $2<=1800) print $1}' GSR-DB_full-16S_filt_seqs.fasta.fai > gsr_filtered_ids

        seqkit faidx -X gsr_filtered_ids GSR-DB_full-16S_filt_seqs.fasta > GSR_final_seqs.fasta

        source $path/bin/activate mmseqs

        mmseqs createdb GSR_final_seqs.fasta GSR_MMSEQS -v 0 --threads $threads

        mmseqs createindex -v 0 --threads $threads --remove-tmp-files 1 GSR_MMSEQS tmp --search-type 3

        rm -r gsr_filtered_ids GSR-DB_full-16S_filt_seqs.fasta* GSR_final_seqs.fasta*

        grep -qF "export MMSEQS_GSR=\"$PWD\"" ~/.bashrc || echo "export MMSEQS_GSR=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/MMSEQS/GSR

grep -qF "export MMSEQS_GSR=\"$PWD\"" ~/.bashrc || echo "export MMSEQS_GSR=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

source $path/bin/activate base
        
cd $base_dir/DATA/MMSEQS

if [ ! -d EMUDB ]; then

        mkdir EMUDB && cd EMUDB

        sed 's/ .*$//g;s/:/_/g' $base_dir/DATA/TMP_DIR/species_taxid.fasta > species_taxid.fasta

        source $path/bin/activate bbtools

        seqkit faidx species_taxid.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' species_taxid.fasta.fai > EMU_filtered_ids

        seqkit faidx -X EMU_filtered_ids species_taxid.fasta > temp && mv temp species_taxid.fasta

        grep ">" species_taxid.fasta | sed 's/>//;s/ .*$//g' | awk 'BEGIN{FS=OFS="\t"}{$2=$1; gsub(/_.*$/,"",$2); print $2,$1}' | sort -k 1b,1 > temp

        awk 'BEGIN{FS=OFS="\t"}{if(NR>1) print $1,$9,$7,$6,$5,$4,$3,$2}' $base_dir/DATA/TMP_DIR/taxonomy.tsv | sort -k 1b,1 > temp2

        join -t $'\t' temp temp2 | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9}' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > EMU_taxa.txt

        source $path/bin/activate mmseqs

        mmseqs createdb species_taxid.fasta EMU_MMSEQS -v 0 --threads $threads

        mmseqs createindex -v 0 --threads $threads --remove-tmp-files 1 EMU_MMSEQS tmp --search-type 3

        rm -r species_taxid.fasta* temp* EMU_filtered_ids

        grep -qF "export MMSEQS_EMU=\"$PWD\"" ~/.bashrc || echo "export MMSEQS_EMU=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/MMSEQS/EMUDB

grep -qF "export MMSEQS_EMU=\"$PWD\"" ~/.bashrc || echo "export MMSEQS_EMU=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc
        
cd $base_dir