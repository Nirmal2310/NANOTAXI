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

if { conda env list |  grep "blast"; } > /dev/null 2>&1; then

        conda list -n blast --explicit > _current_env.txt

        if diff -q $script_dir/blast.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name blast --file $script_dir/blast.txt -y && rm -r _current_env.txt
		conda list -n blast --explicit > $script_dir/blast.txt
        fi

else

        conda create --name blast --file $script_dir/blast.txt
	conda list -n blast --explicit > $script_dir/blast.txt

fi

if { conda env list | grep "nanofilt";} > /dev/null 2>&1; then

        conda list -n nanofilt --explicit > _current_env.txt

        if diff -q $script_dir/nanofilt.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name nanofilt --file $script_dir/nanofilt.txt -y && rm -r _current_env.txt
		conda list -n nanofilt --explicit > $script_dir/nanofilt.txt
        fi

else
        
        conda create -n nanofilt --file $script_dir/nanofilt.txt
	conda list -n nanofilt --explicit > $script_dir/nanofilt.txt
        
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

cd DATA

if [ ! -d TAXONKIT_DATA ]; then
                
        mkdir TAXONKIT_DATA

        cd TAXONKIT_DATA

        wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 

        tar -zxvf taxdump.tar.gz

        rm -r taxdump.tar.gz

        grep -qF "export TAXONKIT_DB=\"$PWD\"" ~/.bashrc || echo "export TAXONKIT_DB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        cd $base_dir

fi

cd $base_dir/DATA/TAXONKIT_DATA

grep -qF "export TAXONKIT_DB=\"$PWD\"" ~/.bashrc || echo "export TAXONKIT_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir


if [ ! -d DATA ]; then
        
	mkdir DATA
fi

cd DATA

if [ ! -d BLAST ]; then
        
        mkdir BLAST
fi

cd BLAST

if [ ! -d REFSEQ ]; then
                
	mkdir REFSEQ

        cd REFSEQ
        
        wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz

        zcat bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz > refseq_16S.fasta && rm -r bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz

        threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)

        source $path/bin/activate taxonkit

        grep ">" refseq_16S.fasta | sed 's/>//g' | awk -F " " '{print $1"\t"$2,$3}' | taxonkit name2taxid --threads $threads --data-dir $TAXONKIT_DB -i 2 | \
        awk -F "\t" '!seen[$1]++' | taxonkit reformat2 --data-dir $TAXONKIT_DB --threads $threads -I 3 | sed 's/;/\t/g' | \
        awk 'BEGIN{FS=OFS="\t"}{if(NR==1) print "REF_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"; else print $1,$4,$5,$6,$7,$8,$9,$10}' > RefSeq_taxa.txt

        source $path/bin/activate seqkit

        seqkit faidx refseq_16S.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' refseq_16S.fasta.fai > refseq_filtered_ids

        seqkit faidx -X refseq_filtered_ids refseq_16S.fasta > refseq_final_seqs.fasta

        source $path/bin/activate blast

        makeblastdb -in refseq_final_seqs.fasta -parse_seqids -blastdb_version 5 -title REFSEQ_BLAST -dbtype nucl -out REFSEQ_BLAST

        rm -r refseq_filtered_ids refseq_16S.fasta* refseq_final_seqs.fasta

        grep -qF "export BLAST_REFSEQ=\"$PWD\"" ~/.bashrc || echo "export BLAST_REFSEQ=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc
fi
        
cd $base_dir/DATA/BLAST/REFSEQ

grep -qF "export BLAST_REFSEQ=\"$PWD\"" ~/.bashrc || echo "export BLAST_REFSEQ=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

source $path/bin/activate base
        
cd $base_dir/DATA/BLAST

if [ ! -d MIMT ]; then

        mkdir MIMT

        cd MIMT

        wget -c https://people.biopolis.pt/bu/mimt/downloads/16S_files/MIMt-16S_M2c_25_10_taxid.fna.gz -O MIMt.fasta.gz && gunzip MIMt.fasta.gz

        sed -i 's/rrna_//g' MIMt.fasta

        source $base/bin/activate taxonkit

        threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)

        grep ">" MIMt.fasta | sed 's/>//g' | taxonkit reformat2 --data-dir $TAXONKIT_DB --threads $threads -I 2 | sed 's/;/\t/g' | \
        awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1) print "REF_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"; else print $1,$3,$4,$5,$6,$7,$8,$9}' > MIMT_taxa.txt

        source $path/bin/activate seqkit

        seqkit faidx MIMt.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' MIMt.fasta.fai > mimt_filtered_ids

        seqkit faidx -X mimt_filtered_ids MIMt.fasta > MIMT_final_seqs.fasta

        source $path/bin/activate blast

        makeblastdb -in MIMT_final_seqs.fasta -parse_seqids -blastdb_version 5 -title MIMT_BLAST -dbtype nucl -out MIMT_BLAST

        rm -r mimt_filtered_ids MIMt.fasta*

        grep -qF "export BLAST_MIMT=\"$PWD\"" ~/.bashrc || echo "export BLAST_MIMT=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/BLAST/MIMT

grep -qF "export BLAST_MIMT=\"$PWD\"" ~/.bashrc || echo "export BLAST_MIMT=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

source $path/bin/activate base
        
cd $base_dir/DATA/BLAST

if [ ! -d GTDB ]; then

        mkdir GTDB && cd GTDB

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/bac120_ssu_reps.fna.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/ar53_ssu_reps.fna.gz

        zcat bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz > GTDB_16S_reps.fasta && rm -r bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz

        grep ">" GTDB_16S_reps.fasta | sed 's/>//g' | sed 's/ /\t/;s/;/\t/g;s/[d,p,c,o,f,g,s]__//g' | sed 's/ \[locus.*$//g' | \
        awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1) print "REF_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"; else print $0}' > GTDB_taxa.txt

        source $path/bin/activate seqkit

        seqkit faidx GTDB_16S_reps.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' GTDB_16S_reps.fasta.fai > gtdb_filtered_ids

        seqkit faidx -X gtdb_filtered_ids GTDB_16S_reps.fasta > GTBD_final_seqs.fasta

        source $path/bin/activate blast

        makeblastdb -in GTBD_final_seqs.fasta -parse_seqids -blastdb_version 5 -title GTDB_BLAST -dbtype nucl -out GTDB_BLAST

        rm -r gtdb_filtered_ids GTDB_16S_reps.fasta* GTBD_final_seqs.fasta*

        grep -qF "export BLAST_GTDB=\"$PWD\"" ~/.bashrc || echo "export BLAST_GTDB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/BLAST/GTDB

grep -qF "export BLAST_GTDB=\"$PWD\"" ~/.bashrc || echo "export BLAST_GTDB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

source $path/bin/activate base
        
cd $base_dir/DATA/BLAST

if [ ! -d GSR ]; then

        mkdir GSR && cd GSR

        wget -c https://manichanh.vhir.org/gsrdb/GSR-DB_full-16S.tar.gz

        tar -xvf GSR-DB_full-16S.tar.gz && rm -r GSR-DB_full-16S.tar.gz GSR-DB_full-16S_filt_taxa.qza GSR-DB_full-16S_filt_seqs.qza

        sed 's/ //g;s/;/\t/g;s/[k,p,c,o,f,g,s]__//g;s/_/ /g' GSR-DB_full-16S_filt_taxa.txt | \
        awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1) print "REF_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"; else print $0}' > GSR_taxa.txt

        source $path/bin/activate seqkit

        seqkit faidx GSR-DB_full-16S_filt_seqs.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' GSR-DB_full-16S_filt_seqs.fasta.fai > gsr_filtered_ids

        seqkit faidx -X gsr_filtered_ids GSR-DB_full-16S_filt_seqs.fasta > GSR_final_seqs.fasta

        source $path/bin/activate blast

        makeblastdb -in GSR_final_seqs.fasta -parse_seqids -blastdb_version 5 -title GSR_BLAST -dbtype nucl -out GSR_BLAST

        rm -r gsr_filtered_ids GSR-DB_full-16S_filt_seqs.fasta* GSR_final_seqs.fasta*

        grep -qF "export BLAST_GSR=\"$PWD\"" ~/.bashrc || echo "export BLAST_GSR=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/BLAST/GSR

grep -qF "export BLAST_GSR=\"$PWD\"" ~/.bashrc || echo "export BLAST_GSR=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

source $path/bin/activate base
        
cd $base_dir