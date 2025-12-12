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

        cp $TAXONKIT_DB/refseq_taxid.txt seqid_taxid.txt

        taxonkit reformat2 --data-dir $TAXONKIT_DB --threads $threads -f "{domain};{phylum};{class};{order};{family};{genus};{species}" -I 2 seqid_taxid.txt | \
        awk 'BEGIN{FS=OFS="\t"}{print $1,$3}' | sort | uniq | sed 's/;/\t/g' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tfamily\tGenus\tSpecies") - > RefSeq_taxa.txt

        source $path/bin/activate seqkit

        seqkit faidx refseq_16S.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' refseq_16S.fasta.fai > refseq_filtered_ids

        seqkit faidx -X refseq_filtered_ids refseq_16S.fasta > refseq_final_seqs.fasta

        source $path/bin/activate blast

        makeblastdb -in refseq_final_seqs.fasta -parse_seqids -blastdb_version 5 -title REFSEQ_BLAST -dbtype nucl -out REFSEQ_BLAST

        rm -r refseq_filtered_ids refseq_16S.fasta* refseq_final_seqs.fasta seqid_taxid.txt

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

        wget -c https://people.biopolis.pt/bu/mimt/downloads/16S_files/MIMt-16S_M2c_25_10.tax.gz -O MIMT_taxa.txt.gz && gunzip MIMT_taxa.txt.gz

        sed -i 's/rrna_//' MIMT_taxa.txt

        sed 's/;/\t/g;s/[K,P,C,O,F,G,S]__//g' MIMT_taxa.txt | awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) for (i=2;i<=NF;i++) gsub(/_/, " ", $i)} 1' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > temp && mv temp MIMT_taxa.txt

        source $path/bin/activate seqkit

        seqkit faidx MIMt.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' MIMt.fasta.fai > mimt_filtered_ids

        seqkit faidx -X mimt_filtered_ids MIMt.fasta > MIMT_final_seqs.fasta

        source $path/bin/activate blast

        makeblastdb -in MIMT_final_seqs.fasta -parse_seqids -blastdb_version 5 -title MIMT_BLAST -dbtype nucl -out MIMT_BLAST

        rm -r mimt_filtered_ids MIMt.fasta* MIMT_final_seqs.fasta*

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
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > GTDB_taxa.txt

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

        sed 's/ //g;s/;/\t/g;s/[k,p,c,o,f,g,s]__//g' GSR-DB_full-16S_filt_taxa.txt | \
        awk 'BEGIN{FS="\t";OFS="\t"}{for (i=2;i<=NF;i++) gsub(/_/, " ", $i)} 1' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > GSR_taxa.txt

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
        
cd $base_dir/DATA/BLAST

if [ ! -d EMUDB ]; then

        source $path/bin/activate emu

        mkdir EMUDB

        osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar.gz

        tar -xvf emu.tar.gz -C EMUDB --strip 1

        rm -r emu.tar.gz

        cd EMUDB

        sed -i 's/ .*$//g' species_taxid.fasta

        join <(grep ">" species_taxid.fasta | sed 's/>//;s/ .*$//g' | awk 'BEGIN{FS=OFS="\t"}{$2=$1; gsub(/:.*$/,"",$2); print $2,$1}' | sort -k1 -n -r) \
        <(awk 'BEGIN{FS=OFS="\t"}{if(NR>1) print $1,$9,$7,$6,$5,$4,$3,$2}' taxonomy.tsv | sort -k1 -n -r) | \
        awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9" "$10}' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > EMU_taxa.txt

        source $path/bin/activate blast

        makeblastdb -in species_taxid.fasta -parse_seqids -blastdb_version 5 -title EMU_BLAST -dbtype nucl -out EMU_BLAST

        rm -r species_taxid.fasta taxonomy.tsv

        grep -qF "export BLAST_EMU=\"$PWD\"" ~/.bashrc || echo "export BLAST_EMU=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/BLAST/EMUDB

grep -qF "export BLAST_EMU=\"$PWD\"" ~/.bashrc || echo "export BLAST_EMU=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc
        
cd $base_dir