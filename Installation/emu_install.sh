#!/bin/env bash

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

if { conda env list |  grep "emu"; } > /dev/null 2>&1; then

        conda list -n emu --explicit > _current_env.txt

        if diff -q $script_dir/emu.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name emu --file $script_dir/emu.txt -y && rm -r _current_env.txt
	
  		source $path/bin/activate emu
	
	        pip install osfclient
	
	        source $path/bin/activate base

  		conda list -n emu --explicit > $script_dir/emu.txt
  		
        fi

else

        conda create --name emu --file $script_dir/emu.txt

        source $path/bin/activate emu

        pip install osfclient

        source $path/bin/activate base

 	conda list -n emu --explicit > $script_dir/emu.txt

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

fi

cd $base_dir/DATA/TAXONKIT_DATA

grep -qF "export TAXONKIT_DB=\"$PWD\"" ~/.bashrc || echo "export TAXONKIT_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA

if [ ! -d EMU ]; then
        
        mkdir EMU
fi

cd EMU

if [ ! -d EMU_DATA ]; then
        
        source $path/bin/activate emu
        
        mkdir EMU_DATA

        osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar.gz

        tar -xvf emu.tar.gz -C EMU_DATA --strip 1

        rm -r emu.tar.gz

        cd EMU_DATA

        awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9}' taxonomy.tsv > temp && mv temp taxonomy.tsv
        
        source $path/bin/activate base

        grep -qF "export EMU_DB=\"$PWD\"" ~/.bashrc || echo "export EMU_DB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/EMU/EMU_DATA

grep -qF "export EMU_DB=\"$PWD\"" ~/.bashrc || echo "export EMU_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/EMU

if [ ! -d GSR ]; then
        
        wget -c https://manichanh.vhir.org/gsrdb/GSR-DB_full-16S.tar.gz

        tar -xvf GSR-DB_full-16S.tar.gz && rm -r GSR-DB_full-16S.tar.gz GSR-DB_full-16S_filt_taxa.qza GSR-DB_full-16S_filt_seqs.qza

        source $path/bin/activate seqkit

        seqkit faidx GSR-DB_full-16S_filt_seqs.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=1200 && $2<=1800) print $1}' GSR-DB_full-16S_filt_seqs.fasta.fai > gsr_filtered_ids

        seqkit faidx -X gsr_filtered_ids GSR-DB_full-16S_filt_seqs.fasta > GSR_final_seqs.fasta

        sed 's/ //g;s/;/\t/g;s/[k,p,c,o,f,g,s]__//g' GSR-DB_full-16S_filt_taxa.txt | \
        awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) for (i=2;i<=NF;i++) gsub(/_/, " ", $i)} 1' > temp && mv temp GSR-DB_full-16S_filt_taxa.txt

        grep -Ff gsr_filtered_ids GSR-DB_full-16S_filt_taxa.txt > temp && mv temp GSR-DB_full-16S_filt_taxa.txt

        source $path/bin/activate nanotaxi-env

        Rscript $script_dir/gsr_emu_db_build.r GSR-DB_full-16S_filt_taxa.txt

        source $path/bin/activate emu

        emu build-database GSR --sequences GSR_final_seqs.fasta --seq2tax seq2tax.map.tsv --taxonomy-list taxonomy.tsv

        rm -r GSR-DB_full-16S_filt_taxa.txt GSR_final_seqs.fasta* seq2tax.map.tsv taxonomy.tsv GSR-DB_full-16S_filt_seqs.fasta* gsr_filtered_ids

        cd GSR

        grep -qF "export EMU_GSR=\"$PWD\"" ~/.bashrc || echo "export EMU_GSR=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/EMU/GSR

grep -qF "export EMU_GSR=\"$PWD\"" ~/.bashrc || echo "export EMU_GSR=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/EMU/

if [ ! -d GTDB ]; then

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/bac120_ssu_reps.fna.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/ar53_ssu_reps.fna.gz

        zcat bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz > GTDB_16S_reps.fasta && rm -r bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz

        zcat bac120_metadata.tsv.gz ar53_metadata.tsv.gz | grep -v "ncbi" awk -F "\t" '{print $1"\t"$81}' > seq2tax.map.tsv

        zcat bac120_metadata.tsv.gz ar53_metadata.tsv.gz | grep -v "ncbi" | awk -F "\t" '{print $81"\t"$82}' | sed 's/;/\t/g;s/[d,p,c,o,f,g,s]__//g' | \
        awk 'BEGIN{FS=OFS="\t"}{print $1,$8,$7,$6,$5,$4,$3,$2}' - | sort -k1 -n -r | uniq | \
        cat <(echo -e "tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom") - > taxonomy.tsv && rm -r bac120_metadata.tsv.gz ar53_metadata.tsv.gz

        source $path/bin/activate seqkit

        seqkit faidx GTDB_16S_reps.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=1200 && $2<=1800) print $1}' GTDB_16S_reps.fasta.fai > gtdb_filtered_ids

        seqkit faidx -X gtdb_filtered_ids GTDB_16S_reps.fasta > temp && mv temp GTDB_16S_reps.fasta

        source $path/bin/activate emu
        
        emu build-database GTDB --sequences GTDB_16S_reps.fasta --seq2tax seq2tax.map.tsv --taxonomy-list taxonomy.tsv

        rm -r GTDB_16S_reps.fasta* gtdb_filtered_ids seq2tax.map.tsv taxonomy.tsv

        cd GTDB

        grep -qF "export EMU_GTDB=\"$PWD\"" ~/.bashrc || echo "export EMU_GTDB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/EMU/GTDB

grep -qF "export EMU_GTDB=\"$PWD\"" ~/.bashrc || echo "export EMU_GTDB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/EMU/

if [ ! -d MIMT ]; then

        wget -c https://people.biopolis.pt/bu/mimt/downloads/16S_files/MIMt-16S_M2c_25_10_taxid.fna.gz -O MIMt.fasta.gz && gunzip MIMt.fasta.gz

        wget -c https://people.biopolis.pt/bu/mimt/downloads/16S_files/MIMt-16S_M2c_25_10.tax.gz -O MIMT_taxa.txt.gz && gunzip MIMT_taxa.txt.gz

        grep ">" MIMt.fasta | sed 's/>//g' | awk -F "\t" '{if($2!=" ") print $0}' > seq2tax.map.tsv

        paste -d "\t" <(awk -F "\t" '{print $2}' seq2tax.map.tsv) <(awk -F "\t" '{print $1}' seq2tax.map.tsv | grep -Ff - MIMT_taxa.txt | awk -F "\t" '{print $2}') | \
        sort -u | sed 's/;/\t/g;s/[K,P,C,O,F,G,S]__//g' | awk 'BEGIN{FS="\t";OFS="\t"}{for (i=2;i<=NF;i++) gsub(/_/, " ", $i) split($8, a, " "); $8=a[1]" "a[2]} 1' | \
        cat <(echo -e "tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom") - | \
        awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1) print $0; else print $1,$8,$7,$6,$5,$4,$3,$2}' > taxonomy.tsv

        source $path/bin/activate seqkit

        seqkit faidx MIMt.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=1200 && $2<=1800) print $1}' MIMt.fasta.fai | grep -Ff - <(awk -F "\t" '{print $1}' seq2tax.map.tsv) > mimt_filtered_ids

        seqkit faidx -X mimt_filtered_ids MIMt.fasta > MIMT_final_seqs.fasta

        source $path/bin/activate emu

        emu build-database MIMT --sequences MIMT_final_seqs.fasta --seq2tax seq2tax.map.tsv --taxonomy-list taxonomy.tsv

        rm -r MIMT_final_seqs.fasta* seq2tax.map.tsv taxonomy.tsv mimt_filtered_ids MIMt.fasta* MIMT_taxa.txt

        cd MIMT

        grep -qF "export EMU_MIMT=\"$PWD\"" ~/.bashrc || echo "export EMU_MIMT=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/EMU/MIMT

grep -qF "export EMU_MIMT=\"$PWD\"" ~/.bashrc || echo "export EMU_MIMT=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/EMU/

if [ ! -d REFSEQ ]; then

        wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz

        zcat bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz > refseq_16S.fasta && rm -r bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz

        cp $TAXONKIT_DB/refseq_taxid.txt seq2tax.map.tsv

        threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)
        
        source $path/bin/activate taxonkit

        taxonkit reformat2 --data-dir $TAXONKIT_DB -f "{domain};{phylum};{class};{order};{family};{genus};{species}" -I 2 seq2tax.map.tsv | \
        awk 'BEGIN{FS=OFS="\t"}{print $2,$3}' | sort | uniq | sed 's/;/\t/g' | \
        cat <(echo -e "tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom") - | \
        awk 'BEGIN{FS=OFS="\t"}{if(NR==1) print $0; else print $1,$8,$7,$6,$5,$4,$3,$2}' > taxonomy.tsv

        source $path/bin/activate seqkit

        seqkit faidx refseq_16S.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=1200 && $2<=1800) print $1}' refseq_16S.fasta.fai > refseq_filtered_ids

        seqkit faidx -X refseq_filtered_ids refseq_16S.fasta > refseq_final_seqs.fasta

        source $path/bin/activate emu

        emu build-database REFSEQ --sequences refseq_final_seqs.fasta --seq2tax seq2tax.map.tsv --taxonomy-list taxonomy.tsv

        rm -r refseq_final_seqs.fasta* refseq_filtered_ids refseq_16S.fasta* seq2tax.map.tsv taxonomy.tsv

        cd REFSEQ

        grep -qF "export EMU_REFSEQ=\"$PWD\"" ~/.bashrc || echo "export EMU_REFSEQ=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/EMU/REFSEQ

grep -qF "export EMU_REFSEQ=\"$PWD\"" ~/.bashrc || echo "export EMU_REFSEQ=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir
