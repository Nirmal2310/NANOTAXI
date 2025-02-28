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

if { conda env list |  grep "kraken2"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name kraken2 --file kraken.txt

fi

if { conda env list | grep "taxonkit";} > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create -n taxonkit --file taxonkit.txt -y

fi

if { conda env list | grep "nanofilt";} > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name nanofilt --file nanofilt.txt -y

fi

if { conda env list | grep "bbtools";} > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name bbtools --file bbtools.txt -y

fi

if { conda env list | grep "seqkit";} > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name seqkit --file seqkit.txt -y

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

fi

cd $base_dir/DATA/TAXONKIT_DATA

grep -qF "export TAXONKIT_DB=\"$PWD\"" ~/.bashrc || echo "export TAXONKIT_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir

if [ ! -d DATA ]; then

        mkdir DATA
fi

cd DATA

if [ ! -d KRAKEN_NCBI ]; then

        source $path/bin/activate kraken2

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/bac120_ssu_reps.fna.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/ar53_ssu_reps.fna.gz

        zcat bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz > GTDB_16S_reps.fasta && rm -r bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz

        zcat bac120_metadata.tsv.gz ar53_metadata.tsv.gz | awk -F "\t" '{if(NR>1) print $1"\t"$81}' > seqid_taxid.txt && rm -r bac120_metadata.tsv.gz ar53_metadata.tsv.gz

        sed -i 's/ .*$//g' GTDB_16S_reps.fasta

        grep ">" GTDB_16S_reps.fasta | sed 's/>//g' | split -l 1000 - ids_chunk_

        source $path/bin/activate seqkit

        threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)

        chunk_number=$(ls ids_chunk_* | wc -l)

        parallel_jobs=$(if [ $chunk_number -gt $threads ]; then echo $threads; else echo $chunk_number; fi)
        
        parallel -j $parallel_jobs "rg -f {} seqid_taxid.txt" ::: ids_chunk_* | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$1"|kraken:taxid|"$2}' > seq_id_replacement.txt && rm -r ids_chunk_*

        seqkit replace -p '^(\S+)' -r '{kv}$2' -k seq_id_replacement.txt GTDB_16S_reps.fasta > GTDB_16S_kraken2_ready.fasta

        source $path/bin/activate kraken2

        kraken2-build --download-taxonomy --db KRAKEN_GTDB --use-ftp --skip-maps

        kraken2-build --add-to-library GTDB_16S_kraken2_ready.fasta --db KRAKEN_GTDB

        grep ">" GTDB_16S_kraken2_ready.fasta | sed 's/>//g' | awk '{split($1,a,"|"); print $1"\t"a[3]}' > KRAKEN_GTDB/seqid2taxid.map

        kraken2-build --build --db KRAKEN_GTDB --threads $threads

        kraken2-build --clean --db KRAKEN_GTDB

        rm -r GTDB_16S_kraken2_ready.fasta GTDB_16S_reps.fasta seq_id_replacement.txt seqid_taxid.txt
        
        cd KRAKEN_GTDB

        grep -qF "export KRAKEN_GTDB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_GTDB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        source $path/bin/activate base

        cd $base_dir

fi

cd $base_dir/DATA/KRAKEN_GTDB

grep -qF "export KRAKEN_GTDB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_GTDB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir

if { conda env list |  grep "minknow_api"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name minknow_api --file minknow_api.txt -y

        conda activate minknow_api

        pip install minknow-api==6.0.4

        pip install grpcio==1.60.1

        pip install pandas==2.2.2

        conda activate base
fi
