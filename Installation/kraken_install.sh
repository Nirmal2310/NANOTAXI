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

script_dir=$(dirname "$(readlink -f "$0")")

source ~/.bashrc

if { conda env list |  grep "kraken2"; } > /dev/null 2>&1; then

        conda list -n kraken2 --explicit > _current_env.txt

        if diff -q $script_dir/kraken.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name kraken2 --file $script_dir/kraken.txt -y && rm -r _current_env.txt
		conda list -n kraken2 --explicit > $script_dir/kraken.txt
        fi

else

        conda create --name kraken2 --file $script_dir/kraken.txt
	conda list -n kraken2 --explicit > $script_dir/kraken.txt

fi

if [ ! -d DATA ]; then
        
	mkdir DATA
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

        conda create --name nanofilt --file $script_dir/nanofilt.txt
	conda list -n nanofilt --explicit > $script_dir/nanofilt.txt

fi

if { conda env list | grep "seqkit";} > /dev/null 2>&1; then

        conda list -n seqkit --explicit > _current_env.txt

        if diff -q $script_dir/seqkit.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name seqkit --file $script_dir/seqkit.txt -y && rm -r _current_env.txt
		conda list -n seqkit --explicit > $script_dir/seqkit.txt
        fi

else

        conda create --name seqkit --file $script_dir/seqkit.txt
	conda list -n seqkit --explicit > $script_dir/seqkit.txt

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

fi

cd $base_dir/DATA/TAXONKIT_DATA

grep -qF "export TAXONKIT_DB=\"$PWD\"" ~/.bashrc || echo "export TAXONKIT_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir

cd DATA

if [ ! -d KRAKEN_DATA ]; then
                
        source $path/bin/activate kraken2

        kraken2-build --download-taxonomy --db KRAKEN_DATA --use-ftp  --skip-maps

 	wget -c https://people.biopolis.pt/bu/mimt/downloads/16S_files/MIMt-16S_M2c_25_10_taxid.fna.gz -O MIMt.fasta.gz && gunzip MIMt.fasta.gz

        grep ">" MIMt.fasta | sed 's/>//g' | awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $1"|kraken:taxid|"$2}' > seq_id_replacement.txt

      	source $path/bin/activate seqkit

        seqkit replace -p '^(\S+)(.*)' -r '{kv}' -k seq_id_replacement.txt MIMt.fasta > MIMt_kraken2_ready.fasta

   	threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)

   	source $path/bin/activate kraken2

    	kraken2-build --add-to-library MIMt_kraken2_ready.fasta --db KRAKEN_DATA

     	grep ">" MIMt_kraken2_ready.fasta | sed 's/>//g' | awk '{split($1,a,"|"); print $1"\t"a[3]}' > KRAKEN_DATA/seqid2taxid.map

     	kraken2-build --build --db KRAKEN_DATA --threads $threads

      	kraken2-build --clean --db KRAKEN_DATA

       	rm -r MIMt.fasta seq_id_replacement.txt MIMt_kraken2_ready.fasta

        cd KRAKEN_DATA

        grep -qF "export KRAKEN_DB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_DB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        source $path/bin/activate base

fi

cd $base_dir/DATA/KRAKEN_DATA

grep -qF "export KRAKEN_DB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir
