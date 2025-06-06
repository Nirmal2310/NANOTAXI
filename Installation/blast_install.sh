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
        fi

        conda activate base

else

        conda create --name blast --file $script_dir/blast.txt

fi


if [ ! -d DATA ]; then
        
	mkdir DATA
fi

cd DATA

# Download 16S rRNA Database (Required by BLAST)

if [ ! -d 16s_DATA ]; then
                
	mkdir 16s_DATA

        cd 16s_DATA
        
        wget https://ftp.ncbi.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
        
        tar -xvf 16S_ribosomal_RNA.tar.gz
        
        rm 16S_ribosomal_RNA.tar.gz

        grep -qF "export BLAST_DB=\"$PWD\"" ~/.bashrc || echo "export BLAST_DB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc
fi
        
cd $base_dir/DATA/16s_DATA

grep -qF "export BLAST_DB=\"$PWD\"" ~/.bashrc || echo "export BLAST_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

source $path/bin/activate base
        
cd $base_dir

if { conda env list | grep "nanofilt";} > /dev/null 2>&1; then

        conda list -n nanofilt --explicit > _current_env.txt

        if diff -q $script_dir/nanofilt.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name nanofilt --file $script_dir/nanofilt.txt -y && rm -r _current_env.txt
        fi

else
        
        conda create -n nanofilt --file $script_dir/nanofilt.txt
        
fi

if { conda env list | grep "taxonkit";} > /dev/null 2>&1; then

        conda list -n taxonkit --explicit > _current_env.txt

        if diff -q $script_dir/taxonkit.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name taxonkit --file $script_dir/taxonkit.txt -y && rm -r _current_env.txt
        fi

else
        
        conda create -n taxonkit --file $script_dir/taxonkit.txt
        
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
