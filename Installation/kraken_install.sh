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

path=$(which conda | sed "s/\b\/conda\b//g")

base_dir=$PWD

source ~/.bashrc

if { conda env list |  grep "kraken2"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name kraken2 --file kraken.txt

fi

if [ ! -d DATA ]; then
        
	mkdir DATA
fi

cd DATA

if [ ! -d KRAKEN_DATA ]; then
                
        source $path/activate kraken2

        kraken2-build --db KRAKEN_DATA --special greengenes

        cd KRAKEN_DATA

        grep -qF "export KRAKEN_DB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_DB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        source $path/activate base


fi

cd $base_dir/DATA/KRAKEN_DATA

grep -qF "export KRAKEN_DB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir

if { conda env list | grep "seqkit";} > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name seqkit --file seqkit.txt

fi

cd $base_dir

if { conda env list | grep "taxonkit";} > /dev/null 2>&1; then

        echo "Environment Exist"

else
        
        conda create -n taxonkit --file taxonkit.txt
        
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