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

if [ ! -d DATA ]; then

        mkdir DATA
fi

cd DATA

if [ ! -d KRAKEN_NCBI ]; then

        source $path/bin/activate kraken2

        wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz -O archaea.16SrRNA.fasta.gz

        wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz -O bacteria.16SrRNA.fasta.gz

        zcat bacteria.16SrRNA.fasta.gz archaea.16SrRNA.fasta.gz > sequences.fasta

        rm -r archaea.16SrRNA.fasta.gz bacteria.16SrRNA.fasta.gz

        kraken2-build --download-taxonomy --db KRAKEN_NCBI --use-ftp

        kraken2-build --add-to-library sequences.fasta --db KRAKEN_NCBI --use-ftp

        kraken2-build --build --db KRAKEN_NCBI --use-ftp && rm -r sequences.fasta

        cd KRAKEN_NCBI

        grep -qF "export KRAKEN_NCBI=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_NCBI=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        source $path/bin/activate base

        cd $base_dir

fi

cd $base_dir/DATA/KRAKEN_NCBI

grep -qF "export KRAKEN_NCBI=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_NCBI=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir

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

if { conda env list |  grep "minknow_api"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name minknow_api --file minknow_api.txt -y

        conda activate minknow_api

        pip install minknow-api==6.0.4

        pip install grpcio==1.60.1

        pip install pandas==2.2.2

        conda activate base
