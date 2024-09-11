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

if { conda env list |  grep "blast"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name blast --file blast.txt

fi


if [ ! -d DATA ]; then
        
	mkdir DATA
fi

cd DATA

if [ ! -d 16s_DATA ]; then
                
	mkdir 16s_DATA
fi
        
# Download 16S rRNA Database (Required by BLAST)
        
cd 16s_DATA
        
wget https://ftp.ncbi.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
        
tar -xvf 16S_ribosomal_RNA.tar.gz
        
rm 16S_ribosomal_RNA.tar.gz

grep -qF "export BLAST_DB=\"$PWD\"" ~/.bashrc || echo "export BLAST_DB=\"$PWD\"" >> ~/.bashrc
                
source $path/bin/activate base
        
cd $base_dir

if { conda env list | grep "nanofilt";} > /dev/null 2>&1; then

        echo "Environment Exist"

else
        
        conda create -n nanofilt --file nanofilt.txt
        
fi

if { conda env list | grep "taxonkit";} > /dev/null 2>&1; then

        echo "Environment Exist"

else
        
        conda create -n taxonkit --file taxonkit.txt
        
fi

cd DATA

if [ ! -d TAXONKIT_DATA ]; then
                
        mkdir TAXONKIT_DATA
fi

cd TAXONKIT_DATA

wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 

tar -zxvf taxdump.tar.gz

rm -r taxdump.tar.gz

grep -qF "export TAXONKIT_DB=\"$PWD\"" ~/.bashrc || echo "export TAXONKIT_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc