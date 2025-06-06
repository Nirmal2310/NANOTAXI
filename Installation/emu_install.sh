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
  		
        fi

        conda activate base

else

        conda create --name emu --file $script_dir/emu.txt

        source $path/bin/activate emu

        pip install osfclient

        source $path/bin/activate base

fi

if { conda env list | grep "nanofilt";} > /dev/null 2>&1; then

	conda list -n nanofilt --explicit > _current_env.txt

        if diff -q $script_dir/nanofilt.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name nanofilt --file $script_dir/nanofilt.txt -y && rm -r _current_env.txt
  		
        fi

        conda activate base        

else
        
        conda create -n nanofilt --file $script_dir/nanofilt.txt
        
fi


if [ ! -d DATA ]; then
        
	mkdir DATA
fi

cd DATA

if [ ! -d EMU_DATA ]; then
        
        source $path/bin/activate emu
        
        mkdir EMU_DATA

        osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar.gz

        tar -xvf emu.tar.gz -C EMU_DATA --strip 1

        rm -r emu.tar.gz

        cd EMU_DATA
        
        source $path/bin/activate base

        grep -qF "export EMU_DB=\"$PWD\"" ~/.bashrc || echo "export EMU_DB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/EMU_DATA

grep -qF "export EMU_DB=\"$PWD\"" ~/.bashrc || echo "export EMU_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir

if { conda env list | grep "taxonkit";} > /dev/null 2>&1; then

        conda list -n taxonkit --explicit > _current_env.txt

        if diff -q $script_dir/taxonkit.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name taxonkit --file $script_dir/taxonkit.txt -y && rm -r _current_env.txt
        fi

        conda activate base

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

fi

cd $base_dir/DATA/TAXONKIT_DATA

grep -qF "export TAXONKIT_DB=\"$PWD\"" ~/.bashrc || echo "export TAXONKIT_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir
