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

        echo "Environment Exist"

else

        conda create --name kraken2 --file $script_dir/kraken.txt

fi

if [ ! -d DATA ]; then
        
	mkdir DATA
fi

if { conda env list | grep "nanofilt";} > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name nanofilt --file $script_dir/nanofilt.txt

fi

if { conda env list | grep "seqkit";} > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name seqkit --file $script_dir/seqkit.txt

fi

if { conda env list | grep "taxonkit";} > /dev/null 2>&1; then

        echo "Environment Exist"

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

if [ ! -d KRAKEN_DATA ]; then
                
        source $path/bin/activate kraken2

        kraken2-build --download-taxonomy --db KRAKEN_DATA --use-ftp  --skip-maps

 	wget -c https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz -O SILVA_138.2_ref.fasta.gz && gunzip SILVA_138.2_ref.fasta.gz

  	source $path/bin/activate taxonkit

   	grep ">" SILVA_138.2_ref.fasta | grep "Bacteria\|Archaea" | sed 's/>//g;s/ /\t/' | sed 's/;/\t/g' | awk -F "\t" '{print $1"\t"$8}' | \
    		taxonkit name2taxid -i 2 --data-dir $TAXONKIT_DB | awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $1"|kraken:taxid|"$3}' > seq_id_replacement.txt

      	source $path/bin/activate seqkit

       	awk -F "\t" '{print $1}' seq_id_replacement.txt | seqkit faidx -X - SILVA_138.2_ref.fasta | \
		awk 'BEGIN{RS=">";FS="\n"}NR>1{printf ">%s\n",$1;for (i=2;i<=NF;i++) {gsub(/U/,"T",$i); printf "%s\n",$i}}' | sed '/^$/d' > SILVA_138.2_16S.fasta

  	seqkit replace -p '^(\S+)' -r '{kv}$2' -k seq_id_replacement.txt SILVA_138.2_16S.fasta > SILVA_138.2_16S_kraken2_ready.fasta

   	threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)

   	source $path/bin/activate kraken2

    	kraken2-build --add-to-library SILVA_138.2_16S_kraken2_ready.fasta --db KRAKEN_DATA

     	grep ">" SILVA_138.2_16S_kraken2_ready.fasta | sed 's/>//g' | awk '{split($1,a,"|"); print $1"\t"a[3]}' > KRAKEN_DATA/seqid2taxid.map

     	kraken2-build --build --db KRAKEN_DATA --threads $threads

      	kraken2-build --clean --db KRAKEN_DATA

       	rm -r SILVA_138.2_16S_kraken2_ready.fasta SILVA_138.2_ref.fasta* SILVA_138.2_16S.fasta* seq_id_replacement.txt

        cd KRAKEN_DATA

        grep -qF "export KRAKEN_DB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_DB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        source $path/bin/activate base

fi

cd $base_dir/DATA/KRAKEN_DATA

grep -qF "export KRAKEN_DB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir
