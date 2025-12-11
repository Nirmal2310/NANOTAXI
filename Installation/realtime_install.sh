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

if { conda env list |  grep -w "kraken2"; } > /dev/null 2>&1; then

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

if { conda env list | grep -w "taxonkit";} > /dev/null 2>&1; then

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

if { conda env list | grep -w "nanofilt";} > /dev/null 2>&1; then

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

if { conda env list | grep -w "bbtools";} > /dev/null 2>&1; then

        conda list -n bbtools --explicit > _current_env.txt

        if diff -q $script_dir/bbtools.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name bbtools --file $script_dir/bbtools.txt -y && rm -r _current_env.txt 
                conda list -n bbtools --explicit > $script_dir/bbtools.txt
        fi

else

        conda create --name bbtools --file $script_dir/bbtools.txt
        conda list -n bbtools --explicit > $script_dir/bbtools.txt

fi

if { conda env list | grep -w "seqkit";} > /dev/null 2>&1; then

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

if { conda env list |  grep -w "minknow_api"; } > /dev/null 2>&1; then

        conda list -n minknow_api --explicit > _current_env.txt

        if diff -q $script_dir/minknow_api.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name minknow_api --file $script_dir/minknow_api.txt -y && rm -r _current_env.txt

                conda init bash

                conda activate minknow_api

                pip install minknow-api==6.0.4
        
                pip install pysam

                pip install osfclient

                conda activate base
                
                conda list -n minknow_api --explicit > $script_dir/minknow_api.txt
                
        fi

else

        conda create --name minknow_api --file $script_dir/minknow_api.txt -y

        conda init bash

        conda activate minknow_api

        pip install minknow-api==6.0.4

        pip install pysam

        pip install osfclient

        conda activate base

        conda list -n minknow_api --explicit > $script_dir/minknow_api.txt
fi

if { conda env list |  grep -w "minimap2"; } > /dev/null 2>&1; then

        conda list -n minimap2 --explicit > _current_env.txt

        if diff -q $script_dir/minimap2.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name minimap2 --file $script_dir/minimap2.txt -y && rm -r _current_env.txt
                conda list -n minimap2 --explicit > $script_dir/minimap2.txt
        fi

else
        conda create --name minimap2 --file $script_dir/minimap2.txt
        conda list -n minimap2 --explicit > $script_dir/minimap2.txt

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

        source ~/.bashrc

        wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

        wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz

        zcat bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz > refseq_16S.fasta && rm -r bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz

        zcat nucl_gb.accession2taxid.gz | grep -w -f <(grep -oP '^>[^\s]+' refseq_16S.fasta | sed 's/^>//') | \
        awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3}' > refseq_taxid.txt

        rm -r nucl_gb.accession2taxid.gz refseq_16S.fasta

fi

cd $base_dir/DATA/TAXONKIT_DATA

grep -qF "export TAXONKIT_DB=\"$PWD\"" ~/.bashrc || echo "export TAXONKIT_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir

cd DATA

if [ ! -d KRAKEN ]; then
        
        mkdir KRAKEN
fi

cd KRAKEN

if [ ! -d GTDB ]; then

        source $path/bin/activate kraken2

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/bac120_ssu_reps.fna.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/ar53_ssu_reps.fna.gz

        zcat bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz > GTDB_16S_reps.fasta && rm -r bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz

        zcat bac120_metadata.tsv.gz ar53_metadata.tsv.gz | awk -F "\t" '{if(NR>1) print $1"\t"$81}' > seqid_taxid.txt && rm -r bac120_metadata.tsv.gz ar53_metadata.tsv.gz

        sed -i 's/ .*$//g' GTDB_16S_reps.fasta

        source $path/bin/activate seqkit

        seqkit faidx GTDB_16S_reps.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' GTDB_16S_reps.fasta.fai > gtdb_filtered_ids

        seqkit faidx -X gtdb_filtered_ids GTDB_16S_reps.fasta > temp && mv temp GTDB_16S_reps.fasta

        grep ">" GTDB_16S_reps.fasta | sed 's/>//g' | split -l 1000 - ids_chunk_

        threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)

        chunk_number=$(ls ids_chunk_* | wc -l)

        parallel_jobs=$(if [ $chunk_number -gt $threads ]; then echo $threads; else echo $chunk_number; fi)
        
        parallel -j $parallel_jobs "rg -f {} seqid_taxid.txt" ::: ids_chunk_* | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$1"|kraken:taxid|"$2}' > seq_id_replacement.txt && rm -r ids_chunk_*

        seqkit replace -p '^(\S+)' -r '{kv}$2' -k seq_id_replacement.txt GTDB_16S_reps.fasta > GTDB_16S_kraken2_ready.fasta

        source $path/bin/activate kraken2

        kraken2-build --download-taxonomy --db GTDB --use-ftp --skip-maps

        kraken2-build --add-to-library GTDB_16S_kraken2_ready.fasta --db GTDB --no-masking

        grep ">" GTDB_16S_kraken2_ready.fasta | sed 's/>//g' | awk '{split($1,a,"|"); print $1"\t"a[3]}' > GTDB/seqid2taxid.map

        kraken2-build --build --db GTDB --threads $threads

        kraken2-build --clean --db GTDB

        rm -r GTDB_16S_kraken2_ready.fasta* GTDB_16S_reps.fasta* seq_id_replacement.txt seqid_taxid.txt gtdb_filtered_ids
        
        cd GTDB

        grep -qF "export KRAKEN_GTDB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_GTDB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        source $path/bin/activate base

        cd $base_dir

fi

cd $base_dir/DATA/KRAKEN/GTDB

grep -qF "export KRAKEN_GTDB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_GTDB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/KRAKEN

if [ ! -d MIMT ]; then

        source $path/bin/activate kraken2

        kraken2-build --download-taxonomy --db MIMT --use-ftp  --skip-maps

 	wget -c https://people.biopolis.pt/bu/mimt/downloads/16S_files/MIMt-16S_M2c_25_10_taxid.fna.gz -O MIMt.fasta.gz && gunzip MIMt.fasta.gz

        grep ">" MIMt.fasta | sed 's/>//g' | awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $1"|kraken:taxid|"$2}' > seq_id_replacement.txt

      	source $path/bin/activate seqkit

        seqkit replace -p '^(\S+)(.*)' -r '{kv}' -k seq_id_replacement.txt MIMt.fasta > MIMt_kraken2_ready.fasta

   	threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)

   	source $path/bin/activate kraken2

    	kraken2-build --add-to-library MIMt_kraken2_ready.fasta --db MIMT --no-masking

     	grep ">" MIMt_kraken2_ready.fasta | sed 's/>//g' | awk '{split($1,a,"|"); print $1"\t"a[3]}' > MIMT/seqid2taxid.map

     	kraken2-build --build --db MIMT --threads $threads

      	kraken2-build --clean --db MIMT

       	rm -r MIMt.fasta* seq_id_replacement.txt MIMt_kraken2_ready.fasta*

        cd MIMT

        grep -qF "export KRAKEN_MIMT=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_MIMT=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        source $path/bin/activate base

fi

cd $base_dir/DATA/KRAKEN/MIMT

grep -qF "export KRAKEN_MIMT=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_MIMT=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/KRAKEN

if [ ! -d GSR ]; then

        wget -c https://manichanh.vhir.org/gsrdb/GSR-DB_full-16S.tar.gz

        tar -xvf GSR-DB_full-16S.tar.gz && rm -r GSR-DB_full-16S.tar.gz GSR-DB_full-16S_filt_taxa.qza GSR-DB_full-16S_filt_seqs.qza

        threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)
        
        source $path/bin/activate seqkit

        seqkit faidx GSR-DB_full-16S_filt_seqs.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' GSR-DB_full-16S_filt_seqs.fasta.fai > gsr_filtered_ids

        seqkit faidx -X gsr_filtered_ids GSR-DB_full-16S_filt_seqs.fasta > temp && mv temp GSR-DB_full-16S_filt_seqs.fasta

        grep -Ff gsr_filtered_ids GSR-DB_full-16S_filt_taxa.txt > temp && mv temp GSR-DB_full-16S_filt_taxa.txt

        source $path/bin/activate nanotaxi-env

        Rscript $script_dir/gsr_kraken_taxa_build.R GSR-DB_full-16S_filt_taxa.txt 16

        mkdir -p GSR GSR/taxonomy GSR/library

        mv nodes.dmp GSR/taxonomy/ && mv names.dmp GSR/taxonomy/

        sed -i -e 's/$/\t|/g' GSR/taxonomy/nodes.dmp

        sed -i -e 's/$/\t|/g' GSR/taxonomy/names.dmp

        source $path/bin/activate seqkit

        seqkit replace -p '^(\S+)' -r '{kv}$2' -k seq_id_replacement.txt GSR-DB_full-16S_filt_seqs.fasta > GSR_kraken2_ready.fasta
        
        source $path/bin/activate kraken2

        kraken2-build --add-to-library GSR_kraken2_ready.fasta --db GSR --no-masking

        kraken2-build --build --db GSR --threads $threads

        kraken2-build --clean --db GSR

        rm -r GSR_kraken2_ready.fasta* GSR-DB_full-16S_filt_seqs.fasta* seq_id_replacement.txt gsr_filtered_ids GSR-DB_full-16S_filt_taxa.txt
        
        cd GSR

        grep -qF "export KRAKEN_GSR=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_GSR=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        source $path/bin/activate base

        cd $base_dir

fi

cd $base_dir/DATA/KRAKEN/GSR

grep -qF "export KRAKEN_GSR=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_GSR=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/KRAKEN

if [ ! -d REFSEQ ]; then

        wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz

        zcat bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz > refseq_16S.fasta && rm -r bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz

        threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}'; fi)

        cp $TAXONKIT_DB/refseq_taxid.txt seqid_taxid.txt

        source $path/bin/activate seqkit

        seqkit faidx refseq_16S.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' refseq_16S.fasta.fai > refseq_filtered_ids

        seqkit faidx -X refseq_filtered_ids refseq_16S.fasta > temp && mv temp refseq_16S.fasta

        grep ">" refseq_16S.fasta | sed 's/>//g' | split -l 1000 - ids_chunk_

        chunk_number=$(ls ids_chunk_* | wc -l)

        parallel_jobs=$(if [ $chunk_number -gt $threads ]; then echo $threads; else echo $chunk_number; fi)
        
        parallel -j $parallel_jobs "rg -f {} seqid_taxid.txt" ::: ids_chunk_* | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$1"|kraken:taxid|"$2}' > seq_id_replacement.txt && rm -r ids_chunk_*

        seqkit replace -p '^(\S+)' -r '{kv}$2' -k seq_id_replacement.txt refseq_16S.fasta > refseq_kraken2_ready.fasta

        source $path/bin/activate kraken2

        kraken2-build --download-taxonomy --db REFSEQ --use-ftp  --skip-maps

        kraken2-build --add-to-library refseq_kraken2_ready.fasta --db REFSEQ --no-masking

        grep ">" refseq_kraken2_ready.fasta | sed 's/>//g' | awk '{split($1,a,"|"); print $1"\t"a[3]}' > REFSEQ/seqid2taxid.map

     	kraken2-build --build --db REFSEQ --threads $threads

      	kraken2-build --clean --db REFSEQ

       	rm -r refseq_16S.fasta* seq_id_replacement.txt refseq_kraken2_ready.fasta* seqid_taxid.txt refseq_filtered_ids

        cd REFSEQ

        grep -qF "export KRAKEN_REFSEQ=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_REFSEQ=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        source $path/bin/activate base

fi

cd $base_dir/DATA/KRAKEN/REFSEQ

grep -qF "export KRAKEN_REFSEQ=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_REFSEQ=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/KRAKEN/

if [ ! -d EMUDB ]; then

        source $path/bin/activate minknow_api

        osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar.gz

        tar -xvf emu.tar.gz --strip 1

        rm -r emu.tar.gz

        threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}'; fi)

        sed -i 's/ .*$//g;s/:/_/g' species_taxid.fasta

        grep ">" species_taxid.fasta | sed 's/>//g' | awk 'BEGIN{FS=OFS="\t"}{$2=$1; gsub(/_.*$/, "",$2); print $1, $1"|kraken:taxid|"$2}' > seq_id_replacement.txt

        source $path/bin/activate seqkit

        seqkit replace -p '^(\S+)' -r '{kv}$2' -k seq_id_replacement.txt species_taxid.fasta > emu_kraken2_ready.fasta

        source $path/bin/activate kraken2

        kraken2-build --download-taxonomy --db EMUDB --use-ftp  --skip-maps

        kraken2-build --add-to-library emu_kraken2_ready.fasta --db EMUDB --no-masking

        grep ">" emu_kraken2_ready.fasta | sed 's/>//g' | awk '{split($1,a,"|"); print $1"\t"a[3]}' > EMUDB/seqid2taxid.map

     	kraken2-build --build --db EMUDB --threads $threads

      	kraken2-build --clean --db EMUDB

       	rm -r species_taxid.fasta emu_kraken2_ready.fasta seq_id_replacement.txt taxonomy.tsv

        cd EMUDB

        grep -qF "export KRAKEN_EMUDB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_EMUDB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

fi

cd $base_dir/DATA/KRAKEN/EMUDB

grep -qF "export KRAKEN_EMUDB=\"$PWD\"" ~/.bashrc || echo "export KRAKEN_EMUDB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA

if [ ! -d MINIMAP2 ]; then
        mkdir MINIMAP2
fi

cd $base_dir/DATA/MINIMAP2

if [ ! -d GSR_DB ]; then

        mkdir GSR_DB

        wget -c https://manichanh.vhir.org/gsrdb/GSR-DB_full-16S.tar.gz

        tar -xvf GSR-DB_full-16S.tar.gz && rm -r GSR-DB_full-16S.tar.gz GSR-DB_full-16S_filt_taxa.qza GSR-DB_full-16S_filt_seqs.qza

        sed 's/ //g;s/;/\t/g;s/[k,p,c,o,f,g,s]__//g' GSR-DB_full-16S_filt_taxa.txt | \
        awk 'BEGIN{FS="\t";OFS="\t"}{for (i=2;i<=NF;i++) gsub(/_/, " ", $i)} 1' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > temp && mv temp GSR-DB_full-16S_filt_taxa.txt

        source $path/bin/activate seqkit

        seqkit faidx GSR-DB_full-16S_filt_seqs.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' GSR-DB_full-16S_filt_seqs.fasta.fai > gsr_filtered_ids

        seqkit faidx -X gsr_filtered_ids GSR-DB_full-16S_filt_seqs.fasta > temp && mv temp GSR-DB_full-16S_filt_seqs.fasta

        mv GSR-DB_full-16S_filt_taxa.txt GSR-DB_full-16S_filt_seqs.fasta GSR_DB/

        rm -r GSR-DB_full-16S_filt_seqs.fasta.fai gsr_filtered_ids

        cd GSR_DB

        grep -qF "export GSR_DB=\"$PWD\"" ~/.bashrc || echo "export GSR_DB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        cd $base_dir

fi

cd $base_dir/DATA/MINIMAP2/GSR_DB

grep -qF "export GSR_DB=\"$PWD\"" ~/.bashrc || echo "export GSR_DB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/MINIMAP2

if [ ! -d GTDB ]; then

        mkdir GTDB

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/bac120_ssu_reps.fna.gz

        wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/ar53_ssu_reps.fna.gz

        zcat bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz > GTDB_16S_reps.fasta && rm -r bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz

        grep ">" GTDB_16S_reps.fasta | sed 's/>//g' | sed 's/ /\t/;s/;/\t/g;s/[d,p,c,o,f,g,s]__//g' | sed 's/ \[locus.*$//g' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > GTDB_taxa.txt

        source $path/bin/activate seqkit

        seqkit faidx GTDB_16S_reps.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' GTDB_16S_reps.fasta.fai > gtdb_filtered_ids

        seqkit faidx -X gtdb_filtered_ids GTDB_16S_reps.fasta > GTBD_final_seqs.fasta

        rm -r gtdb_filtered_ids GTDB_16S_reps.fasta*

        mv GTBD_final_seqs.fasta GTDB_taxa.txt GTDB/

        cd GTDB

        grep -qF "export GTDB=\"$PWD\"" ~/.bashrc || echo "export GTDB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        cd $base_dir

fi

cd $base_dir/DATA/MINIMAP2/GTDB

grep -qF "export GTDB=\"$PWD\"" ~/.bashrc || echo "export GTDB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/MINIMAP2

if [ ! -d MIMT ]; then

        mkdir MIMT

        wget -c https://people.biopolis.pt/bu/mimt/downloads/16S_files/MIMt-16S_M2c_25_10_taxid.fna.gz -O MIMt.fasta.gz && gunzip MIMt.fasta.gz

        wget -c https://people.biopolis.pt/bu/mimt/downloads/16S_files/MIMt-16S_M2c_25_10.tax.gz -O MIMT_taxa.txt.gz && gunzip MIMT_taxa.txt.gz

        sed 's/;/\t/g;s/[K,P,C,O,F,G,S]__//g' MIMT_taxa.txt | awk 'BEGIN{FS="\t";OFS="\t"}{if(NR>1) for (i=2;i<=NF;i++) gsub(/_/, " ", $i)} 1' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > temp && mv temp MIMT_taxa.txt

        source $path/bin/activate seqkit

        seqkit faidx MIMt.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' MIMt.fasta.fai > mimt_filtered_ids

        seqkit faidx -X mimt_filtered_ids MIMt.fasta > MIMT_final_seqs.fasta

        rm -r mimt_filtered_ids MIMt.fasta*

        mv MIMT_final_seqs.fasta MIMT_taxa.txt MIMT/

        cd MIMT

        grep -qF "export MIMT=\"$PWD\"" ~/.bashrc || echo "export MIMT=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        cd $base_dir

fi

cd $base_dir/DATA/MINIMAP2/MIMT

grep -qF "export MIMT=\"$PWD\"" ~/.bashrc || echo "export MIMT=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/MINIMAP2

if [ ! -d REFSEQ ]; then

        mkdir REFSEQ

        wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz

        zcat bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz > refseq_16S.fasta && rm -r bacteria.16SrRNA.fna.gz archaea.16SrRNA.fna.gz

        threads=$(if [ $(nproc) -gt 16 ]; then echo 16; else echo $(nproc) | awk '{print $1/2}' ; fi)

        source $path/bin/activate taxonkit

        cp $TAXONKIT_DB/refseq_taxid.txt seqid_taxid.txt

        taxonkit lineage --data-dir $TAXONKIT_DB --threads $threads -i 2 seqid_taxid.txt | sed 's/;/\t/g' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - | \
        awk 'BEGIN{FS=OFS="\t"}{if(NR==1) print $0; else print $1,$4,$5,$6,$7,$8,$9,$10,$11}' > RefSeq_taxa.txt

        source $path/bin/activate seqkit

        seqkit faidx refseq_16S.fasta

        awk 'BEGIN{FS="\t";OFS="\t"}{if($2>=900 && $2<=1800) print $1}' refseq_16S.fasta.fai > refseq_filtered_ids

        seqkit faidx -X refseq_filtered_ids refseq_16S.fasta > refseq_final_seqs.fasta

        rm -r refseq_filtered_ids refseq_16S.fasta* seqid_taxid.txt

        mv refseq_final_seqs.fasta RefSeq_taxa.txt REFSEQ/

        cd REFSEQ

        grep -qF "export REFSEQ=\"$PWD\"" ~/.bashrc || echo "export REFSEQ=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        cd $base_dir

fi

cd $base_dir/DATA/MINIMAP2/REFSEQ

grep -qF "export REFSEQ=\"$PWD\"" ~/.bashrc || echo "export REFSEQ=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir/DATA/MINIMAP2

if [ ! -d EMUDB ]; then

        source $path/bin/activate minknow_api

        mkdir EMUDB

        osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar.gz

        tar -xvf emu.tar.gz -C EMUDB --strip 1

        rm -r emu.tar.gz

        cd EMUDB

        sed -i 's/ .*$//g' species_taxid.fasta && mv species_taxid.fasta EMU.fasta

        join <(grep ">" EMU.fasta | sed 's/>//;s/ .*$//g' | awk 'BEGIN{FS=OFS="\t"}{$2=$1; gsub(/:.*$/,"",$2); print $2,$1}' | sort -k1 -n -r) \
        <(awk 'BEGIN{FS=OFS="\t"}{if(NR>1) print $1,$9,$7,$6,$5,$4,$3,$2}' taxonomy.tsv | sort -k1 -n -r) | \
        awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9" "$10}' | \
        cat <(echo -e "REF_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies") - > EMU_taxa.txt

        rm -r taxonomy.tsv

        grep -qF "export EMUDB=\"$PWD\"" ~/.bashrc || echo "export EMUDB=\"$PWD\"" >> ~/.bashrc

        source ~/.bashrc

        cd $base_dir

fi

cd $base_dir/DATA/MINIMAP2/EMUDB

grep -qF "export EMUDB=\"$PWD\"" ~/.bashrc || echo "export EMUDB=\"$PWD\"" >> ~/.bashrc

source ~/.bashrc

cd $base_dir