#!/usr/bin/bash

set -ex

if [ ! -d "atac-seq-pipeline" ]; then
wget https://github.com/ENCODE-DCC/atac-seq-pipeline/archive/refs/tags/v2.2.3.tar.gz
tar xvfz v2.2.3.tar.gz 
rm v2.2.3.tar.gz
mv atac-seq-pipeline-2.2.3 atac-seq-pipeline
cd atac-seq-pipeline/scripts
mkdir -p ../data
bash download_genome_data.sh hg38 ../data/hg38

conda env create -f atac-seq-pipeline.yaml
conda activate atac-seq-pipeline
caper init local
fi

if ! conda env list | grep -q "^QTLtools\s"; then
if [ -d "QTLtools" ]; then
cd QTLtools
conda build r-base
conda build qtltools --use-local
conda create -n QTLtools
conda install -n QTLtools qtltools --use-local -y
fi
fi

if [ ! -d "pre_imputation_check" ]; then
mkdir -p pre_imputation_check
cd pre_imputation_check
wget https://www.chg.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
unzip HRC-1000G-check-bim-v4.3.0.zip
wget https://www.chg.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz
gunzip 1000GP_Phase3_combined.legend.gz
fi

