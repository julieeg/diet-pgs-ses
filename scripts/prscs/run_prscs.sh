#!/bin/bash

#$ -l h_vmem=30G
#$ -l h_rt=5:00:00

#$ -j y
#$ -o reports/
#$ -cwd



CHR=$SGE_TASK_ID
y=$1



prscs_dir=../data/processed/prscs
opt=../../opt
PRScs=$opt/PRScs/PRScs.py


bim_file=ukb_chr${CHR}_MA_regenie
ss_file=ukb_chr${CHR}_Y${y}_MA_regenie.txt


source /broad/software/scripts/useuse
reuse -q Anaconda3



## Run PRScs (using EUR reference panel)

${PRScs} \
 --ref_dir=${opt}/PRScs/ldblk_1kg_eur \
 --bim_prefix=${prscs_dir}/${bim_file} \
 --sst_file=${prscs_dir}/${ss_file} \
 --n_gwas=108038 \
 --phi=1e-2 \
 --chrom=${CHR} \
 --out_dir=${prscs_dir}/ukb_chr${CHR}_Y${y}_MA_regenie 


mv ${prscs_dir}/ukb_chr${CHR}_Y${y}_MA_regenie_pst_eff_a1_b0.5_phi1e-02_chr${CHR}.txt ${prscs_dir}/ukb_chr${CHR}_Y${y}_MA_regenie.posterior

echo Done runnings PRScs!


#EOF

