#!/bin/bash

#$ -l h_vmem=30G
#$ -l h_rt=0:30:00
#$ -o reports/

#$ -j y
#$ -cwd


ANC=$1



ukb_sample=/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample
ukb_bgen_dir=/broad/ukbb/imputed_v3
gwas_dir=../data/processed/gwas


scratch=/broad/hptmp/gervis2
opt=../../opt


source /broad/software/scripts/useuse
use Anaconda3
source activate $opt/bgen




############################################
## Create combined variant list for step1 ## 
############################################

## Write empty file & print chr* with for loop

> ${gwas_dir}/ukb_allchr_hapmap_${ANC}_step1.prune.in

for i in {1..22} ; do 
cat ${gwas_dir}/ukb_chr${i}_hapmap_${ANC}_step1.prune.in | awk '!seen[$0]++' >> $gwas_dir/ukb_allchr_hapmap_${ANC}_step1.prune.in
done



## Combine chr-stratified & qc-ed bgen files

cat-bgen -g ${scratch}/chr*_hapmap_${ANC}_step1.bgen -og ${scratch}/ukb_allchr_hapmap_${ANC}_step1.bgen -clobber
bgenix -g ${scratch}/ukb_allchr_hapmap_${ANC}_step1.bgen -index -clobber



## Delete breadcrumbs

# chr bgen files
for i in {1..22} ; do rm ${scratch}/chr${i}_hapmap_${ANC}_step1.bgen* ; done 

# chr snplists
rm ${gwas_dir}/ukb_chr*_hapmap_${ANC}_step1.snplist



##EOF

