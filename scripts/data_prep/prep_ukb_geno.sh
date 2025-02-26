#!/bin/bash

#$ -l h_vmem=55G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -j y
#$ -cwd



CHR=$SGE_TASK_ID

ANC=$1



ukb_sample=/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample
ukb_bgen_dir=/broad/ukbb/imputed_v3
gwas_dir=../data/processed/gwas
scratch=/broad/hptmp/gervis2


source /broad/software/scripts/useuse
opt=../../opt
reuse -q Anaconda3



## Compile list of common variants from HapMap3, with imputation INFO>0.8 and maf>0.01 & mac>200 among samples with complete pheno data

# Filter variants based on INFO (mfi) using R
use R-4.1
R --vanilla <<EOF
library(tidyverse);library(data.table)
highqual=fread("${ukb_bgen_dir}/ukb_mfi_chr${CHR}_v3.txt") %>% filter(V8>0.8) %>% pull(V2)
fread("hapmap_b36_ma.snplist", header=F) %>% filter(V1 %in% highqual) %>% fwrite("${gwas_dir}/ukb_chr${CHR}_hapmap.snplist")
EOF



# Create filtered bgen file to HapMap3 common variants with INFO>0.8

source activate $opt/bgen
bgenix -g ${ukb_bgen_dir}/ukb_imp_chr${CHR}_v3.bgen -incl-rsids ${gwas_dir}/ukb_chr${CHR}_hapmap.snplist > ${scratch}/chr${CHR}_hapmap.bgen
bgenix -g ${scratch}/chr${CHR}_hapmap.bgen -index -clobber



# Further restrict snplist to maf>0.01, mac>200, geno>0.1, mind>0.1, hwe p>1e-15 & independent variants (1000kb window; 100 snp step-size, LDr2>0.9)
# Following guidance from regenie paper for UKB genetic data qc for step1: https://www.nature.com/articles/s41588-021-00870-7#Sec23

$opt/plink2 \
--bgen ${scratch}/chr${CHR}_hapmap.bgen ref-first \
--sample ${ukb_sample} \
--keep ${gwas_dir}/ukb_phenos_gwas_macros_${ANC}_samples.txt \
--maf 0.01 --mac 200 --geno 0.1 \
--hwe 1e-15 \
--mind 0.1 \
--indep-pairwise 1000 100 0.9 \
--write-snplist \
--rm-dup force-first \
--memory 50000 \
--out ${gwas_dir}/ukb_chr${CHR}_hapmap_${ANC}_step1


# Filter bgen file to variants passing regenie step1 qc 
bgenix -g ${scratch}/chr${CHR}_hapmap.bgen -incl-rsids ${gwas_dir}/ukb_chr${CHR}_hapmap_${ANC}_step1.prune.in > ${scratch}/chr${CHR}_hapmap_step1.bgen
bgenix -g ${scratch}/chr${CHR}_hapmap_step1.bgen -index -clobber
rm ${scratch}/chr${CHR}_hapmap.bgen*


##EOF



