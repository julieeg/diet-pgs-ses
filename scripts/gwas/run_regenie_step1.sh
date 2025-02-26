#!/bin/bash

#$ -l h_vmem=100G
#$ -l h_rt=36:00:00
#$ -o reports/

#$ -j y
#$ -cwd


ANC=$1



pheno_file=../data/processed/gwas/ukb_phenos_gwas_macros_${ANC}.txt
covar_file=../data/processed/gwas/ukb_phenos_gwas_covariates_${ANC}.txt
sample_file=../data/processed/gwas/ukb_phenos_gwas_macros_${ANC}_samples.txt


ukb_sample=/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample
ukb_bgen_dir=/broad/ukbb/imputed_v3 


scratch=/broad/hptmp/gervis2
gwas_dir=../data/processed/gwas
opt=../../opt



source /broad/software/scripts/useuse
use .regenie-3.2.2



##################
## regenie GWAS ##
##################

## run REGENIE Step 1 
regenie \
  --step 1 \
  --bgen ${scratch}/ukb_allchr_hapmap_step1.bgen \
  --ref-first \
  --sample ${ukb_sample} \
  --keep ${sample_file} \
  --force-step1 \
  --phenoFile ${pheno_file} \
  --covarFile ${covar_file} \
  --maxCatLevels 35 \
  --catCovarList geno_array,ac,smoke_current,smoke_former,alcohol,health \
  --bsize 1000 \
  --no-split \
  --lowmem \
  --lowmem-prefix ../data/temp \
  --out ../data/processed/gwas/ukb_gwas_macros_${ANC}_step1 \
  --threads 8


#EOF


