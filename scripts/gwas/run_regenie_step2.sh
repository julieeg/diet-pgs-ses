#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=6:00:00
#$ -o reports/

#$ -j y
#$ -cwd



CHR=$SGE_TASK_ID


ANC=$1




pheno_file=../data/processed/gwas/ukb_phenos_gwas_macros_${ANC}.txt
covar_file=../data/processed/gwas/ukb_phenos_gwas_covariates_${ANC}.txt
sample_file=../data/processed/gwas/ukb_phenos_gwas_macros_${ANC}_samples.txt


ukb_sample=/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample
ukb_bgen_dir=/broad/ukbb/imputed_v3 


scratch=/broad/hptmp/gervis
gwas_dir=../data/processed/gwas
opt=../../opt



source /broad/software/scripts/useuse
use .regenie-3.2.2





##################
## regenie GWAS ##
##################

## run REGENIE Step 1 
regenie \
  --step 2 \
  --bgen ${scratch}/chr${CHR}_maf001_imp04.bgen \
  --ref-first \
  --sample ${ukb_sample} \
  --keep ${sample_file} \
  --phenoFile ${pheno_file} \
  --covarFile ${covar_file} \
  --catCovarList geno_array,ac,smoke_current,smoke_former,alcohol,health \
  --maxCatLevels 35 \
  --test additive \
  --bsize 1000 \
  --no-split \
  --lowmem --lowmem-prefix ../data/temp/ukb_chr${CHR}_macros_MA_step2 \
  --pred ${gwas_dir}/ukb_gwas_macros_MA_step1_pred.list \
  --out ../data/processed/gwas/ukb_chr${CHR}_macros_MA_step2 \
  --threads 8


#EOF


