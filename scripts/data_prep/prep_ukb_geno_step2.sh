#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=2:00:00
#$ -o reports/

#$ -j y
#$ -cwd



CHR=$SGE_TASK_ID



ukb_sample=/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample
ukb_bgen_dir=/broad/ukbb/imputed_v3
gwas_dir=../data/processed/gwas
scratch=/broad/hptmp/gervis2


source /broad/software/scripts/useuse
opt=../../opt

reuse -q Anaconda3
source activate $opt/bgen




################################################################
## Simple genotyping QC (maf>0.001  & info >0.4) for ukb GWAS ##
################################################################

# Filter variants using mfi file

awk -v CHR=$CHR '{ if ($6 > 0.001 && $8 >0.4 ) { print $2 } }' ${ukb_bgen_dir}/ukb_mfi_chr${CHR}_v3.txt > ${scratch}/chr${CHR}_maf001_imp04.snplist



# Create filtered bgen file

bgenix -g ${ukb_bgen_dir}/ -incl-rsids ${scratch}/chr${CHR}_maf001_imp04.snplist > ${scratch}/chr${CHR}_maf001_imp04.bgen
bgenix -g ${scratch}/chr${CHR}_maf001_imp04.bgen index -clobber



##EOF

