#!/bin/bash

#$ -l h_vmem=5G
#$ -l h_rt=0:30:00

#$ -j y
#$ -o reports/
#$ -cwd



CHR=$SGE_TASK_ID



gwas_dir=../data/processed/gwas
prscs_dir=../data/processed/prscs



opt=../../opt

source /broad/software/scripts/useuse
use R-4.1



regenie_file=ukb_chr${CHR}_macros_MA_step2.regenie




## Prepare bim file

echo Preparing bim file ...

R --vanilla << EOF
library(tidyverse) ; library(data.table)
fread("${gwas_dir}/${regenie_file}") %>%
mutate(CHR=ifelse(as.numeric(CHROM)<10, paste0("0", CHROM), paste(CHROM))) %>%
mutate(ALT=ifelse(A1FREQ<0.5, ALLELE1, ALLELE0), REF=ifelse(A1FREQ<0.5, ALLELE0, ALLELE1)) %>%
mutate(POS=0) %>%
select(CHR, ID, GENPOS, POS, ALT, REF) %>%
fwrite("${prscs_dir}/ukb_chr${CHR}_MA_regenie.bim", sep="\t", row.names=F, col.names=F)
EOF



#Prepare summary_stats for EACH macronutrient (Y1=carb; Y2=fat; Y3=pro)

echo Preparing summary stats file ...

for y in {1..3} ; do
R --vanilla <<EOF
library(tidyverse) ; library(data.table)
fread("${gwas_dir}/${regenie_file}") %>%
select(SNP=ID, A1=ALLELE1, A2=ALLELE0, BETA=BETA.Y${y}, SE=SE.Y${y}) %>%
fwrite(paste0("${prscs_dir}/ukb_chr${CHR}_Y${y}_MA_regenie.txt"),sep="\t",row.names=F)
EOF
done


echo Done preparing files for PRScs!


#EOF

