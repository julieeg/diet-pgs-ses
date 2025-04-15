#!/bin/bash

#$ -l h_vmem=5G
#$ -l h_rt=0:30:00

#$ -j y
#$ -o reports/
#$ -cwd



CHR=$SGE_TASK_ID

PHENOS=$1
ANC=MA



gwas_dir=../data/processed/gwas
prscs_dir=../data/processed/prscs



opt=../../opt

source /broad/software/scripts/useuse
use R-4.1


regenie_file=${gwas_dir}/ukb_chr${CHR}_${ANC}_${PHENOS}_step2.regenie



#Note, regenie: ALLELE0=REF/ALLELE1=ALT ; PRScs: A1=GWAS allele/A2=GWAS reference

## Prepare bim file

echo Preparing bim file ...

R --vanilla << EOF
library(data.table) ; library(tidyverse) 

mfi=fread("/broad/ukbb/imputed_v3/ukb_mfi_chr${CHR}_v3.txt") %>% select(ID=V2,REF=V4, ALT=V5)
Ydict=fread("${regenie_file}.Ydict", header=F) %>% pull(V2) ; Ydict

ss=fread("${regenie_file}") %>% arrange(desc(A1FREQ)) %>%
mutate(CHR=ifelse(as.numeric(CHROM)<10, paste0("0", CHROM), paste(CHROM))) %>%
distinct(ID, .keep_all=T)

joined=left_join(ss, mfi) %>% arrange(desc(A1FREQ)) %>% distinct(ID, .keep_all=T) %>%
mutate(POS=0) %>%
select(CHR, ID, GENPOS, POS, REF, ALT) %>%
fwrite("${prscs_dir}/ukb_chr${CHR}_${ANC}_${PHENOS}.bim", sep="\t", row.names=F, col.names=F)


## Prepare summary_stats for EACH phenotype

cat("Preparing summary stats file for PRScs... \n ")

lapply(Ydict, function(y) { 
fread(paste0("${gwas_dir}/ukb_gwas_${ANC}_",y,"_regenie_merged")) %>%
select(SNP, A1=EA, A2=NEA, BETA, SE) %>%
fwrite(paste0("${prscs_dir}/ukb_chr${CHR}_MA_",y,"_regenie.prscsInput.txt"),sep="\t",row.names=F)
})

EOF


echo Done preparing files for PRScs!


#EOF

