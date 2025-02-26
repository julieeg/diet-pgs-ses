#!/bin/bash
#$ -l h_vmem=30G
#$ -l h_rt=0:10:00
#$ -o reports/prep_ukb_phenos.report

#$ -cwd
#$ -j y

> reports/prep_ukb_phenos.report


source /broad/software/scripts/useuse
use R-4.1


Rscript --no-save ../scripts/data_prep/prep_ukb_phenos.R


##END_OF_SYNTAX
