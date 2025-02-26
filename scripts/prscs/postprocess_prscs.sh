#!/bin/bash

#$ -l h_vmem=15G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y



y=$1



source /broad/software/scripts/useuse
use R-4.1


Rscript --vanilla ../scripts/prscs/postprocess_prscs.R $y



##EOF



