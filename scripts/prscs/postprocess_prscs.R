## regenie postprocess

library(tidyverse) 
library(data.table)


# set command args
args <- commandArgs(trailingOnly = TRUE)
y <- args[1]

prscs_dir="../data/processed/prscs"


Y=paste0("Y",y)



# Set summary stats columns
sumstats_cols <- c(CHR="V1", SNP="V2", POS="V3", EA="V4", NEA="V5", BETA="V6")

## Combine prscs sumstats
cat("Gathering prscs summary stats files: ... \n")

filepaths <- grep(Y, list.files(prscs_dir, pattern="posterior", full.names = T), value=T)
sumstats <- reduce(lapply(filepaths, fread), bind_rows) %>%
  select(all_of(sumstats_cols)) %>%
  arrange(CHR)


## Save merged sumstats file
cat("Saving merged summary stats file ... \n")

sumstats %>%
  fwrite(paste0(prscs_dir, "/ukb_gwas_",Y, "_MA_prscs_merged"), sep="\t")


#EOF

