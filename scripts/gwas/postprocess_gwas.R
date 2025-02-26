## regenie postprocess

library(tidyverse) 
library(data.table)


# set command args
args <- commandArgs(trailingOnly = TRUE)
y <- args[1]

gwas_dir="../data/processed/gwas"
prscs_dir="../data/processed/prscs"
pgs_dir="../data/processed/pgs"


Y=paste0("Y",y)


## Write functions (from kwesterman)

calc_lambda <- function(x, p=0.5){
  # Calculate genomic inflation lambda value
  x = x[!is.na(x)]
  x.quantile <- quantile(x, p)
  round(qchisq(1 - x.quantile, 1) / qchisq(p, 1), 2)
}

make_qq <- function(data=play, pval_col="p", main=""){
  # Make a quantile-quantile plot
  #data <- filter(data, data[[pval_col]] > 0)  # In case extremely low p-values are stored as zero
  
  # Process p-values
  y_vals <- sort(-log10(data[[pval_col]]))
  x_vals <- -log10(rev(ppoints(length(y_vals))))  # ppoints generates a uniform probability distribution
  
  # Trim points at higher p-values (credit to RaMWAS package for code snippet)
  levels = as.integer((x_vals - x_vals[1]) / (tail(x_vals, 1) - x_vals[1]) * 2000)
  keep = c(TRUE, diff(levels) != 0)
  levels = as.integer((y_vals - y_vals[1])/(tail(y_vals, 1) - y_vals[1]) * 2000)
  keep = keep | c(TRUE, diff(levels) != 0)
  keep = which(keep)
  
  par(ps = 18)
  plot(x = x_vals[keep], y = y_vals[keep], 
       xlab = expression(-log[10](italic(p)) * " (Expected)"), 
       ylab = expression(-log[10](italic(p)) * " (Observed)"),
       main = main, cex = 0.8, 
       cex.lab = 0.8, cex.main = 0.9, 
       pch = 16, ylim = c(0, ceiling(max(y_vals))))
  abline(0, 1, lty = 2)
  legend(x = 'topleft', y = 'topleft',
         bquote(lambda == .(calc_lambda(data[[pval_col]]))), 
         cex = 0.9, bty = "n")
}

sumstats_cols <- c(CHR="CHROM", POS="GENPOS", SNP="ID", EA="ALLELE1", NEA="ALLELE0",
                   EAF="A1FREQ", N="N", BETA=paste0("BETA.",Y), SE=paste0("SE.",Y),
                   P_LOG10=paste0("LOG10P.",Y))



## Combine regenie sumstats
cat("Gathering regenie summary stats files: ... \n")

filepaths <- grep("Ydict", list.files(gwas_dir, pattern="step2.regenie", full.names = T), invert=T, value=T)
sumstats <- reduce(lapply(filepaths, fread), bind_rows) %>%
  select(all_of(sumstats_cols)) %>%
  mutate(P=10^(-P_LOG10))



## Make Q-Q plot
cat("Making QQ plot ... \n")

qq_dir <- paste0(gwas_dir, "/qq_plots/")
system(paste0("mkdir -p ", qq_dir))

write(calc_lambda(sumstats$P), paste0(qq_dir, paste0("ukb_gwas_",Y, "_MA_lambda")))
plot_filepath <- paste0(qq_dir, paste0("ukb_gwas_",Y, "_MA_qq.pdf"))
pdf(file = plot_filepath)
make_qq(sumstats, "P")
dev.off()



## Save merged sumstats file
cat("Saving merged summary stats file ... \n")

sumstats %>%
  arrange(P) %>%
  fwrite(paste0(gwas_dir, "/ukb_gwas_",Y, "_MA_regenie_merged"), sep="\t")


#EOF

