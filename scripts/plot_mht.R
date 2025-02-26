
#manhattans

library(data.table)
library(tidyverse)
library(qqman)
library(topr)



args <- commandArgs(trailingOnly = TRUE)
y <- args[1]
sst_file=paste0("../data/processed/gwas/ukb_gwas_Y",y,"_MA_regenie_merged")


sst <- fread(sst_file)
sst05 <- sst %>% filter(P<0.05)
y_label <- ifelse(y==1, "Carbohydrate (%)", ifelse(y==2, "Fat (%)", ifelse(y==3, "Protein (%)", NA)))




## Make simple mht plot with P<0.05 (qqman)
plot_mht_simple_05 <- function(sst_05) {
  qqman::manhattan(sst_05, chr = "CHR", bp = "POS", p = "P", snp = "SNP", 
                   main = paste("Manhattan Plot:", y_label, "(P <0.05)")) 
  }


## Make simple mht plot with gene annotations with P<0.05 
plot_mht_genes_05 <- function(sst05) {
  topr::manhattan(
    sst05, 
    annotate = 5e-08, 
    color = "black",
    ntop=4,
    sign_thresh_color = "darkred",
    sign_thresh_label_size = 2,
    label_fontface = "bold",
    title = paste("Manhattan Plot: ", y_label), 
    angle = 90,
    label_size = 3)
}



#Save plot as pdf
pdf(paste0("../data/processed/gwas/qq_plots/ukb_gwas_Y",y,"_mht_basic.pdf"), height=6, width=10)
plot_mht_simple_05(sst05)
dev.off()


#Save plot as pdf
pdf(paste0("../data/processed/gwas/qq_plots/ukb_gwas_Y",y,"_mht_basic_annotated.pdf"), height=6, width=10)
plot_mht_genes_05(sst05)
dev.off()



#sst5e5_genes <- sst05 %>% filter(P<5e-5) %>% rename(CHROM=CHR) %>% annotate_with_nearest_gene()
#sst5e5_annotated <- sst5e5_genes %>% filter(Gene_Symbol %in% sst5e5_genes$Gene_Symbol)

#pdf(paste0("../data/processed/gwas/qq_plots/ukb_gwas_Y",y,"_mht_annotated.pdf"), height=6, width=10)
#manhattan(list(sst05, sst5e5_annotated), 
#          color=c("darkgrey", c("#C6373C75")), 
#          annotate = c(5e-50), 
#          even_no_chr_lightness = c(0.8,0.5,0.5), size=1.2, 
#          sign_thresh_color = "black",
#          axis_text_size=8, axis_title_size=9,
#          label_fontface = "bold",
#          ntop = 2)
#dev.off()


#EOF






