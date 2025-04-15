# GWAS comparison (SAIGE vs. REGENIE)

## load packages & pantry scripts
library(tidyverse) 
library(data.table)

lapply(list.files("../../pantry/functions/", full.names = T), source)



#####################################
## Set-up variable arrays & themes ##
#####################################

# set color palettes
palette_macros <- c("carb"=palettes$NatExt$Oranges[2], "fat"=palettes$NatExt$Blues[2], "prot"= palettes$NatExt$Greens[2])
macros.labs <- c(carb="Carbohydrate (%kcal)", fat="Fat (%kcal)", prot="Protein (%kcal)")
macros_pgs.labs <- c(carb="Carbohydrate PGS", fat="Fat PGS", prot="Protein PGS")

# set variable arrays
macro_phenos <- c("carb_pheno", "fat_pheno", "prot_pheno")
names(macro_phenos) <- macros.labs

macro_pgs <- c("carb_pgs", "fat_pgs", "prot_pgs")
names(macro_pgs) <- macros_pgs.labs

# set plot themes
plot_theme <- theme_bw() + theme(
  axis.text = element_text(size=6), 
  axis.title = element_text(size=8)
) 


## Load analysis data ----------------------

analysis <- readRDS("../data/processed/phenos/ukb_analysis_MA.rda")



##################################################
## Histograms of macronutrient phenotypes & PGS ##
##################################################

# macronutrient phenotypes ----------------------

plot_macro_hists.l <- lapply(1:3, function(m) {
  plot_continuous(macro_phenos[[m]], data=analysis, fill=palette_macros[m]) + ggtitle(" ") + ylim(0,31000) +
    xlab(names(macro_phenos[m])) + plot_theme
}) ; ggarrange(plotlist = plot_macro_hists.l, nrow=1) %>% 
ggsave(filename="../output/plot_macro_pheno_histograms.pdf", height=3.5, width=8)


# macronutrient pgs ----------------------------

plot_pgs_hists.l <- lapply(1:3, function(m) {
  plot_continuous(macro_pgs[[m]], data=analysis, fill=palette_macros[m]) + ggtitle(names(macro_pgs[m])) + 
    xlab(names(macro_pgs[m])) + ylim(0,31000) + plot_theme + ggtitle(" ") #theme(axis.text.y = element_text(angle=90, hjust=0.5, vjust=0)) + 
}) ; ggarrange(plotlist = plot_pgs_hists.l, nrow=1) %>% 
ggsave(filename="../output/plot_macro_pgs_histograms.pdf", height=3.5, width=8)



####################################################
## Comparison of SAIGE (JM) vs REGENIE (JG) GWAS  ##
####################################################

library(qqman)
library(topr)

# set column names for gwas summary stats
sumstats_cols <- c(
  CHR="CHROM", POS="GENPOS", SNP="ID", EA="ALLELE1", NEA="ALLELE0", EAF="A1FREQ", 
  N="N", BETA=paste0("BETA.",Y), SE=paste0("SE.",Y), P_LOG10=paste0("LOG10P.",Y))

# simple manhattan plot function
plot_mht <- function(sumstats, y_label,...) {
  qqman::manhattan(sumstats, chr = "CHR", bp = "POS", p = "P", snp = "SNP", 
                   main = y_label, ...) 
}


## Load & format gwas summary stats ------------------------

# Load & format JM gwas (chr21)
cols <- c("CHR", "POS", SNP="rsid", A1="Allele1", A2="Allele2", A2FREQ="AF_Allele2", "BETA", "SE", "zscore", P="p.value")
chr21_jm.l <- list(
  carb=fread("../data/sumstats/CARB_chr21_UKBB_results.txt") %>% select(all_of(cols)), 
  fat=fread("../data/sumstats/FAT_chr21_UKBB_results.txt") %>% select(all_of(cols)),
  prot=fread("../data/sumstats/PROT_chr21_UKBB_results.txt") %>% select(all_of(cols))
)


# Load & format JG gwas (chr21)
chr21_jg <- fread("../data/processed/gwas/test_ukb_chr21_macros_EUR_step2.regenie")
cols.l <- lapply(1:3, function(m) {
  c(CHR="CHROM", POS="GENPOS", SNP="ID", A1="ALLELE0", A2="ALLELE1",
    A2FREQ="A1FREQ", BETA=paste0("BETA.Y",m), SE=paste0("SE.Y",m), LOGP=paste0("LOG10P.Y",m))
}) ; chr21_jg.l <- lapply(1:3, function(m) {chr21_jg %>% select(cols.l[[m]]) %>% 
    mutate(zscore=BETA/SE) %>% mutate(P=10^(-LOGP))})
names(chr21_jg.l) <- c("carb", "fat", "prot")


# Remove duplicate (multi-allelioc) SNPs (prioritizing smaller P-value)
remove_dupl_snps <- function(sumstats) {
  dupl <- which(sumstats %>% pull(SNP) %>% duplicated())
  sumstats <- sumstats[-dupl] 
  return(sumstats)
}

chr21_jm_filtered.l <- lapply(chr21_jm.l, remove_dupl_snps)
chr21_jg_filtered.l <- lapply(chr21_jg.l, remove_dupl_snps)



## Pearson correlations & dot plots --------------------------------

plot_chr21_mthd.l <- lapply(c("carb", "fat", "prot"), function(macro) {
  
  temp <- left_join(chr21_jm_filtered.l[[macro]], chr21_jg_filtered.l[[macro]], by = "SNP") %>%
    mutate(zscore.y=ifelse(A2.x==A2.y,zscore.y,-zscore.y)) %>% 
    select(SNP, saige=zscore.x, regenie=zscore.y)
  
  cor.mthd <- cor.test(temp$saige, temp$regenie)
  cor_mthd <- sprintf("%s (%s, %s)", round(cor.mthd$estimate,4), 
                      round(cor.mthd$conf.int[1],4), round(cor.mthd$conf.int[2],4))
  
  temp %>% 
    ggplot(aes(x=saige, y=regenie)) + plot_theme +
    geom_point(alpha=0.55, color=palette_macros[[macro]]) +
    geom_abline(aes(slope=1, intercept=0)) + 
    xlab("SAIGE, beta/SE") + ylab("REGENIE, beta/SE") + 
    ggtitle(paste0(macros.labs[[macro]], ", chr21: cor = ", round(cor.mthd$estimate,2))) + 
    theme(axis.text = element_text(size=7),
          plot.title = element_text(face="bold", size=8))
  
}) ; ggarrange(plotlist = plot_chr21_mthd.l, nrow=1) %>% 
  ggsave(filename="../output/plot_compare_gwas_byMethod_macros_chr21.pdf", height=3.5, width=8.5)


## Plot mini manhattan plots for chr 21 -----------------------------
library(gridExtra)
lapply(1:3, function(m) {
  pdf(paste0("../output/plot_compare_gwas_manhattan_byMethod_", macro_pgs[m], "_macros_chr21.pdf"), height=3.5, width=8)
  plot_mht(chr21_jm_filtered.l[[m]], col=palette_macros[m], y_label = paste0("SAIGE: ", macros_pgs.labs[m], ", EUR only"))
  plot_mht(chr21_jg_filtered.l[[m]], col=palette_macros[m], y_label = paste0("REGENIE: ", macros_pgs.labs[m], ", EUR only"))
  dev.off()
  })




############################################
## Comparison of REGENIE EUR vs M-A GWAS  ##
############################################

## Pearson correlations & dot plots --------------------------------

chr21_ma=fread("../data/processed/gwas/ukb_chr21_macros_MA_step2.regenie")
chr21_ma.l <- lapply(1:3, function(m) {
  remove_dupl_snps( chr21_ma %>% select(cols.l[[m]]) %>% 
                      mutate(zscore=BETA/SE) %>% mutate(P=10^(-LOGP)) )
}) ; names(chr21_ma.l) <- c("carb", "fat", "prot")


plot_chr21_anc.l <- lapply(c("carb", "fat", "prot"), function(macro) {
  
  temp <- left_join(chr21_jg_filtered.l[[macro]], chr21_ma.l[[macro]], by = "SNP") %>%
    mutate(zscore.y=ifelse(A2.x==A2.y,zscore.y,-zscore.y)) %>% 
    select(SNP, EUR=zscore.x, MA=zscore.y)
  
  cor_anc <- round(cor.test(temp$EUR, temp$MA)$estimate, 2)
  
  temp %>% 
    ggplot(aes(x=EUR, y=MA)) + plot_theme +
    geom_point(alpha=0.55, color=palette_macros[[macro]]) +
    geom_abline(aes(slope=1, intercept=0)) + 
    xlab("European only, beta/SE (N=190801)") + ylab("Multi-ancetry, beta/SE (N=192402)") + 
    ggtitle(paste0(macros_pgs.labs[[macro]], ", chr21: cor = ", cor_anc)) + 
    theme(axis.text = element_text(size=7),
          plot.title = element_text(face="bold", size=8))
  
}) ; ggarrange(plotlist = plot_chr21_mthd.l, nrow=1) %>% 
  ggsave(filename="../output/plot_compare_gwas_byAnc_macros_chr21.pdf", height=3.5, width=8.5)


## EOF 
# Last updated: 0-20-2025







