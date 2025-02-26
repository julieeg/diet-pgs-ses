## Script to prepare macronutrient phenotype files for UKBB GWAS, mirroring phenotypes 
## prepared for the Merino et al GWAS of macronutrient preference

library(data.table)
library(tidyverse)


#######################
## Load raw datasets ##
#######################

# open main, alch & diet UKBB phenotype files 
cat("Reading in phenotype files ... \n")
main=fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz",
           data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=f.eid,
         ac=f.54.0.0,
         sex=f.31.0.0, 
        # ethnicity=f.21000.0.0, 
         age=f.21022.0.0, 
         bmi=f.21001.0.0, 
         age2=f.21003.0.0,
         smoke=f.20116.0.0,
         health=f.2178.0.0,
         townsend=f.189.0.0) %>% 
  mutate(sex.lab = ifelse(sex==0, "Female", "Male")) #n=502618


# alcohol use
alch_id=fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb669173.tab.gz",
                 data.table=FALSE, stringsAsFactors = FALSE) %>%
  select(id=f.eid, alcohol=f.1558.0.0)


# macronutrient intake
diet=fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb22861.tab.gz", 
           data.table=FALSE, stringsAsFactors = FALSE) %>%
  select(id=f.eid, starts_with("f.100002."), 
         starts_with("f.100005."), 
         starts_with("f.100003."), 
         starts_with("f.100004.")) %>%
  # rename with nutrient labels
  rename_with(., ~gsub("f.100002.", "energy", gsub("f.100005.", "carb", 
                       gsub("f.100003.", "prot", gsub("f.100004.", "fat", .))))) %>% 
  rename_with(., ~gsub("[.]0", "", .)) 


# assessment center
ac_labs <- c("11012"="Barts", "11021"="Birmingham", "11011"="Bristol", "11008"="Bury", "11003"="Cardiff", 
	     "11024"="Cheadle.(revisit)", "11020"= "Croydon", "11005"="Edinburgh", "11004"="Glasgow", "11018"="Hounslow",
	     "11010"="Leeds", "11016"="Liverpool", "11001"="Manchester", "11017"="Middlesborough", "11009"="Newcastle", 
	     "11013"="Nottingham", "11002"="Oxford", "11007"="Reading", "11014"="Sheffield", "10003"="Stockport.(pilot)",
	     "11006"="Stoke", "11022"="Swansea", "11023"="Wrexham", "11025"="Cheadle.(imaging)",
	     "11026"="Reading.(imaging)", "11027"="Newcastle.(imaging)", "11028"="Bristol.(imaging)")

main$ac <- ac_labs[as.character(main$ac)]
table(main$ac)


## Get ids for participants who withdrew consent
withdraw=scan("/humgen/florezlab/UKBB_app27892/withdraw/w27892_20241217.csv", what=character()) #488


## merge phenotype files & remove participants that withdrew consent
dat = full_join(main, diet, by = "id") %>%
  left_join(alch_id, by = "id") %>%
  filter(!id %in% withdraw) 

#N=502130



################################
## Create phenotypes for GWAS ##
################################
cat("Preparing macronutrient & covariate phenotypes....\n")

# winsorize
winsorize <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x<bounds[1], bounds[1], ifelse(x>bounds[2], bounds[2], x))
  x
}


## Prepare lifestyle phenotypes
dat <- dat %>% 
  ## Convert energy (kj -> kcal) & nutrients (g -> kcal) intakes
  mutate(across(starts_with("energy"), ~./4.184)) %>%
  mutate(across(starts_with("carb"), ~.*4)) %>%
  mutate(across(starts_with("prot"), ~.*4)) %>%
  mutate(across(starts_with("fat"), ~.*9)) %>%

  # Prepare covariates
  mutate(across(c("smoke", "health", "alcohol"), ~ifelse(. <0, NA, .))) %>%
  mutate(
    #Coding: smoke == 0 ~ "Never", smoke == 1 ~ "Previous", smoke == 2 ~ "Current")),
    smoke_former = case_when(smoke == 0 ~ 0, smoke == 1 ~ 1, smoke == 2 ~ 0),
    smoke_current = case_when(smoke == 0 ~ 0, smoke == 1 ~ 0, smoke == 2 ~ 1),
    #Coding: alcohol = (alcohol == 6 ~ "Never", alcohol == 5 ~ "Special occasions only",
    #  alcohol == 4 ~ "One to three times a month", alcohol == 3 ~ "Once or twice a week",
    #  alcohol == 2 ~ "Three or four times a week", alcohol == 1 ~ "Daily or almost daily")),
    alcohol = case_when(
      alcohol == 6 ~ 0, alcohol == 5 ~ 1, alcohol == 4 ~ 2, alcohol == 3 ~ 3,
      alcohol == 2 ~ 4, alcohol == 1 ~ 5))
    #Coding: health = as.factor(case_when(
    #  health == 1 ~ "Poor", health == 2 ~ "Fair", health == 3 ~ "Good", health == 4 ~ "Excellent"))) %>%
  

# Prepare macronutrient phenotypes

valid_macro_mean <- function(macro, df) {
  # For a given nutrient:
  # - Select all columns for that nutrient
  # - Tabulate number of instances (not missing & with valid total energy)
  # - calculate mean over all valid incidences; if not, code as NA
  diet_fields_df <- lapply(0:4, function(i) {
    macro_id <- paste0(macro, i) # diet field 
    energy_field <- paste0("energy", i) # ukb field for valid 24hr (energy & sex)
    valid_id <- df %>% select(id, macro=macro_id, energy=energy_field, sex=sex.lab) %>%
      mutate(macro=ifelse(energy<500 & sex=="Female" | energy>3500 & sex=="Female" |  
                            energy<800 & sex=="Male" | energy>4200 & sex=="Male", NA, macro)) %>% # replace with NA, if invalid macro
      select(id, macro_id=macro) %>% rename_with(., ~gsub("macro_id", macro_id, .))
    valid_id
  }) %>% reduce(full_join, by = "id") %>%
    mutate(n_valid_id = rowSums(!is.na(across(paste0(macro,0:4))))) %>%  # count n of valid 24hr
    mutate("macro_mean"=ifelse(n_valid_id>0, rowSums(across(paste0(macro,0:4)), na.rm=T)/n_valid_id, NA)) %>%
    select(id, macro_mean) %>%
    rename_with(., ~gsub("macro_mean", paste0(macro, "_pheno"), .)) %>%
    pull(ends_with("pheno"))
}


macros <- c("energy", "carb", "fat", "prot")
macro_phenos <- lapply(macros, function(m) {
  dat %>% mutate(
    macro_pheno=valid_macro_mean(m, .)) %>%
    rename_with(., ~gsub("macro", m, .)) %>%
    select(id, ends_with("pheno"))
}) %>% reduce(full_join, by = "id")


dat <- dat %>% left_join(macro_phenos, by = "id") %>%
  # Select phenotypes
  select(id, age, ac, sex, bmi, age2, 
         energy_pheno, carb_pheno, fat_pheno, prot_pheno,
         smoke_former, smoke_current, alcohol, health, townsend)

dim(dat)
head(dat)


###################################
## Run genotying & genetic PC QC ##
###################################
cat("Running genotype QC: restricting to european genetic ancestry ....\n")

# function to identify outliers
find_outliers.fun <- function(x, SDs=6) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- as.factor(ifelse(x<bounds[1] | x>bounds[2], 1, NA))
  x
}

## open UKBB QC file (get PCs and genotyping array info)
geno_qc <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb669177.tab.gz",
            data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=f.eid,
         het_miss_outliers=f.22027.0.0,
         sex_chrom_aneuploidy=f.22019.0.0,
         reported_sex=f.31.0.0,
         genetic_sex=f.22001.0.0,
         ethnicity=f.21000.0.0,
         geno_array=f.22000.0.0) %>%
  mutate_at("geno_array", ~ifelse(. <0, "UKBL", "UKBB"))

geno_ancestry <- fread("/humgen/florezlab/UKBB_app27892/ukbreturn2442/all_pops_non_eur_pruned_within_pop_pc_covs_app27892.csv",
                       data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id = f.eid, ancestry=pop_return2442) 

geno_pcs <- fread("/humgen/florezlab/users/pschroeder/UKBB/TOPMedImputation/GWAS_UKBB_UKBL/analyze_with_regenie/covariates_all_with_PRS",,
                  data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=IID, paste0("PC",1:20))


qc <- full_join(geno_qc, geno_ancestry, by = "id") %>%
  full_join(geno_pcs, by="id")

#Note: 


## Run QC
qc1 <- qc %>% mutate(
  # Exclude participants not passing genetic QC (missing values indicate PASS)
  qc_exclude = ifelse(het_miss_outliers == 1  | sex_chrom_aneuploidy == 1 | 
                        reported_sex != genetic_sex, 1, 0)) %>% 
  filter(is.na(qc_exclude)) %>% 
  
  # Add indicator for EUR ethnicity (matching Merino et al.)
  mutate_at("ethnicity", ~case_when(
    . == 1001 ~ "British", . == 1 ~ "White", . == 1002 ~ "Irish", 
    . == 1003 ~ "Any other white background")) %>%
  mutate(eur = ifelse(ethnicity == "British" | ethnicity == "White" | 
                        ethnicity == "Irish" | ethnicity == "Any other white background", 1, 0))
#N=500701 


## Add indicator to exclude outliers at +/- 6SD from the mean based on the first 4 PCs (use PanUKB)
geno_id <- qc1 %>%
  mutate(PC_outliers = ifelse(find_outliers.fun(PC1) == 1 | find_outliers.fun(PC2) == 1 | 
                                find_outliers.fun(PC3) == 1 | find_outliers.fun(PC4) == 1, 1, NA)) %>%
  #filter(is.na(PC_outliers)) %>%
  select(id, starts_with("PC"), ethnicity, ancestry, geno_array, PC_outliers, eur)
 


#######################################
## Create phenotype dataset for GWAS ##
#######################################
cat("Writing .txt files for analysis .... \n")

processed <- dat %>% filter(id %in% geno_id$id) %>%
  left_join(geno_id, by = "id") %>%
  mutate(FID=id, IID=id, .before=id) %>%
  distinct(IID, .keep_all=T)

head(processed)

# ==================================
## EUR only (for JM replication)
# ==================================

# Winsorize macronutrients (5 SD)
eur_processed <- processed %>%
  filter(eur == 1 & is.na(PC_outliers)) %>%
  select(-PC_outliers) %>%
  mutate(across(c("carb_pheno", "prot_pheno", "fat_pheno"), ~winsorize(.)))

# all phenos
eur_processed %>%
  fwrite("../data/processed/ukb_phenos_gwas_EUR.txt", sep="\t", na="NA") #N=467490
eur_processed <- eur_processed %>% select(-c(ethnicity, ancestry, eur))

# macronutrient phenos
eur_processed %>%
  select(FID, IID, carb_pheno, fat_pheno, prot_pheno) %>%
  filter(complete.cases(.)==T) %>%
  fwrite("../data/processed/gwas/ukb_phenos_gwas_macros_EUR.txt", sep="\t", na="NA") #196434

# covariate phenos
eur_processed %>%
  select(FID, IID, geno_array, ac, age, sex, bmi, paste0("PC", 1:20), smoke_former, smoke_current, alcohol, health, townsend) %>%
  filter(complete.cases(.) == T) %>%
  fwrite("../data/processed/gwas/ukb_phenos_gwas_covariates_EUR.txt", sep="\t", na="NA")

# sample file (with complete data)
eur_processed %>%
  select(FID, IID, carb_pheno, fat_pheno, prot_pheno, geno_array, ac, age, sex, bmi, 
         paste0("PC", 1:20), smoke_former, smoke_current, alcohol, health, townsend) %>%
  filter(complete.cases(.)) %>% select(FID, IID) %>%
  unique(.) %>%
  fwrite("../data/processed/gwas/ukb_phenos_gwas_macros_EUR_samples.txt", sep="\t")


# ========================================
## MA for primary GWAS & PGS generation
# ========================================

# Winsorize macronutrients (5 SD)
ma_processed <- processed %>%
  select(-PC_outliers, -eur) %>%
  mutate(across(c("carb_pheno", "prot_pheno", "fat_pheno"), ~winsorize(.)))

# all phenos
ma_processed %>%
  fwrite("../data/processed/ukb_phenos_gwas_MA.txt", sep="\t", na="NA")

# macronutrient phenos
ma_processed %>%
  select(FID, IID, carb_pheno, fat_pheno, prot_pheno) %>%
  filter(complete.cases(.)==T) %>%
  fwrite("../data/processed/gwas/ukb_phenos_gwas_macros_MA.txt", sep="\t", na="NA")

# covariate phenos
ma_processed %>% 
  select(FID, IID, geno_array, ac, age, sex, bmi, paste0("PC", 1:20), smoke_former, smoke_current, alcohol, health, townsend) %>%
  filter(complete.cases(.) == T) %>%
  fwrite("../data/processed/gwas/ukb_phenos_gwas_covariates_MA.txt", sep="\t", na="NA")


# sample file (with complete data)
ma_processed %>%
  filter(complete.cases(.)) %>% select(FID, IID) %>%
  unique(.) %>%
  fwrite("../data/processed/gwas/ukb_phenos_gwas_macros_MA_samples.txt", sep="\t")



cat("DONE preparing ukb phenotype files for REGENIE GWAS")

# Single-trait genetic association analyses were performed separately for carbohydrate, 
# fat and protein as percentages of total energy in 191,157 UKBB participants using a 
# generalized mixed model implemented in SAIGE. (Scalable and Accurate Implementation 
# of GEneralized mixed model) v.0.35.8.8 (https://github.com/weizhouUMICH/SAIGE/).
# Models were adjusted for age, sex, BMI, 20 principal components of ancestry, 
# genotyping array and accounted for sample relatedness. We used a minor allele 
# count threshold of 20 to select the genetic variants. 



##EOF


