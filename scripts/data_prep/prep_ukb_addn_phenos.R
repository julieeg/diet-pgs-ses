## Script to prepare macronutrient phenotype files for UKBB GWAS, mirroring phenotypes 
## prepared for the Merino et al GWAS of macronutrient preference

library(data.table)
library(tidyverse)


# winsorize
winsorize <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x<bounds[1], bounds[1], ifelse(x>bounds[2], bounds[2], x))
  x
}


################################################
## Compile additional phenotypes for analysis ##
################################################

cat("Preparing macronutrient & covariate phenotypes....\n")

dat <- fread("../data/processed/phenos/ukb_phenos_")


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
    #derive mean macronutrient intake variable (kcal)
    macro_pheno=valid_macro_mean(m, .)) %>%
    rename_with(., ~gsub("macro", m, .)) %>%
    select(id, ends_with("pheno"))
}) %>% reduce(full_join, by = "id")

macro_kcal <- macro_phenos

macro_phenos <- macro_phenos %>% mutate(
  carb_pheno = (carb_pheno/energy_pheno)*100,
  fat_pheno = (fat_pheno/energy_pheno)*100,
  prot_pheno = (prot_pheno/energy_pheno)*100
)

dat <- dat %>% left_join(macro_phenos, by = "id") %>%
  # Select phenotypes
  select(id, age, ac, sex, bmi, age2, 
         energy_pheno, carb_pheno, fat_pheno, prot_pheno, 
         smoke_former, smoke_current, alcohol, health, townsend, sleep)

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
 


###########################################
## Prepare other phenotypes for analysis ##
###########################################

## Education level ------------------------------

# coding based on: Ge T., et al. Cerebral Cortex 2019;29(8): 3471-3481.
educ_level_labs <- list(
  "None of the above" = -7, 
  "Prefer not to answer"= -3, 
  "College or university degree" = 1, 
  "A/AS levels or equivalent" = 2, 
  "O/GCSE levels or equivalent" = 3, 
  "CSEs or equivalent" = 4,
  "NVQ/HND or equivalent" = 5, 
  "Other professional qualifications" = 6)

educ_isced_level_labs <- list("Level 5" = 1, "Level 3" = 2, "Level 2" = 3, 
                              "Level 2" = 4, "Level 5" = 5, "Level 4" = 6, "Level 1" = -7)  # NA = -3 or missing

educ_years_labs <- list("20"=1, "13"=2, "10"=3, "10"=4, "19"=5, "15"=6, "7"=-7)


educ_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672670.tab.gz",
                 data.table = FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, educ_level = f.6138.0.0) %>%
  mutate(educ_level.lab = case_when(
    educ_level == -7 ~ "None of the above", educ_level == -3 ~ "Prefer not to answer",
    educ_level == 1 ~ "College or university degree", educ_level == 2 ~ "A/AS levels or equivalent",
    educ_lelve == 3 ~ "O/GCSE levels or equivalent", educ_level == 4 ~ "CSEs or equivalent",
    educ_level == 5 ~ "NVQ/HND or equivalent", educ_level == 6 ~ "Other professional qualifications",
    TRUE ~ as.character(NA)),
    educc_isced.lab = case_when(
      educ_level == 1 ~ "Level 5", educ_level == 2 ~ "Level 3", educ_level == 3 ~ "Level 2", 
      educ_level == 4 ~ "Level 2", educ_level == 5 ~ "Level 5", educ_level == 6 ~ "Level 4",
      educ_level == -7 ~ "Level 1", educ_level == -3 ~ is.na(educ_level) == as.character(NA)))



## Biochemical parameters ---------------------------------

biomark_fields <- c(
  chol = 30690, tg = 30870, ldl = 30780, hdl = 30760, glu = 30740, hba1c = 30750, crp = 30710
) ; biomark_vars <- setNames(paste0("f.", biomark_fields, ".0.0"), names(biomark_fields))

biomark_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb28679.tab.gz",
                    data.table=TRUE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, 
         all_of(biomark_vars))


## T2D (Eastwood algorithm + HbA1c < 5.7) ---------------------------------

t2d <- fread("../data/raw/UKB_Diabetes.csv", data.table=FALSE, stringsAsFactors=FALSE) %>% 
  rename(id = f.eid) %>%
  select("id", starts_with(c("probable", "possible")), "agedm_ts_or_ni_all", 
         "meds_any_sr_ni_all", "meds_any_sr_ni_ts_all", "dm_unlikely_all") %>% 
  left_join(fread("../data/raw/UKB_HbA1c.csv") %>% rename(id=f.eid) %>%
              select("id", "hba1c.30750.NGSP.max"),
            by = "id") %>% 
  left_join(fread("../data/raw/UKB_ICD10DM.csv") %>% rename(id=f.eid), 
            by = "id")

t2d_id <- t2d %>% mutate(
  t2d_med_any = case_when(
    meds_any_sr_ni_all == 1 ~ 1, 
    meds_any_sr_ni_ts_all == 1 ~ 1,
    meds_any_sr_ni_all == 0 & meds_any_sr_ni_ts_all == 0 ~ 0
  )) %>% 
  mutate(
    t2d_case = case_when(
      possible_t2dm_all == 1 ~ 1, 
      probable_t2dm_all == 1 ~ 1,
      hba1c.30750.NGSP.max >= 6.5 & possible_t1dm_all != 1 & probable_t1dm_all != 1 & dm_main != 1 & dm_secondary != 1 ~ 1,
      hba1c.30750.NGSP.max < 5.7 & dm_unlikely_all == 1 & dm_main != 1 & dm_secondary != 1 & t2d_med_any != 1 ~ 0,
      TRUE ~ as.numeric(NA)
    )) %>%
  mutate(
    t2d_case.f = case_when(
      as.numeric(t2d_case) == 0 ~ "Control",
      as.numeric(t2d_case) == 1 ~ "Case"),
    t2d_age_diagnosis = agedm_ts_or_ni_all
  ) %>%
  select("id", "t2d_case", "t2d_case.f", "hba1c_max"="hba1c.30750.NGSP.max", 
         "t2d_med_any", "t2d_age_diagnosis") 


print(paste0("DONE: Biochemical & anthropometric phenotypes prepared for ", nrow(t2d_id), " participants: ",
             "T2D Cases = ",  paste0(table(t2d_id$t2d_case.f))[1], ";",
             "T2D Controls = ", paste0(table(t2d_id$t2d_case.f))[2]) )


## Incident T2D -----------------------------------------

# grab incident t2d data
t2d_incid <- fread("/humgen/florezlab/UKBB_app27892/longitudinal_phenotypes/ t2d_longitudinal_phenotype_07022024_SH.csv") %>%
  select(id=eid, baseline_date, t2d, t2d_date, followup_date, censor_date, death_date) %>%
  mutate(across(ends_with("date"), ~as.Date(.))) %>%
  mutate(
    time_to_t2d = t2d_date - baseline_date, # 
    time_to_censor = censor_date - baseline_date) %>%
  mutate(time_to_event = ifelse(!is.na(time_to_t2d), time_to_t2d, time_to_censor))



##################################################################
##  Prep diet variables from FFQs & nutrient intakes from 24HR  ## 
##################################################################

print("Preparing dietary data ...")

### Functions to prepare 24HR data across multiple measurements (4/1 year) 
## Adapted with permission from KEW 

fetch_diet_fields <- function(fieldIDs, df, coding=FALSE) {
  # Given a list of fields constituting a food group:
  # - Determine the set of 24HR that are valid for that food group
  # - Recode the relevant variables based on their codings if necessary
  # - Sum over all fields for that food group within each instance
  # - Take the mean food group quantity over all instances
  diet_field_df <- lapply(0:4, function(i) {
    tcals_field <- paste0("f.26002.", i, ".0")  # Variable name for total calories in instance "i" 
    typical_diet_field <- paste0("f.100020.", i, ".0")  # Variable name for typical diet in instance "i" 
    valid_24hr <- (findInterval(df[[tcals_field]] / 4.18, c(600, 4800)) == 1) &
      df[[typical_diet_field]] == 1
    instance_fields <- paste0("f.", fieldIDs, ".", i, ".0")  # Variable names for all fields in instance "i"
    instance_df <- df[, instance_fields, drop=FALSE]
    if (coding) {  # Recode the variable if necessary (for food groups)
      instance_df <- mutate_all(instance_df, ~codings[as.character(.)])
    }
    ifelse(valid_24hr,  # Sum over fields if valid 24HR, else NA
           rowSums(instance_df, na.rm=TRUE), NA) } ) %>%
    setNames(paste0("instance", 0:4)) %>%
    bind_cols()
  diet_mean <- rowMeans(diet_field_df, na.rm=TRUE)
  ifelse(is.nan(diet_mean), NA, diet_mean)
}

check_num_valid_24hr <- function(df) {
  valid_24hr_df <- lapply(0:4, function(i) {
    tcals_field <- paste0("f.26002.", i, ".0")  # Variable name for total calories in instance "i" 
    typical_diet_field <- paste0("f.100020.", i, ".0")  # Variable name for typical diet in instance "i" 
    valid_24hr <- (findInterval(df[[tcals_field]] / 4.18, c(600, 4800)) == 1) &
      df[[typical_diet_field]] == 1
    valid_24hr
  }) %>%
    setNames(paste0("instance", 0:4)) %>%
    bind_cols()
  rowSums(valid_24hr_df, na.rm=TRUE)
}

#### food frequency questionnaires -------------------------------------

intake_fields <- c(cooked_veg = 1289, raw_veg = 1299,
                   fresh_fruit = 1309, dried_fruit = 1319,
                   bread_intake = 1438, bread_type = 1448, 
                   water = 1528, milk_type = 1418, spread_type = 1428,
                   spread_type_nonbutter=2654,
                   cereal_intake = 1458, cereal_type = 1468,
                   addsalt=1478, tea=1488, coffee = 1498, coffee_type = 1508,
                   hotdrink_temp = 1518) #**non_butter_spread_type

freq_fields <- c(oily_fish = 1329, nonoily_fish = 1339,
                 procmeat = 1349, poultry = 1359, cheese = 1408,
                 beef = 1369, lamb = 1379, pork = 1389)

ffq_fields<-c(intake_fields, freq_fields)
ffq_vars <- setNames(paste0("f.", ffq_fields, ".0.0"), names(ffq_fields))


## function to convert frequency values to servings/day (from KEW)
ffq_freq_to_sev <- function(x) {
  case_when(  # Data-coding 100377
    x == 5 ~ 1,  # "Once or more daily"
    x == 4 ~ 5.5 / 7,  # "5-6 times a week"
    x == 3 ~ 3 / 7,  # "2-4 times a week"
    x == 2 ~ 1 / 7,  # "Once a week"
    x == 1 ~ 0.5 / 7,  # "Less than once a week"
    x == 0 ~ 0,  # "Never"
    TRUE ~ as.numeric(NA)
  )
}


## function to recode negative values as meaninginful
neg_to_num <- function(x) {
  #x <- as.double(x) #"double" required to add values with decimals (previously, integer)
  case_when(
    x >= 0 ~ as.numeric(x), # -1 = "Do not know" ; -3 = "Prefer not to answer"
    x == -10 ~ 0.5, # -10 = "Less than 1 serving/day"
    TRUE ~ as.numeric(NA)
  )
}


## merge ffq variables & add total food groups
ffq_id <- main %>% select(id=f.eid, ffq_vars) %>%
  mutate(whole_bread = case_when(
    bread_type == 3 ~ 1,
    bread_type %in% c(1, 2, 4) ~ 0,
    TRUE ~ as.numeric(NA)
  )) %>%
  mutate(across(names(freq_fields), ffq_freq_to_sev)) %>%
  mutate(across(names(intake_fields), neg_to_num)) %>%
  mutate(total_veg = cooked_veg + raw_veg,
         total_fruit = fresh_fruit + dried_fruit,
         total_fish = oily_fish + nonoily_fish,
         red_meat = beef + lamb + pork,
         bread_intake = bread_intake / 7, # bread intake was provided in slices/week
         cereal_intake = cereal_intake / 7 # cereal intake was provided in bowls/week)
  ) %>%
  
  # Add FFQ vars for PCA analysis
  mutate(
    bread_type_white_vs_brown_or_whole = case_when(      
      bread_type == 1 ~ 1, bread_type == 2 | bread_type == 3 | 
        bread_type == 4 ~ 0, TRUE ~ as.numeric(NA)),
    milk_type_full_vs_low_or_nonfat = case_when(       
      milk_type == 1 ~ 1, milk_type == 2 | milk_type == 3 ~ 0, TRUE ~ as.numeric(NA)),
    milk_type_rare_never_BIN = case_when(
      milk_type == 6 ~ 1, milk_type != 6 ~ 0, TRUE ~ as.numeric(NA)),
    spread_type_butter_vs_any_other = case_when(          
      spread_type == 1 ~ 1, spread_type == 2 | spread_type == 3 ~ 0, TRUE ~ as.numeric(NA)),
    spread_type_rare_never_BIN = case_when(
      spread_type == 0 ~ 1, spread_type != 0 ~ 0, TRUE ~ as.numeric(NA)),
    cereal_type_sugar_vs_any_bran = case_when(
      cereal_type == 5 ~ 1, cereal_type != 5 ~ 0, TRUE ~ as.numeric(NA)),
    coffee_type_decaf_vs_regular = case_when(              
      coffee_type == 1 ~ 1, coffee_type == 2 | coffee_type == 3 | 
        coffee_type == 4 ~ 0, TRUE ~ as.numeric(NA)),
    addsalt_freq_QT = addsalt,        # INCLUDE
    addsalt_always_often_vs_nrs = case_when(
      addsalt == 3 | addsalt == 4 ~ 1, addsalt == 1 | addsalt == 2 ~ 0, TRUE ~ as.numeric(NA)),
    hotdrink_temp_hot_or_vhot_vs_warm = case_when(
      hotdrink_temp == 1  ~ 1, hotdrink_temp == 3 | hotdrink_temp == 2 ~ 0, TRUE ~ as.numeric(NA))
  )


## Recode categorical diet variables with descriptive levels

ffq_id <- ffq_id %>% 
  mutate(
    bread_type.lab = case_when(
      bread_type == 1 ~ "White", bread_type == 2 ~ "Brown", bread_type == 3 ~ "Wholemeal/Wholegrain",
      bread_type == 4 ~ "Other", bread_type == -1 ~ "Do not know", bread_type == -3 ~ "Prefer not to answer",
      TRUE ~ as.character(NA)),
    
    milk_type.lab = case_when(
      milk_type == 1 ~ "Full cream", milk_type == 2 ~ "Semi-skimmed", milk_type == 3 ~ "Skimmed",
      milk_type == 4 ~ "Soy", milk_type == 5 ~ "Other", milk_type == 6 ~ "Never/rarely have milk",
      milk_type == -1 ~ "Do not know", milk_type == -3 ~ "Prefer not to answer", 
      TRUE ~ as.character(NA)),
    
    spread_type.lab = case_when(
      spread_type == 1 ~ "Butter/spreadable butter", spread_type == 2 ~ "Flora Pro-Active/Benecol",
      spread_type == 3 ~ "Other spread/margarine", spread_type == 0 ~ "Never/rarely use spread",
      spread_type == -1 ~ "Do not know", spread_type == -3 ~ "Prefer not to answer"),
    
    spread_type_nonbutter.lab = case_when(
      spread_type_nonbutter == 4 ~ "Soft (tub) margarine", spread_type_nonbutter == 5 ~	"Hard (block) margarine",
      spread_type_nonbutter == 6 ~ "Olive oil based spread (eg: Bertolli)",
      spread_type_nonbutter == 7 ~ "Polyunsaturated/sunflower oil based spread (eg: Flora)",
      spread_type_nonbutter == 2 ~ "Flora Pro-Active or Benecol",
      spread_type_nonbutter == 8 ~ "Other low or reduced fat spread",
      spread_type_nonbutter == 9 ~ "Other type of spread/margarine", spread_type_nonbutter == -1 ~	"Do not know",
      spread_type_nonbutter == -3 ~ "Prefer not to answer"),
    
    cereal_type.lab = case_when(
      cereal_type == 1 ~ "Bran cereal (e.g. All Bran, Branflakes)", cereal_type == 2 ~ "Biscuit cereal (e.g. Weetabix)",
      cereal_type == 3 ~ "Oat cereal (e.g. Ready Brek, porridge)", cereal_type == 4 ~ "Muesli",
      cereal_type == 5 ~ "Other (e.g. Cornflakes, Frosties)", cereal_type == -1 ~ "Do not know",
      cereal_type == -3 ~ "Prefer not to answer"),
    
    coffee_type.lab = case_when(
      coffee_type == 1 ~ "Decaffeinated coffee (any type)", coffee_type == 2 ~ "Instant coffee",
      coffee_type == 3 ~ "Ground coffee (include espresso, filter etc)", coffee_type == 4 ~ "Other type of coffee",
      coffee_type == -1 ~ "Do not know", coffee_type == -3 ~ "Prefer not to answer"),
    
    addsalt.lab = case_when(
      addsalt == 1 ~ "Never/Rarely", addsalt == 2 ~ "Sometimes", addsalt == 3 ~ "Often",
      addsalt == 4 ~ "Always", TRUE ~ as.character(NA)),
    
    hotdrink_temp.lab = case_when(
      hotdrink_temp == 1 ~ "Very hot", hotdrink_temp == 2 ~ "Hot", hotdrink_temp == 3 ~ "Warm",
      TRUE ~ as.character(NA))
  )





#######################################
## Create phenotype dataset for GWAS ##
#######################################
cat("Writing .txt files for analysis .... \n")

allphenos <-dat %>% filter(id %in% geno_id$id) %>%
  left_join(geno_id, by = "id") %>%
  mutate(FID=id, IID=id, .before=id) %>% 
  left_join()
  distinct(IID, .keep_all=T)


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
  select(FID, IID, geno_array, ac, age, sex, bmi, paste0("PC", 1:20), smoke_former, smoke_current, alcohol, health, townsend, sleep) %>%
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
  select(FID, IID, geno_array, ac, age, sex, bmi, paste0("PC", 1:20), smoke_former, smoke_current, alcohol, health, townsend, sleep) %>%
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


