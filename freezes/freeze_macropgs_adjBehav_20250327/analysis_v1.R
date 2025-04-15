## Preliminary pgs-ses analysis

library(tidyverse)
library(data.table)

lapply(list.files("../../pantry/functions/", full.names = T), source)

macros <- c("carb", "fat", "prot")
macro_phenos<-paste0(macros, "_pheno")
macro_pgs<-paste0(macros, "_pgs")

basic_vars <- c(
  carb_pheno="Carb intake (%kcal)", fat_pheno="Fat intake (%kcal)", prot_pheno="Protein (%kcal)", 
  ancestry="Genetic ancestry", sex.lab="Sex, female", age="Age, years", bmi= "BMI, kg/m2",  
  smoking.lab="Smoking status", alcohol.lab="Alcohol use", health.lab="Self-reported health", 
  sleep="Sleep, hrs/day", townsend="Townsend deprivation index", 
  educ_level.lab="Education level (UKB)", educ_isced.lab="Education level (ISCED)",
  educ_4lvl.lab="Education level (4-level)", income_level="Income level")


cmrf_vars <- c(
  sbp="SBP, mmHg", dbp="DBP, mmHg", waist="Waist circumference, cm", 
  chol="Total cholesterol, mmol/L", tg_log="log(TG)", ldl="LDL, mmol/L",
  hdl="HDL, mmol/L", glu="Glucose, mmol/L", hba1c="HbA1c, %", crp="CRP, mg/dL",
  t2d="Incident T2D")


diet_vars <- c(cooked_veg="Cooked veg", raw_veg="Raw veg", fresh_fruit="Fresh fruit", 
               dried_fruit="Dried fruit", oily_fish="Oily fish", nonoily_fish="Non-oily fish", 
               bread_intake="Bread intake", cheese="Cheese intake", cereal_intake="Cereal intake",
               poultry="Poultry intake", procmeat="Processed meat intake", red_meat="Red meat intake",
               beef="Beef intake", lamb="Lamb intake", pork="Pork intake", 
               tea="Tea intake", coffee="Coffee intake", water="Water intake",
               bread_type.lab="Prefer bread type", milk_type.lab="Prefer milk type", 
               cereal_type.lab="Prefer cereal type", spread_type.lab="Prefer spread type", 
               coffee_type.lab="Prefer coffee type", hotdrink_temp.lab="Prefer drink temperature", 
               addsalt.lab="Frequency of adding salt to food", pref_sweet_food="Prefer sweet food (1-9)", 
               pref_bitter_food="Prefer bitter foods (1-9)", pref_fatty_food="Prefer fatty food", 
               pref_salty_food="Prefer salty food")

diet24hr_vars <- c(fg_bevs_alcoholic = "Alcoholic beverages", fg_cereals="Cereals and cereal products",
                    df_dairy_products = "Dairy and dairy-free products", fg_eggs = "Egg and egg-dishes",
                    fg_fats_spreads = "Fat and spreads", fg_fish_dishes="Fish and fish dishes",
                    fg_fruits="Fruits", fg_meat_products="Meat and meat products", 
                    fg_meat_substitutes="Meat substitutes", fg_mixed_dishes="Mixed dishes",                                             
                    fg_bevs_nonalch="Non-alcoholic beverages", fg_nuts_seeds="Nuts and seeds",
                    fg_condiments="Sauces and condiments", fg_sweets_snacks="Sugar, preserves, cakes and confectionary, snacks",
                    fg_vegetables_potatoes="Vegetables and potatoes")



################################################################################
## Prepare analysis dataframe
################################################################################

# ===================================
## Prepare FULL data 
# ===================================

analysis_full <- readRDS("../data/processed/phenos/ukb_analysis_MA.rda") %>%
  mutate(educ_4lvl.lab = case_when(
    educ_level %in% c(1,5,6) ~ "College graduate",
    educ_level == 2 ~ "Some college",
    educ_level %in% c(3,4) ~ "Some HS",
    educ_level %in% c(-7) ~ "Less than HS",
    TRUE ~ as.character(NA))
  ) %>%
  mutate() %>%
  mutate_at("ancestry", ~ifelse(.=="","Missing",.))

# 24hr data  
diet24hr <- readRDS("../data/processed/phenos/ukb_diet24hr.rda")

## pgs_full
pgs_full <- analysis_full %>% select(id, contains("pgs")) #N=15967

# merge
analysis_full <- analysis_full %>% left_join(diet24hr %>% select(id, starts_with("fg"), "df_dairy_products"), by="id")


# =====================================
## Prepare TESTING data
# =====================================

# restrict to testing samples
training_samples <- fread("../data/processed/gwas/training/ukb_phenos_gwas_macros_MA_samples_training.txt")
testing_samples <- analysis_full %>% filter(!id %in% training_samples$IID) %>% pull(id)

## Load testing PGS
pgs_testing <- left_join(fread("../data/processed/pgs/ukb_allchr_Y1_MA_testing_pgs") %>% 
                      rename(carb_pgs=score_sum), 
                    fread("../data/processed/pgs/ukb_allchr_Y2_MA_testing_pgs") %>% 
                      rename(fat_pgs=score_sum), by="id") %>%
  left_join(fread("../data/processed/pgs/ukb_allchr_Y3_MA_testing_pgs") %>% 
              rename(prot_pgs=score_sum), by="id")

# Add pgs to analysis & restrict to testing samples 
analysis <- inner_join(analysis_full %>% select(-contains("pgs")), pgs_testing, by = "id") 
dim(analysis) #N=45848


## Reconstruct quintiles & deciles for macro_pgs
analysis <- analysis %>% 
  mutate(carb_pgs_decile = cut(carb_pgs, breaks = quantile(carb_pgs, probs = seq(0, 1, 0.1), na.rm = TRUE), 
                               include.lowest = TRUE, labels = 1:10),
         fat_pgs_decile = cut(fat_pgs, breaks = quantile(fat_pgs, probs = seq(0, 1, 0.1), na.rm = TRUE), 
                              include.lowest = TRUE, labels = 1:10),
         prot_pgs_decile = cut(prot_pgs, breaks = quantile(prot_pgs, probs = seq(0, 1, 0.1), na.rm = TRUE), 
                               include.lowest = TRUE, labels = 1:10) ) %>%
  mutate(carb_pgs_quintile = cut(carb_pgs, breaks = quantile(carb_pgs, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                                 include.lowest = TRUE, labels = 1:5),
         fat_pgs_quintile = cut(fat_pgs, breaks = quantile(fat_pgs, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                                include.lowest = TRUE, labels = 1:5),
         prot_pgs_quintile = cut(prot_pgs, breaks = quantile(prot_pgs, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                                 include.lowest = TRUE, labels = 1:5)) %>%
  mutate(across(contains("decile"), ~as.factor(.))) %>%
  mutate(across(contains("quintile"), ~as.factor(.)))



# =====================================
## Prepare TRAINING data
# =====================================

# Add pgs to analysis & restrict to testing samples 
analysis_training <- analysis_full %>% select(-contains("pgs")) %>% filter(id %in% training_samples$IID)
dim(analysis_training) #N=107011

analysis_training %>% fwrite("../data/processed/phenos/ukb_analysis_MA_training.csv")

print_summary_table(data=analysis_training, vars_to_summarize = basic_vars, 
                    var_strata = F, p_adjust = "none") %>%
  rename(Training=Total) %>%
  fwrite("../data/processed/results/tab_descr_basic_training.csv", row.names = T)


#########################################################################
## Summarize participant characteristics by macronutrient PGS quintile ##
#########################################################################

tab_descr_macros.l <- lapply(c("carb", "fat", "prot"), function(macro) {
  
  # Quintiles
  #strata <- paste0(macro, "_pgs_quintile")
  #strata_order <- 1:5
  
  # Deciles
  strata <- paste0(macro, "_pgs_decile")
  strata_order <- 1:10

  ## Basic descriptive -----------------------------
  descr <- print_summary_table(data=analysis, vars_to_summarize = basic_vars, 
                      var_strata = strata, var_strata_order = strata_order,
                      digits = c(1,1), p_adjust = "none") 
  
  ## Cardiometabolic risk factors -------------------------
  cmrf <- print_summary_table(data=analysis, vars_to_summarize = cmrf_vars, 
                              var_strata = strata, var_strata_order = strata_order,
                              digits = c(1,1), p_adjust = "none") 
  
  ## FFQ & food preferences -------------------------
  diet <- print_summary_table(data=analysis, vars_to_summarize = diet_vars, 
                              var_strata = strata, var_strata_order = strata_order,
                              digits = c(1,1), p_adjust = "none") 
  
  ## 24hr dietary recall ------------------------------
  diet24hr <- print_summary_table(data=analysis, vars_to_summarize = diet24hr_vars, 
                              var_strata = strata, var_strata_order = strata_order,
                              digits = c(1,1), p_adjust = "none") 
  
  
  return(list(tab_descr=descr, tab_cmrf=cmrf, tab_diet=diet, tab_diet24hr=diet24hr))
  
  }) ; names(tab_descr_macros.l)<-macros

#tab_descr_macros.l %>% saveRDS("../data/processed/tab_descr_macros_allPhenos_byDecile_testing.rda")

tab_descr_macros.l <- readRDS("../data/processed/tab_descr_macros_allPhenos_testing.rda")
names(tab_descr_macros.l)<-macros

################################################################################
############################## by PGS QUINTILE #################################

## Combine tables for Q1/Q3/Q5 for each set of descriptive characteristics

# Basic descriptives -------------
tab_descr_allMacros <- do.call(cbind.data.frame, lapply(macros, function(macro) {
  tab_descr_macros.l[[macro]]$tab_descr %>% rename_at(2:6, ~paste0(macro, "_Q", .)) %>% 
    rename_at(c("P_value", "P_trend"), ~paste0(macro,"_",.)) %>%
    filter(!startsWith(rownames(.), " NA")) %>%
    select("Total", paste0(macro,"_Q", c(1,3,5)), contains("_P_"))
})) %>% fwrite("../output/tab_descr_basic_allMacros_testing.csv", row.names = T)


# Basic cmrf -------------
tab_cmrf_allMacros <- do.call(cbind.data.frame, lapply(macros, function(macro) {
  tab_descr_macros.l[[macro]]$tab_cmrf %>% rename_at(2:6, ~paste0(macro, "_Q", .)) %>% 
    rename_at(c("P_value", "P_trend"), ~paste0(macro,"_",.)) %>%
    filter(!startsWith(rownames(.), " NA")) %>%
    select("Total", paste0(macro,"_Q", c(1,3,5)), contains("_P_"))
})) %>% fwrite("../output/tab_descr_cmrf_allMacros_testing.csv", row.names = T)


# Basic diet -------------
tab_diet_allMacros <- do.call(cbind.data.frame, lapply(macros, function(macro) {
  tab_descr_macros.l[[macro]]$tab_diet %>% rename_at(2:6, ~paste0(macro, "_Q", .)) %>% 
    rename_at(c("P_value", "P_trend"), ~paste0(macro,"_",.)) %>%
    filter(!startsWith(rownames(.), " NA")) %>%
    select("Total", paste0(macro,"_Q", c(1,3,5)), contains("_P_"))
})) %>% fwrite("../output/tab_descr_diet_allMacros_testing.csv", row.names = T)


# 24hr dietary diet -------------
tab_diet_allMacros <- do.call(cbind.data.frame, lapply(macros, function(macro) {
  tab_descr_macros.l[[macro]]$tab_diet %>% rename_at(2:6, ~paste0(macro, "_Q", .)) %>% 
    rename_at(c("P_value", "P_trend"), ~paste0(macro,"_",.)) %>%
    filter(!startsWith(rownames(.), " NA")) %>%
    select("Total", paste0(macro,"_Q", c(1,3,5)), contains("_P_"))
})) %>% fwrite("../output/tab_descr_diet_allMacros_testing.csv", row.names = T)

################################################################################
################################################################################


################################################################################
############################### by PGS DECILE ##################################

## Combine tables for Q1/Q5/Q10 for each set of descriptive characteristics

# Basic descriptives -------------
do.call(cbind.data.frame, lapply(macros, function(macro) {
  tab_descr_macros.l[[macro]]$tab_descr %>% rename_at(2:11, ~paste0(macro, "_D", .)) %>% 
    rename_at(c("P_value", "P_trend"), ~paste0(macro,"_",.)) %>%
    filter(!startsWith(rownames(.), " NA")) %>%
    select("Total", paste0(macro,"_D", c(1,5,10)), contains("_P_"))
})) %>% fwrite("../output/tab_descr_basic_allMacros_byDecile_testing.csv", row.names = T)


# Basic cmrf -------------
do.call(cbind.data.frame, lapply(macros, function(macro) {
  tab_descr_macros.l[[macro]]$tab_cmrf %>% rename_at(2:11, ~paste0(macro, "_D", .)) %>% 
    rename_at(c("P_value", "P_trend"), ~paste0(macro,"_",.)) %>%
    filter(!startsWith(rownames(.), " NA")) %>%
    select("Total", paste0(macro,"_D", c(1,5,10)), contains("_P_"))
})) %>% fwrite("../output/tab_descr_cmrf_allMacros_byDecile_testing.csv", row.names = T)


# Basic diet -------------
do.call(cbind.data.frame, lapply(macros, function(macro) {
  tab_descr_macros.l[[macro]]$tab_diet %>% rename_at(2:11, ~paste0(macro, "_D", .)) %>% 
    rename_at(c("P_value", "P_trend"), ~paste0(macro,"_",.)) %>%
    filter(!startsWith(rownames(.), " NA")) %>%
    select("Total", paste0(macro,"_D", c(1,5,10)), contains("_P_"))
})) %>% fwrite("../output/tab_descr_diet_allMacros_byDecile_testing.csv", row.names = T)


# 24hr dietary diet -------------
do.call(cbind.data.frame, lapply(macros, function(macro) {
  tab_descr_macros.l[[macro]]$tab_diet %>% rename_at(2:11, ~paste0(macro, "_D", .)) %>% 
    rename_at(c("P_value", "P_trend"), ~paste0(macro,"_",.)) %>%
    filter(!startsWith(rownames(.), " NA")) %>%
    select("Total", paste0(macro,"_D", c(1,5,10)), contains("_P_"))
})) %>% fwrite("../output/tab_descr_diet_allMacros_byDecile_testing.csv", row.names = T)

################################################################################
################################################################################


######################################################
## Correlations of macro PGS with dietary variables ##
######################################################

## Pearson correlations with all dietary vars -------------------------------
cor_pgs_diet <- cor(analysis %>% select(
  macro_pgs, energy_pheno, names(diet_vars), names(diet24hr_vars)) %>%
    select(-ends_with(".lab")) %>% 
    mutate_all(., ~as.numeric(.)) %>%
    filter(complete.cases(.))) %>%
  as.data.frame() %>%
  filter(rownames(.) %in% macro_pgs) %>%
  as.matrix()

cor_pgs_diet %>% fwrite("../data/processed/results/tab_cormat_pgs_diet_testing.csv")

## Sex, Age, gPC adjusted associations with PREF traits -------------------------------
prefs_id <- fread("../../tasteprefs/data/raw//ukbb_app27892_diet_preference_11062024.csv") %>%
  rename(id=eid, survey_date=p20750, survey_duration=p20751) %>%
  filter(!is.na(survey_date)) %>% # filter out participants without preference data (N=320071)
  mutate(across(starts_with("p"), ~as.numeric(
    case_when(
      . %in% c(1:10) ~ as.character(.),
      . == "extremely dislike" ~ "1",
      . == "neither like nor dislike" ~ "5",
      . == "extremely like" ~ "9",
      . %in% c("never tried", "do not wish to answer") ~ "NA",
      TRUE ~ as.character(NA) ) ))
  ) 

analysis_prefs <- left_join(analysis, prefs_id)

codebook_prefs <- readxl::read_xlsx("../../tasteprefs/data/raw/food_pref_vars.xlsx")
tastes <- c("bitter", "sweet", "salt", "sour", "umami", "fatty")

m_base <- paste0(c("age", "sex", paste0("PC",1:10), "ancestry"),collapse = "+")
pref_vars <- names(prefs_id %>% select(-c("id", "survey_date", "survey_duration", "p20599")))

  

#############################################################
## Generalized linear models of assoc. with dietary intake ##
#############################################################

m_base <- paste0(c("age", "sex", paste0("PC",1:10), "ancestry"),collapse = "+")
m_bmi <- paste0(m_base,"bmi", collapse = "+")
m_life <- paste0(c(m_bmi,"smoke_current", "smoke_former", "alcohol.lab", "energy_pheno", 
                   "health.lab"), collapse = "+")


## mbase: lm of macro PGS with food preference traits 
lm_macroPGS_allPrefs <- do.call(rbind.data.frame, lapply(macro_pgs, function(pgs) {
  do.call(rbind.data.frame, lapply(pref_vars, function(var) {
    print_lm(exposure=pgs, outcome=var, covariates = m_base, data=analysis_prefs, 
             label=codebook_prefs$description[which(codebook_prefs$varID == var)]) 
    })) %>% mutate(macronutrient_pgs = pgs)
})) ; lm_macroPGS_allPrefs

lm_macroPGS_allPrefs %>% fwrite("../data/processed/results/lm_macroPGS_allPrefs_mBase_testing.csv")


## Add adjustment for macronutrint intake (KEW suggestions)
lm_macroPGS_allPrefs_adjMacro <- do.call(rbind.data.frame, lapply(macro_pgs, function(pgs) {
  do.call(rbind.data.frame, lapply(pref_vars, function(var) {
    cov <- paste0(m_base, "+", ifelse(pgs=="carb_pgs", "carb_pheno", 
                                     ifelse(pgs=="fat_pgs", "fat_pheno", "prot_pheno"))) 
    print_lm(exposure=pgs, outcome=var, covariates = cov, data=analysis_prefs, 
             label=codebook_prefs$description[which(codebook_prefs$varID == var)]) 
  })) %>% mutate(macronutrient_pgs = pgs) %>% mutate(adjusted="T")
})) ; lm_macroPGS_allPrefs_adjMacro

lm_macroPGS_allPrefs_adjMacro %>% fwrite("../data/processed/results/lm_macroPGS_allPrefs_mBase_adjMacro_testing.csv")



########################
## Save analysis data ##
########################

analysis %>% saveRDS("../data/processed/phenos/ukb_analysis_testing.rda")


