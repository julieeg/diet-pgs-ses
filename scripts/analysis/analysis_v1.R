## Preliminary pgs-ses analysis

library(tidyverse)
library(data.table)

lapply(list.files("../../pantry/functions/", full.names = T), source)


## Preference phenos/PGS vars & labels 
pheno_labs <- c(pref_sweet_food="Prefer sweet food", pref_bitter_food="Prefer bitter food", 
               pref_salty_food="Prefer salty food", pref_fatty_food="Prefer fatty food")
phenoPalettes <- list(
  pref_sweet_food=c(palettes$NatExt$Reds[c(2,5)]), pref_bitter_food=c(palettes$NatExt$Blues[c(2,5)]),
  pref_salty_food=c(palettes$NatExt$Oranges[c(2,5)]), pref_fatty_food=c(palettes$NatExt$Purples[c(2,5)]) )

pgs_labs <- paste0(as.vector(pheno_labs), " PGS")
names(pgs_labs) <- paste0(names(pheno_labs),"_pgs")
pgsPalettes <- phenoPalettes ; names(pgsPalettes)<-paste0(names(phenoPalettes),"_pgs")

pgs <- names(pgs_labs)

################################################################################
## Define variable groups
################################################################################

basic_vars <- c(
  carb_pheno="Carb intake (%kcal)", fat_pheno="Fat intake (%kcal)", prot_pheno="Protein (%kcal)", 
  ancestry="Genetic ancestry", sex.lab="Sex, female", age="Age, years", bmi= "BMI, kg/m2",  
  smoke_level.lab="Smoking status", alch_level.lab="Alcohol use", health.lab="Self-reported health", 
  sleep="Sleep, hrs/day", townsend="Townsend deprivation index", 
  educ_level.lab="Education level (UKB)", educ_isced.lab="Education level (ISCED)",
  educ_4level.lab="Education level (4-level)", income_level.lab="Income level")

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
               addsalt.lab="Frequency of adding salt to food")

diet24hr_vars <- c(fg_bevs_alcoholic = "Alcoholic beverages", fg_cereals="Cereals and cereal products",
                    df_dairy_products = "Dairy and dairy-free products", fg_eggs = "Egg and egg-dishes",
                    fg_fats_spreads = "Fat and spreads", fg_fish_dishes="Fish and fish dishes",
                    fg_fruits="Fruits", fg_meat_products="Meat and meat products", 
                    fg_meat_substitutes="Meat substitutes", fg_mixed_dishes="Mixed dishes",                                             
                    fg_bevs_nonalch="Non-alcoholic beverages", fg_nuts_seeds="Nuts and seeds",
                    fg_condiments="Sauces and condiments", fg_sweets_snacks="Sugar, preserves, cakes and confectionary, snacks",
                    fg_vegetables_potatoes="Vegetables and potatoes")


################################################################################
## Load in analysis data 
################################################################################

## Load testing data
analysis <- readRDS("../data/processed/ukb_MA_analysis_testing.rda")

## Full (training+testing) data
fullsample <- readRDS("../data/processed/ukb_MA_allphenos.rda")

training_ids <- fread("../data/processed/phenos/ukb_phenos_gwas_MA_training.txt") %>% pull(IID)
training <- fullsample %>% filter(id %in% training_ids)
testing <- fullsample %>% filter(!id %in% training_ids)

# descriptive table for training set ------
print_summary_table(data=training, vars_to_summarize = basic_vars, var_strata = F, p_adjust = "none") %>%
  rename(Training=Total) %>%
  fwrite("../data/results/descr_tab_basic_training.csv", row.names = T)

# descriptive table for testing set ------
print_summary_table(data=testing, vars_to_summarize = basic_vars, var_strata = F, p_adjust = "none") %>%
  rename(Testing=Total) %>%
  fwrite("../data/results/descr_tab_basic_testing.csv", row.names = T)

## TEMPORARY: analysis
print_summary_table(data=analysis, vars_to_summarize = basic_vars, var_strata = F, p_adjust = "none") %>%
  rename(Testing=Total) %>%
  fwrite("../data/results/descr_tab_basic_testing_retire_20250413.csv", row.names = T)



#########################################################################
## Summarize participant characteristics by preference PGS quintile ##
#########################################################################

tab_descr_pgs.l <- lapply(names(pgs_labs), function(pgs) {
  
  # Quintiles
  #strata <- paste0(pgs, "_quintile")
  #strata_order <- 1:5
  
  # Deciles
  strata <- paste0(pgs, "_decile")
  strata_order <- 1:10

  ## Basic descriptive -----------------------------
  descr <- print_summary_table(data=analysis, vars_to_summarize = basic_vars, 
                      var_strata = strata, var_strata_order = strata_order,
                      digits = c(1,1,3), p_adjust = c("agesex")) 
  
  ## Cardiometabolic risk factors -------------------------
  cmrf <- print_summary_table(data=analysis, vars_to_summarize = cmrf_vars, 
                              var_strata = strata, var_strata_order = strata_order,
                              digits = c(1,1,3), p_adjust = "agesex") 
  
  ## FFQ & food preferences -------------------------
  diet <- print_summary_table(data=analysis, vars_to_summarize = diet_vars, 
                              var_strata = strata, var_strata_order = strata_order,
                              p_adjust = "agesex", digits = c(1,1,3)) 
  
  ## 24hr dietary recall ------------------------------
  diet24hr <- print_summary_table(data=analysis, vars_to_summarize = diet24hr_vars, 
                              var_strata = strata, var_strata_order = strata_order,
                              p_adjust = "agesex", digits = c(1,1,3)) 
  
  
  return(list(tab_descr=descr, tab_cmrf=cmrf, tab_diet=diet, tab_diet24hr=diet24hr))
  
  }) ; names(tab_descr_pgs.l)<-names(pgs_labs)

#tab_descr_pgs.l %>% saveRDS("../data/results/descr_tablist_pgs_allPhenos_byQuintile_testing.rda")
#tab_descr_pgs.l %>% saveRDS("../data/results/descr_tablist_pgs_allPhenos_byDecile_testing.rda")

# ===========================================
## Summary tables by quintiles / deciles 
# ===========================================

## Combine tables for Q1/Q3/Q5 for each set of descriptive characteristics
tab_descr_pgs.l <- readRDS("../data/results/descr_tablist_pgs_allPhenos_byQuintile_testing.rda")
lapply(names(tab_descr_pgs.l[[1]]), function(tab) {
  
  do.call(cbind.data.frame, lapply(pgs, function(p) {
    tab_descr_pgs.l[[p]][[tab]] %>% rename_at(2:6, ~paste0(gsub("_food","", p), "_Q", .)) %>% 
      rename_at(c("P_value", "P_trend"), ~paste0(gsub("_food","", p) ,"_",.)) %>%
      filter(!startsWith(rownames(.), " NA")) %>%
      select("Total", contains(paste0("_Q", c(1,3,5))), contains("_P_"))
  })) %>% 
    
    fwrite(paste0("../data/results/descr_tab_allPGS_",gsub("tab_","",tab),"_byQuintile_testing.csv"), row.names = T)
})


## Combine tables for Q1/Q5/Q10 for each set of descriptive characteristics
tab_descr_pgs.l <- readRDS("../data/results/descr_tablist_pgs_allPhenos_byDecile_testing.rda")
lapply(names(tab_descr_pgs.l[[1]]), function(tab) {
  
  do.call(cbind.data.frame, lapply(pgs, function(p) {
    tab_descr_pgs.l[[p]][[tab]] %>% rename_at(2:11, ~paste0(gsub("_food","", p), "_D", .)) %>% 
      rename_at(c("P_value", "P_trend"), ~paste0(gsub("_food","", p) ,"_",.)) %>%
      filter(!startsWith(rownames(.), " NA")) %>%
      select("Total", contains(paste0("_Q", c(1,5,10))), contains("_P_"))
  })) %>% 
    
    fwrite(paste0("../data/results/descr_tab_allPGS_",gsub("tab_","",tab),"_byDecile_testing.csv"), row.names = T)
})
  


######################################################
## Correlations of macro PGS with dietary variables ##
######################################################

## Pearson correlations with all dietary vars -------------------------------
cor_pgs_diet <- cor(analysis %>% select(
  pgs, names(diet_vars), names(diet24hr_vars)) %>%
    select(-ends_with(".lab")) %>% 
    mutate_all(., ~as.numeric(.)) %>%
    filter(complete.cases(.))) %>%
  as.data.frame() %>%
  filter(rownames(.) %in% pgs) %>%
  as.matrix()

cor_pgs_diet %>% fwrite("../data/results/descr_cormat_pgs_allDiet_testing.csv")


## Sex, Age, gPC adjusted associations with PREF traits -------------------------------
codebook_prefs <- readxl::read_xlsx("../../tasteprefs/data/raw/food_pref_vars.xlsx")
codebook_24hr <- readxl::read_xlsx("../../rediMR/run/ukb_24hr_codebook.xlsx")
tastes <- c("bitter", "sweet", "salt", "sour", "umami", "fatty")

m_base <- paste0(c("age", "sex", paste0("PC",1:10), "ancestry"),collapse = "+")
pref_vars <- analysis %>% select(-c("id", "survey_date", "survey_duration", "p20599"))
pref_vars <- names(pref_vars)[names(pref_vars) %in%  codebook_prefs$varID]


#############################################################
## Generalized linear models of assoc. with dietary intake ##
#############################################################

m_base <- paste0(c("age", "sex", paste0("PC",1:10), "ancestry"),collapse = "+")
m_bmi <- paste0(m_base,"bmi", collapse = "+")
m_life <- paste0(c(m_bmi,"smoke_current", "smoke_former", "alcohol.lab", "energy_pheno", 
                   "health.lab"), collapse = "+")

## mbase: lm of macro PGS with food preference traits 
lm_pgs_foodprefs <- do.call(rbind.data.frame, lapply(1:length(pgs_labs), function(p) {
  do.call(rbind.data.frame, lapply(pref_vars, function(var) {
    print_lm(exposure=names(pgs_labs)[p], outcome=var, covariates = m_base, data=analysis, 
             label=codebook_prefs$description[which(codebook_prefs$varID == var)]) 
    }))
})) ; lm_pgs_foodprefs

lm_pgs_foodprefs %>% fwrite("../data/results/lm_prefpgs_foodprefs_mBase_testing.csv")

##NOTE: RE-DID THIS WITH PREFERENCS, WHICH MAY OR MAY NOT BE HELPFUL ...... 

## Add adjustment for macronutrint intake (KEW suggestions)
lm_pgs_foodprefs_adjMacro <- do.call(rbind.data.frame, lapply(macro_pgs, function(pgs) {
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









