
## load packages & pantry scripts
library(tidyverse) 
library(data.table)

lapply(list.files("../../opt/pantry/functions/", full.names = T), source)


#######################
##  Basic functions  ##
#######################

## Add descriptive labels 
descr_label.fun <- function(data, base_var, labs_vals) {
  base <- data %>% select(all_of(base_var)) 
  temp <- rep(NA, length(base))
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; return(temp)
}

descr_label.fun <- function(data, base_var, labs_vals, ordered = F) {
  base <- data %>% select(all_of(base_var)) 
  x <- rep(NA, length(base))
  for(i in 1:length(labs_vals)) {
    x[base == labs_vals[i] ] <- names(labs_vals)[i]
  } ; if(ordered == T) {
    x <- factor(x, levels=unique(names(labs_vals)))
  } ; return(x)
}



## winsorize
winsorize <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x<bounds[1], bounds[1], ifelse(x>bounds[2], bounds[2], x))
  x
}


################################################
## Load parent phenotypes from UKB ## 
################################################

## Base phenotypes -----------------------------

base <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz", 
              data.table=FALSE, stringsAsFactors=FALSE)

base_id <- base %>% 
  select(id = f.eid,
         ac = f.54.0.0,
         ac_date = f.53.0.0,
         sbp = f.4080.0.0,
         dbp = f.4079.0.0,
         waist = f.48.0.0,
         fasting_hrs = f.74.0.0,
         med_code = f.20003.0.0,
         med_mets = f.6177.0.0)
  
withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/withdraw/withdraw27892_459_20240527.txt", what=character())


## Education level ---------------------------------------------------

## Coding based on: Ge T., et al. Cerebral Cortex 2019;29(8): 3471-3481.
educ_level_labs <- list("None of the above" = -7, "Prefer not to answer"= -3, 
  "College or university degree" = 1, "A/AS levels or equivalent" = 2, 
  "O/GCSE levels or equivalent" = 3, "CSEs or equivalent" = 4,
  "NVQ/HND or equivalent" = 5, "Other professional qualifications" = 6)

educ_isced_level_labs <- list("Level 5" = 1, "Level 5" = 5,  "Level 4" = 6, 
                              "Level 3" = 2, "Level 2" = 3, "Level 2" = 4, 
                              "Level 1" = -7)  # NA = -3 or missing

educ_4level_labs = list( "College graduate" = 1, "College gradute" = 5, "Some college" = 6,
                      "Some college" = 2, "Some HS" = 3, "Some HS" = 4, "Less than HS" = -7)


educ_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672670.tab.gz",
                 data.table = FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, educ_level = f.6138.0.0) %>%
  mutate(educ_level.lab = descr_label.fun(., "educ_level", educ_level_labs),
         educ_isced.lab = descr_label.fun(., "educ_level", educ_isced_level_labs, ordered = T),
         educ_4level.lab = descr_label.fun(., "educ_level", educ_4level_labs, ordered = T))



# Income level ---------------------------------------------------

income_labs <- list("Less than 18,000"=1, "18,000 to 30,999"=2, "31,000 to 51,999"=3, 
  "52,000 to 100,000"=4, "Greater than 100,000"=5, "Do not know"=-1,
  "Prefer not to answer"=-3)

income_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672750.tab.gz",
                   data.table = FALSE, stringsAsFactors = FALSE) %>%
  select(id=f.eid, income=f.738.0.0) %>%
  mutate(income_level.lab = descr_label.fun(., "income", income_labs, ordered=T))



## Biochemical parameters ---------------------------------

biomark_fields <- c(
  chol = 30690, tg = 30870, ldl = 30780, hdl = 30760, glu = 30740, hba1c = 30750, crp = 30710
) ; biomark_vars <- setNames(paste0("f.", biomark_fields, ".0.0"), names(biomark_fields))

biochem_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb28679.tab.gz",
                    data.table=TRUE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, 
         all_of(biomark_vars)) %>%
  mutate(tg_log=log(tg))



## T2D (Eastwood algorithm + HbA1c < 5.7) ---------------------------------

t2d <- fread("../../taste2d/data/raw/UKB_Diabetes.csv", data.table=FALSE, stringsAsFactors=FALSE) %>% 
  rename(id = f.eid) %>%
  select("id", starts_with(c("probable", "possible")), "agedm_ts_or_ni_all", 
         "meds_any_sr_ni_all", "meds_any_sr_ni_ts_all", "dm_unlikely_all") %>% 
  left_join(fread("../../taste2d/data/raw/UKB_HbA1c.csv") %>% rename(id=f.eid) %>%
              select("id", "hba1c.30750.NGSP.max"),
            by = "id") %>% 
  left_join(fread("../../taste2d/data/raw/UKB_ICD10DM.csv") %>% rename(id=f.eid), 
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
t2d_incid <- fread("/humgen/florezlab/UKBB_app27892/incident_phenotypes/t2d_incident_phenotype_07022024_SH.csv") %>%
  select(id=eid, baseline_date, t2d, t2d_date, followup_date, censor_date, death_date) %>%
  mutate(across(ends_with("date"), ~as.Date(.))) %>%
  mutate(
    time_to_t2d = t2d_date - baseline_date, # 
    time_to_censor = censor_date - baseline_date) %>%
  mutate(time_to_event = ifelse(!is.na(time_to_t2d), time_to_t2d, time_to_censor))



## food frequency questionnaires -------------------------------------

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
ffq_id <- base %>% select(id=f.eid, ffq_vars) %>%
  mutate(whole_bread = case_when(
    bread_type == 3 ~ 1,
    bread_type %in% c(1, 2, 4) ~ 0,
    TRUE ~ as.numeric(NA) )) %>%
  mutate(across(names(freq_fields), ffq_freq_to_sev)) %>%
  mutate(across(names(intake_fields), neg_to_num)) %>%
  mutate(total_veg = cooked_veg + raw_veg,
         total_fruit = fresh_fruit + dried_fruit,
         total_fish = oily_fish + nonoily_fish,
         red_meat = beef + lamb + pork,
         bread_intake = bread_intake / 7, # bread intake was provided in slices/week
         cereal_intake = cereal_intake / 7 # cereal intake was provided in bowls/week)
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


## Load in 24hr data (prepared with: prep_ukb_phenos_24hr.R)

diet24hr_id <-readRDS("../data/processed/ukb_diet24hr.rda")


##  Food preference traits --------------------------

#Food preference traits processing
#-Food preference traits from the UKB were assessed on a likert scale from 0 (extremely dislike) to 10 (extremely like); 5 was considered “neither like nor dislike”. 
#-Values of “never tried” and “do not wish to answer” were recoded as missing. 
#-Variables were then recoded as numeric values, for use in PCA

prefs_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_jul_2024/ukbb_app27892_diet_preference_11062024.csv") %>%
  rename(id=eid, survey_date=p20750, survey_duration=p20751) %>%
  filter(!is.na(survey_date)) %>% # filter out participants without preference data (N=320071)
  mutate(across(starts_with("p"), ~as.numeric(
    case_when(
      . %in% c(1:10) ~ as.character(.),
      . == "extremely dislike" ~ "1",
      . == "neither like nor dislike" ~ "5",
      . == "extremely like" ~ "9",
      . %in% c("never tried", "do not wish to answer") ~ "NA",
      TRUE ~ as.character(NA) ) ))) #%>%
  #select(id, pref_sweet_food=p20732, 
   #      pref_bitter_food = p20616,
    #     pref_salty_food=p20716,
     #    pref_fatty_food=p20659)


## Food preference phenotypes ------------------------------

phenos_id <- fread("../data/processed/gwas/ukb_phenos_gwas_MA_prefs.txt") %>%
  rename(id=IID) %>% select(-FID)

allphenos_id <- left_join(phenos_id, base_id, by="id") %>% 
  left_join(educ_id, by="id") %>%
  left_join(income_id, by="id") %>%
  left_join(biochem_id, by="id") %>%
  left_join(t2d_id, by="id") %>%
  left_join(t2d_incid, by="id") %>%
  left_join(ffq_id, by="id") %>%
  left_join(prefs_id, by="id") %>% 
  left_join(diet24hr_id, by="id") %>%
  mutate(id=as.integer(id))


################################################
## Re-code phenotypes with descriptive labels ## 
################################################

smoke_labs <- list("Poor"=1, "Fair"=2, "Good"=3, "Excellent"=4)
health_labs <- list("Poor"=1, "Fair"=2, "Good"=3, "Excellent"=4)
alcohol_labs <- list(
  "Never"=0, "Special occasions only"=1,"One to three times a month"=2, "Once or twice a week"=3,
  "Three or four times a week"=4, "Daily or almost daily"=5)

covars_id <- fread("../data/processed/ukb_phenos_gwas_MA.txt") %>% mutate(
  sex.lab = case_when(sex==0~"Female", sex==1~"Male"),
  health.lab = descr_label.fun(., "health", health_labs, ordered=T),
  smoke_level.lab = factor(case_when(
    smoke_former==1~"Previous", smoke_current==1~"Current",
    smoke_former==0 & smoke_current==0 ~ "Never",
    TRUE~as.character(NA)), levels=c("Current", "Previous", "Never")),
  alch_level.lab = descr_label.fun(., "alcohol", alcohol_labs, ordered=T)) %>%
  select(-c(starts_with("pref_"), "ac")) %>%
  mutate_at("ancestry", ~ifelse(.=="","Missing",.))

allphenos_id <- left_join(allphenos_id, covars_id, by="id")


#########################
## Save as .rda & .csv ## 
#########################

#allphenos_id %>% fwrite("../data/processed/ukb_MA_allphenos.csv")
allphenos_id %>% saveRDS("../data/processed/ukb_MA_allphenos.rda")



###################################################
## Load PGS variables & subset to complete cases ## 
###################################################

# PGS in testing set (n=34599)
pgs_id <- reduce(lapply(
  grep("merged", list.files("../data/processed/pgs/freeze_prefs_adjBehav_20250327", full.names = T),value=T), fread),
  full_join, by="id")

analysis <- left_join(pgs_id, allphenos_id, by="id")

define_quantile <- function(pgs_var, quants=10, data=analysis) {
  x=(100/quants)/100
  data %>% select("pgs"=pgs_var) %>%
    mutate(pgs_quant=cut(pgs, breaks = quantile(pgs, probs = seq(0,1,x), na.rm = TRUE), 
                         include.lowest = F, labels = 1:quants)) %>%
    pull(pgs_quant)
}

analysis <- analysis %>% mutate(
  pref_sweet_food_pgs_decile=define_quantile("pref_sweet_food_pgs"),
  pref_bitter_food_pgs_decile=define_quantile("pref_bitter_food_pgs"),
  pref_salty_food_pgs_decile=define_quantile("pref_salty_food_pgs"),
  pref_fatty_food_pgs_decile=define_quantile("pref_fatty_food_pgs"),
  pref_sweet_food_pgs_quintile=define_quantile("pref_sweet_food_pgs", quants=5),
  pref_bitter_food_pgs_quintile=define_quantile("pref_bitter_food_pgs", quants=5),
  pref_salty_food_pgs_quintile=define_quantile("pref_salty_food_pgs", quants=5),
  pref_fatty_food_pgs_quintile=define_quantile("pref_fatty_food_pgs", quants=5)) %>%
  mutate(across(contains("decile"), ~as.factor(.))) %>%
  mutate(across(contains("quintile"), ~as.factor(.)))


analysis %>% fwrite("../data/processed/ukb_MA_analysis_testing.csv")
analysis %>% saveRDS("../data/processed/ukb_MA_analysis_testing.rda")


## End of file 

