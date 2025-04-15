# load packages

library(tidyverse)
library(data.table)

opt="../../opt"

# load functions from pantry package (cloned from GitHub :])
lapply(list.files(paste0(opt,"../../pantry/functions"), full.names=T), source)


## winsorize
winsorize <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x<bounds[1], bounds[1], ifelse(x>bounds[2], bounds[2], x))
  x
}


## From Cole et al., Heritability:
#A detailed 24HR questionnaire in which a subset of participants answered over 200 questions 
#on specific foods and beverages consumed (with quantities) in the preceding 24-hour day. 
#The 24HR was implemented as a questionnaire for the final 70 K in-person baseline assessment
#center participants from 2009–2010 and emailed four times to 320 K participants who consented 
#to re-contact via email between February 2011 and April 2012. Approximately 200 K individuals 
#have at least one and up to five recorded 24HR questionnaires.

#Each questionnaire was filtered for credible estimates of total energy intake 
#[≥1,000 kJ (UKB field 100002) and ≤20 MJ for males and ≤18 MJ for females (UKB field 100026)],
#typical dietary intake (UKB fields 100020 and 20085), completion duration greater than or 
#equal to 5 min (UKB field 20082), and overall completion (UKB field 20081). Additionally, 
#the participant could not be pregnant within 1 year of taking the 24HR nor have a cancer 
#diagnosis within the previous year (UKB fields 3,140 and 40005). All 24HR questions were 
#converted into 1/0 for yes/no to consumption; each categorical response was coded similarly
#[e.g., UKB field 20086 for special diet was converted into six binary variables for each response
#(gluten-free, lactose-free, low calorie, vegetarian, vegan, and a combined vegetarian or vegan field)]. 
#24HR questions pertaining to quantity consumed were also included as continuous variables.


#####################################
##  Load 24HR & helper data files  ##
#####################################

cat("Loading 24HR data ...\n")

path_to_24hr="/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_jul_2024/ukbb_app27892_diet_07242024.csv"

# 24hr data
diet24hr=fread(path_to_24hr) %>% rename(id = eid)

## ids for participants who withdrew consent
withdraw=scan("/humgen/florezlab/UKBB_app27892/withdraw/w27892_20241217.csv", what=character()) #488

diet24hr <- diet24hr %>% filter(!id %in% withdraw)

codebook=readxl::read_excel("../run/ukb_24hr_codebook.xlsx")



#########################
##  24HR data cleaning ##
#########################

## Used pre-calculated food group weights (g) by UKB, from 24-hour recalls, as described
# in https://link.springer.com/article/10.1007/s00394-021-02535-x (N=93 food groups)

## Data cleaning
# -restrict plausible total energy intakes (Cole: 1000 kj/day & <20,000kj [M] or <18000 kj [F])
# -restrict to typical_diet_yesterday
# -restrict to >= 5 min completion time
# -restrict to overall completion (yes)

# ======================================================================
## Restrictions: typical diet, plausible total energy, complete data
# ======================================================================

diet_valid_mean <- function(field, df) {
  # For a given nutrient:
  # - Select all columns for that nutrient
  # - Tabule number of instances (not missing)
  # - if >=2 instances, calculate mean; if not, code as NA
  diet_fields_df <- lapply(0:4, function(i) {
    field_id <- paste0("p", field,"_i", i) # diet field 
    valid_field <- paste0("p100026_i", i) # ukb field for valid 24hr (energy & sex)
    typical_field <- paste0("p100020_i", i) # ukb field for typical diet yesterday
    reason_field <- paste0("p20085_i", i) # ukb field for reason for atypical diet
    duration_field <- paste0("p20082_i", i) # ukb field for time duration of 24hr
    valid_24hr <- df %>% select(id, field=field_id, valid=valid_field, typical=typical_field, 
                                reason=reason_field, duration=duration_field) %>%
      mutate(field=ifelse(valid == "" & typical != "No" & reason == "" & !duration<5, field, NA)) %>% # replace with NA, if invalid 24hr 
      select(id, field_id=field) %>% rename_with(., ~gsub("field_id", field_id, .))
    valid_24hr
  }) %>% reduce(full_join, by = "id") %>%
    mutate(n_valid_24hr = rowSums(!is.na(across(paste0("p", field, "_i",0:4))))) %>%  # count n of valid 24hr
    mutate(field_mean=ifelse(n_valid_24hr>0, rowSums(across(paste0("p", field, "_i",0:4)), na.rm=T)/n_valid_24hr, NA)) %>%
    select(id, field_mean) %>%
    rename_with(., ~gsub("field_mean", paste0("p", field, "_mean"), .)) %>%
    pull(ends_with("mean"))
}

# Run over all 93 diet group fields
diet_fields <- codebook %>% pull(Field.ID)
diet_means_raw <- lapply(diet_fields, function(f) {
  diet24hr %>% mutate(
    field_mean=diet_valid_mean(f, .)) %>%
    rename_with(., ~gsub("field", paste0("p", f), .)) %>%
    select(id, ends_with("mean"))
  }) %>% reduce(full_join, by = "id")


#diet_means %>% fwrite("../data/processed/ukb_diet_mean24hr_raw.txt", sep="\t")

diet_means <- diet_means_raw %>%
  mutate(across(-id, ~winsorize(., SDs=5)))


# ======================================================================
## Calculate intake 13 food (g) and 2 beverage (mL) groups 
# ======================================================================

FoodGroups <- unique(codebook$FoodGroup.Name)


build_foodgroup <- function(food_group, df) {
  # For a given nutrient:
  # - Select all columns for that nutrient
  # - Tabule number of instances (not missing)
  # - if >=2 instances, calculate mean; if not, code as NA
  diet_vars_id <- codebook$Var.ID[codebook$FoodGroup.Name==food_group]
  diet_vars_id <- paste0(diet_vars_id, "_mean") # diet_id_mean 
  vars_dat <- df %>% select(id, all_of(diet_vars_id)) %>%
    mutate(fg_dat = rowSums(across(all_of(diet_vars_id)), na.rm=T)) %>%
    select(id, fg_dat) %>%
    rename_with(., ~gsub("fg_dat", paste0(food_group, "_mean"), .)) %>%
    pull(ends_with("mean"))
  return(vars_dat)
}

diet_means <- diet_means %>% mutate(
  fg_bevs_alcoholic=build_foodgroup("fg_bevs_alcoholic", .),
  fg_cereals=build_foodgroup("fg_cereals", .),
  df_dairy_products=build_foodgroup("df_dairy_products", .),
  
  fg_eggs=build_foodgroup("fg_eggs", .),
  fg_fats_spreads=build_foodgroup("fg_fats_spreads", .),
  fg_fish_dishes=build_foodgroup("fg_fish_dishes", .),
  
  fg_fruits=build_foodgroup("fg_fruits", .),
  fg_meat_products=build_foodgroup("fg_meat_products", .),
  fg_meat_substitutes=build_foodgroup("fg_meat_substitutes", .),
  
  fg_mixed_dishes=build_foodgroup("fg_mixed_dishes", .),
  fg_bevs_nonalch=build_foodgroup("fg_bevs_nonalch", .),
  fg_nuts_seeds=build_foodgroup("fg_nuts_seeds", .),
  
  fg_condiments=build_foodgroup("fg_condiments", .),
  fg_sweets_snacks=build_foodgroup("fg_sweets_snacks", .),
  fg_vegetables_potatoes=build_foodgroup("fg_vegetables_potatoes", .)
)


diet_means %>% fwrite("../data/processed/ukb_diet24hr.csv")
diet_means %>% saveRDS("../data/processed/ukb_diet24hr.rda")





  
# Write function to build food groups
build_foodgroup <- function(food_group_data, food_group, codebook=codebook) {
  # For a sub_group in the DHD15,
  # - Use the diet_codebook to determine list of included food items (foodIDs)
  # - Sum over all food items to create a single food group
  food_items <- c(paste0(with(codebook, Var.ID[which(FoodGroup.ID == food_group)]), "_mean"))
  food_items.df <- diet_means %>% select(c(all_of(food_items)))
  food_group.df <- rowSums(mutate_all(food_items.df, ~as.numeric(.)))
  return(food_group.df + food_group_data)
}

# Build food groups 
food_groups.dat <- diet_means %>% bind_cols(
  as.data.frame(matrix(0, nrow(diet_means), length(FoodGroups), dimnames = list(NULL, FoodGroups)))) %>%
  mutate(across(c(FoodGroups), ~build_foodgroup(., cur_column()))) %>%
  #rename_at(c(FoodGroups), ~paste0(., "_gday")) %>% 
  select("id", all_of(starts_with("fg_")))


## Convert fats/oils from sev/wk (1 sev = 1 gram) to GRAMS per DAY
#foodgroups_swk.dat <- foodgroups_swk.dat %>% mutate(across(starts_with("dhdfg_fat"), ~ ./7)) %>%
# rename_with(., ~gsub("_swk", "_gday", .), starts_with("dhdfg_fat"))

food_groups_swk <- names(foodgroups_swk.dat %>% select(-"id"))


         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         