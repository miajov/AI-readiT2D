
##### step. 1 matching ppts and propensity scoring

library(dplyr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(MatchIt)

####################################
###load ppts
ppts <- read.csv("~/Downloads/participants.tsv", sep="\t", header=TRUE) #change path 
length(ppts$participant_id)

###only T2D patients
ppts <- ppts[ppts$study_group %in% c("oral_medication_and_or_non_insulin_injectable_medication_controlled", "insulin_dependent"), ]
length(ppts$participant_id) 

##check how many on insluin
table(ppts$study_group)                       

########################################################

### merge with HBA1C blood data
ppts$id = ppts$participant_id
blood <- read.csv("~/Downloads/measurement.csv", header=TRUE)
blood$id = blood$person_id
merged = merge(blood,ppts, by="id")
biomarkers_df <- merged %>%
  filter(grepl("HbA1c", measurement_source_value))

biomarkers_df$value_source_value = as.numeric(biomarkers_df$value_source_value)
biomarkers_df$hbA1c = biomarkers_df$value_source_value
length(unique(biomarkers_df$person_id)) 
hist(biomarkers_df$hbA1c)
summary(biomarkers_df$hbA1c)
biomarkers_df <- biomarkers_df[!is.na(biomarkers_df$hbA1c), ]

#########

### merge with demo 
age <- read.csv("~/Downloads/person(in).csv", header=TRUE)
age$year_of_birth
age$id = age$person_id 

biomarkers_df= merge(biomarkers_df, age, by="id")

# Add age column
biomarkers_df <- biomarkers_df %>%
  mutate(age = 2025 - year_of_birth)

# Group by HbA1c control and summarize age
age_summary <- biomarkers_df %>%
  group_by(hbA1c_control) %>%
  summarise(
    count = n(),
    mean_age = mean(age, na.rm = TRUE),
    median_age = median(age, na.rm = TRUE),
    sd_age = sd(age, na.rm = TRUE),
    min_age = min(age, na.rm = TRUE),
    max_age = max(age, na.rm = TRUE),
    q25 = quantile(age, 0.25, na.rm = TRUE),
    q75 = quantile(age, 0.75, na.rm = TRUE)
  )

print(age_summary)

### merge with demo (medical conditions)

## conditions
condition_occ <- read.csv("~/Downloads/condition_occurrence(in).csv", header=FALSE)
# Remove header row first 
data_no_header <- condition_occ[-1, , drop=FALSE]
# Split V1 into multiple columns
condition_occ_sep <- data_no_header %>%
  separate(col = V1, 
           into = paste0("V", 1:16),  # Adjust number of columns expected
           sep = ",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)",  # Regex to split on commas outside quotes
           extra = "merge",
           fill = "right")
condition_occ_sep$id = condition_occ_sep$V2
condition_occ_sep$comorbidity = condition_occ_sep$V14

cond = condition_occ_sep %>% 
  select(id,comorbidity)

cond_clean <- cond %>%
  # Remove the surrounding quotes from comorbidity first
  mutate(comorbidity = str_replace_all(comorbidity, '^"|"$', '')) %>%
  # Then separate by the first comma into two new columns
  separate(comorbidity, into = c("code", "description"), sep = ",", extra = "merge", fill = "right") %>%
  # Trim whitespace in the new columns
  mutate(across(c(code, description), str_trim))
table(cond_clean$code)

# Vector of codes that require HbA1c target adjustment (less stringent target)
adjust_codes <- c(
  "mhoccur_mi",   # Heart attack
  "mhoccur_oa",   # Other heart issues
  "mhoccur_strk", # Stroke
  "mhoccur_circ", # Circulation problems
  "mhoccur_hbp",  # High blood pressure
  "mhoccur_pd",   # Parkinson's disease
  "mhoccur_ded",  # Dementia
  "mhoccur_cogn", # Mild cognitive impairment
  "mhoccur_ms",   # Multiple sclerosis
  "mhoccur_cns",  # Other neurological conditions
  "mhoccur_ca",   # Cancer
  "mhoccur_plm",  # Chronic pulmonary (lung) problems
  "mhoccur_rnl"   # Kidney problems
)

# Create binary variable 1 = adjust HbA1c, 0 = no adjustment
cond_clean <- cond_clean %>%
  mutate(
    hba1c_adjustment = ifelse(code %in% adjust_codes, 1, 0)
  )


# Severe comorbidities: more relaxed target (<8.0%)
severe_codes <- c(
  "mhoccur_ded",    # Dementia
  "mhoccur_ca",     # Cancer
  "mhoccur_plm",    # Chronic pulmonary problems
  "mhoccur_rnl",    # Kidney problems
  "mhoccur_mi",     # Heart attack
  "mhoccur_strk",   # Stroke
  "mhoccur_pd",     # Parkinson's disease
  "mhoccur_cns"     # Other neurological conditions
)

# Mild comorbidities: moderate adjustment (<7.5%)
mild_codes <- c(
  "mhoccur_hbp",    # High blood pressure
  "mhoccur_cogn",   # Mild cognitive impairment
  "mhoccur_oa",     # Other heart issues
  "mhoccur_circ",   # Circulation problems
  "mhoccur_ms"      # Multiple sclerosis
)

# Label each condition severity
cond_clean <- cond_clean %>%
  mutate(
    severity = case_when(
      code %in% severe_codes ~ "severe",
      code %in% mild_codes ~ "mild",
      TRUE ~ "none"
    )
  )

# Summarize counts of total, severe, and mild comorbidities per person
comorbidity_counts <- cond_clean %>%
  #filter(severity != "none") %>%  # only count relevant comorbidities
  group_by(id) %>%
  summarize(
    total_comorbidities = n(),
    severe_comorbidities = sum(severity == "severe"),
    mild_comorbidities = sum(severity == "mild")
  )

# Join back to original data if needed:
cond_clean <- cond_clean %>%
  left_join(comorbidity_counts, by = "id")

# Get all unique person ids
all_ids <- cond_clean %>%
  distinct(id)

# Summarize counts as before, but from full data (not filtered)
comorbidity_counts <- cond_clean %>%
  group_by(id) %>%
  summarize(
    total_comorbidities = sum(severity != "none"),
    severe_comorbidities = sum(severity == "severe"),
    mild_comorbidities = sum(severity == "mild")
  )

# Join with all_ids to ensure everyone is included
# Get all unique person ids
all_ids <- cond_clean %>%
  distinct(id)

# Summarize counts as before, but from full data (not filtered)
comorbidity_counts <- cond_clean %>%
  group_by(id) %>%
  summarize(
    total_comorbidities = sum(severity != "none"),
    severe_comorbidities = sum(severity == "severe"),
    mild_comorbidities = sum(severity == "mild")
  )

# Join with all_ids to ensure everyone is included
comorbidity_counts_full <- all_ids %>%
  left_join(comorbidity_counts, by = "id") %>%
  mutate(
    total_comorbidities = replace_na(total_comorbidities, 0),
    severe_comorbidities = replace_na(severe_comorbidities, 0),
    mild_comorbidities = replace_na(mild_comorbidities, 0)
  )

#### merge back final dfs

df= merge(comorbidity_counts_full, biomarkers_df, by="id")
df = df %>% select(id, total_comorbidities,
                   severe_comorbidities,
                   mild_comorbidities, age, 
                   hbA1c, measurement_date)

length(unique(df$id))
head(df)

df <- df %>%
  mutate(
    hba1c_target = case_when(
      severe_comorbidities > 0 ~ "< 8.0%",  # any severe comorbidity -> relax target
      age >= 65 & mild_comorbidities > 0  ~ "< 7.5%",  # older + mild comorbidity
      age < 65 & mild_comorbidities > 0 ~ "< 7.5%",  # optionally treat younger with mild also mild
      TRUE ~ "< 7.0%"  # default tighter target (no severe or mild comorbidity)
    )
  )

df <- df %>%
  mutate(
    hba1c_target = case_when(
      severe_comorbidities >= 2 ~ "< 8.5%",  # very complex health
      severe_comorbidities == 1 ~ "< 8.0%",  # complex health
      age >= 65 & mild_comorbidities >= 2 ~ "< 7.5%",  # older + multiple mild
      age >= 65 & mild_comorbidities == 1 ~ "< 7.0%",  # older + one mild
      age < 65 & mild_comorbidities >= 1 ~ "< 7.0%",  # younger + any mild
      age < 65 & mild_comorbidities == 0 ~ "< 7.0%",  # tightest for healthy younger
      age >= 65 & mild_comorbidities == 0 ~ "< 7.0%",  # healthy older
      TRUE ~ "< 7.0%"  # fallback
    )
  )

  df <- df %>%
  dplyr::mutate(
    # Extract numeric target threshold from hba1c_target string, e.g. "< 7.0%" -> 7.0
    target_numeric = as.numeric(str_extract(hba1c_target, "[0-9.]+")),
    meets_target = (hbA1c < target_numeric))

df$target_numeric = as.numeric(str_extract(df$hba1c_target, "[0-9.]+"))
df$meets_target = df$hbA1c < df$target_numeric

# Summary counts
summary_table <- df %>%
  group_by(hba1c_target) %>%
  dplyr::summarise(
    n = n(),
    met_target = sum(meets_target),
    not_met_target = n - met_target,
    pct_met = round(100 * met_target / n, 1)
  )

df <- df %>%
  mutate(
    glycemic_control = ifelse(meets_target, "good", "poor")
  )

table(df$glycemic_control)


### demographics 
ppts_med_group <- read.csv("~/Downloads/participants.tsv", sep="\t", header=TRUE) %>% 
  select(participant_id, study_group)
ppts_med_group$id = ppts_med_group$participant_id
df_m = merge(df, ppts_med_group, by="id")

df_m %>%
  filter(glycemic_control == "poor") %>%
  group_by(study_group) %>%
  summarise(
    n = n(),
    percent = round(100 * n / sum(n), 1),
    mean_age = round(mean(age, na.rm = TRUE), 1),
    sd_age = round(sd(age, na.rm = TRUE), 1),
    min_age = min(age, na.rm = TRUE),
    max_age = max(age, na.rm = TRUE)
  )


#### filtering to only those who have available wearable data 

#change path 
ppts_ids <- read.csv("~/Downloads/participants.tsv", sep="\t", header=TRUE) %>% 
  select(participant_id, wearable_activity_monitor)
ppts_ids$id =  ppts_ids$participant_id
filtere_w_have =merge(ppts_ids, df_m, by="id")
filtered_true <- filtere_w_have %>%
  filter(wearable_activity_monitor == TRUE)
sum(filtered_true$glycemic_control == "poor", na.rm = TRUE) #54



###### propensity matching


filtered_true <- filtered_true %>%
  mutate(glycemic_control_binary = ifelse(glycemic_control == "poor", 1, 0))

match_model <- matchit(
  glycemic_control_binary ~ age + study_group + total_comorbidities + severe_comorbidities + mild_comorbidities,
  data = filtered_true,
  method = "nearest",    # nearest-neighbor matching
  ratio = 1              # 1:1 matching
)

matched_df <- match.data(match_model)
table(matched_df$glycemic_control)  
summary(matched_df)

length(unique(matched_df$id))

fisher.test(table(matched_df$study_group, matched_df$glycemic_control))
table(matched_df$study_group, matched_df$glycemic_control)

wilcox.test(age ~ glycemic_control, data = matched_df)
wilcox.test(total_comorbidities ~ glycemic_control, data = matched_df)
wilcox.test(severe_comorbidities ~ glycemic_control, data = matched_df)
wilcox.test(mild_comorbidities ~ glycemic_control, data = matched_df)

summary(match_model, standardize = TRUE)


length(match_model$)
#write.csv(matched_df, "matched108_w_wearable_data_ppts.csv")


