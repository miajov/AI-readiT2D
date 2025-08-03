library(dplyr)
library(tidyr)
library(stringr)
library(tidyLPA)

##Loading PAID survey

matched_ppts_ids <- read.csv("~/matched107_w_wearable_data_ppts.csv", header = TRUE) 
length(matched_ppts_ids$id)
survey <- read.csv("~/Downloads/observation(in) (2).csv", header = TRUE)  %>%
  filter(str_detect(observation_id, "paid")) %>% 
  separate(
    col = observation_id,
    into = paste0("V", 1:25),  # You can increase 25 if you need more columns
    sep = ",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)",  # Splits commas outside quotes only
    fill = "right",
    extra = "merge"
  )

survey_split$id =as.numeric(survey_split$V2)
survey_split$response =  as.numeric(survey_split$V8)
survey_split$survey = survey_split$V15

survey_paid = survey_split %>% 
  select(id, response, survey)

# Assuming your dataframe is called df
df_clean <- survey_paid %>%
  mutate(variable = str_extract(survey, "^[^,]+")) %>%
  filter(!is.na(variable), !variable %in% c("paidstartts", "paidcmpts")) %>%
  # Pivot to wide format: one row per person, one column per PAID variable
  pivot_wider(
    id_cols = id,                # One row per participant
    names_from = variable,       # Column names = PAID question shortcodes
    values_from = response       # Values = response scores
  )


#######

edgeweights <- read.csv("~/df_edgeweights_clustering.csv", header = TRUE) 
edgeweights$id = edgeweights$subject_id
matched_ppts_ids <- read.csv("~/matched107_w_wearable_data_ppts.csv", header = TRUE) 
length(matched_ppts_ids$id)

##### avg steps 
daily_steps <- pa %>%
  group_by(participant_id, day) %>%
  summarise(daily_total_steps = sum(total_steps, na.rm = TRUE), .groups = "drop")
# Step 2: Compute average daily steps per person
avg_steps_per_person <- daily_steps %>%
  group_by(participant_id) %>%
  summarise(avg_daily_steps = mean(daily_total_steps, na.rm = TRUE), .groups = "drop")
avg_steps_per_person$id = avg_steps_per_person$participant_id

### avg sleep 
sleep <- read.csv("~/all_participants_sleep_data.csv", header = TRUE) 
head(sleep)
mean(sleep$avg_total_sleep_hours)
sleep$id= sleep$user_id
sleep$id <- as.numeric(gsub("AIREADI-", "", sleep$id))

sleep_df  = sleep %>%
  distinct(id, .keep_all = T)
sleep_df = sleep_df %>% 
  select(id, avg_total_sleep_hours)

### avg stress 
stress <- read.csv("~/stress-data_hourly-scirep.csv", header = TRUE) 
head(stress)
daily_stress <- stress %>%
  group_by(participant_id, day) %>%
  summarise(daily_mean_stress = mean(mean_stress, na.rm = TRUE), .groups = "drop")

# Step 2: Compute average daily stress per person
avg_stress_per_person <- daily_stress %>%
  group_by(participant_id) %>%
  summarise(avg_daily_stress = mean(daily_mean_stress, na.rm = TRUE), .groups = "drop")

avg_stress_per_person$id =avg_stress_per_person$participant_id

### merge all dfs 
df= merge(edgeweights, df_clean, by="id")
length(df$id)
df= merge(matched_ppts,df, by="id")
length(df$id)
df= merge(df, avg_steps_per_person, by="id")
length(df$id)
df= merge(df, sleep_df, by="id")
length(df$id)
# [1] 91

df= merge(avg_stress_per_person, df, by="id")

####### creating LPA profiles 

lpa = df %>%
  select(
    age,
    hbA1c,
    mean_stress_body_cgm_value,
    body_cgm_value_mean_stress,
    total_steps_body_cgm_value,
    body_cgm_value_total_steps
  ) %>%
  scale() %>%                         # Standardize variables
  as.data.frame() %>%                # Convert back to data frame
  single_imputation() %>%            # Impute missing values (if needed)
  estimate_profiles(2)              #estimate 2 profiles 

lpadata =get_data(lpa)$Class
paiddata =get_data(paid)$Class


lpadata_df <- data.frame(
  id = df$id,
  lpa_class = get_data(lpa)$Class
)

###recode DV complications variable
df <- df %>%
  mutate(paid_cml_binary = if_else(X.paid_cml >= 1, 1, 0))

##run glm model 

library(pscl)
model <- glm(paid_cml_binary ~ glycemic_control, data = df, family = binomial)
model1 <- glm(paid_cml_binary ~ glycemic_control + lpadata_df$lpa_class, data = df, family = binomial)

pR2(model)
pR2(model1)

anova(model1, model, test = "Chisq")

# Fit the two models
model1 <- glm(df$paid_cml_binary ~ df$glycemic_control, family = binomial)
model2 <- glm(df$paid_cml_binary ~ lpadata_df$lpa_class + df$glycemic_control, family = binomial)

# Function to extract odds ratios and CIs
get_or_ci <- function(model) {
  OR <- exp(coef(model))
  CI <- exp(confint(model))  # Profile likelihood CIs
  result <- cbind(OR, CI)
  colnames(result) <- c("OR", "2.5 %", "97.5 %")
  round(result, 3)
}

# Display ORs and CIs for both models
get_or_ci(model1)
get_or_ci(model2)

# --- 2. Compare AICs ---
aic_model1 <- AIC(model1)
aic_model2 <- AIC(model2)

cat("\nModel 1 AIC:", aic_model1, "\n")
cat("Model 2 AIC:", aic_model2, "\n")

# --- 3. Likelihood Ratio Test ---
lrt_result <- anova(model1, model2, test = "LRT")
cat("\nLikelihood Ratio Test:\n")
print(lrt_result)\


ibrary(pscl)

# Fit both models again (if not already fitted)
model1 <- glm(df$paid_cml_binary ~ df$glycemic_control, family = binomial)
model2 <- glm(df$paid_cml_binary ~ lpadata_df$lpa_class + df$glycemic_control, family = binomial)

# Get McFadden's pseudo R²
r2_model1 <- pR2(model1)["McFadden"]
r2_model2 <- pR2(model2)["McFadden"]

# Print results
cat("McFadden's pseudo R²:\n")
cat("Model 1 (Glycemic Control Only):", round(r2_model1, 4), "\n")
cat("Model 2 (Glycemic Control + LPA Profile):", round(r2_model2, 4), "\n")

coefs <- summary(model2)$coefficients

# Calculate odds ratios
OR <- exp(coefs[, "Estimate"])
# Calculate 95% confidence intervals for OR
lower_CI <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
upper_CI <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])

# Combine into a data frame
results <- data.frame(
  Predictor = rownames(coefs),
  Estimate = coefs[, "Estimate"],
  Std_Error = coefs[, "Std. Error"],
  OR = round(OR, 3),
  CI_Lower = round(lower_CI, 3),
  CI_Upper = round(upper_CI, 3),
  p_value = coefs[, "Pr(>|z|)"]
)

print(results)
