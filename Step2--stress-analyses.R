library(jsonlite)
library(dplyr)
library(tidyr)
library(lubridate)
library(lme4)

# Load metadata
matched_ppts <- read.csv("~/matched107_w_wearable_data_ppts.csv", header = TRUE)

# List participant folders from directory
base_dir <- "~/Downloads/OneDrive_stress_works/"
participant_folders <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
participant_ids <- participant_folders

# Process each participant and save individual hourly stress
for (participant_id in participant_ids) {
  message("Processing participant: ", participant_id)
  
  participant_folder <- file.path(base_dir, participant_id)
  json_file <- file.path(participant_folder, paste0(participant_id, "_stress.json"))
  
  if (!file.exists(json_file)) {
    message("  → JSON file not found - skipping.")
    next
  }
  
  json_data <- fromJSON(json_file, flatten = TRUE)
  
  if (is.null(json_data$body$stress)) {
    message("  → No stress data - skipping.")
    next
  }
  
  stress_data <- json_data$body$stress %>%
    mutate(
      user_id = json_data$header$user_id,
      uuid = json_data$header$uuid,
      creation_date = json_data$header$creation_date_time,
      stress_value = stress.value,
      stress_time = ymd_hms(effective_time_frame.date_time)
    ) %>%
    select(user_id, uuid, creation_date, stress_time, stress_value) %>%
    filter(!is.na(stress_value))
  
  hourly_stress <- stress_data %>%
    mutate(hour = floor_date(stress_time, unit = "hour")) %>%
    group_by(user_id, hour) %>%
    summarise(mean_stress = mean(stress_value, na.rm = TRUE), .groups = "drop") %>%
    arrange(hour) %>%
    mutate(hour_number = row_number())
  
  saveRDS(hourly_stress, file = file.path(participant_folder, paste0(participant_id, "_hourly_stress.rds")))
  message("  → Processed and saved hourly stress for participant ", participant_id)
}

# Combine all participants
all_hourly_stress_list <- list()

for (participant_id in participant_ids) {
  participant_folder <- file.path(base_dir, participant_id)
  rds_file <- file.path(participant_folder, paste0(participant_id, "_hourly_stress.rds"))
  
  if (file.exists(rds_file)) {
    df <- readRDS(rds_file)
    df$participant_id <- participant_id
    all_hourly_stress_list[[participant_id]] <- df
  }
}

combined_hourly_stress <- bind_rows(all_hourly_stress_list)

## inspect missing data
## to plot descriptives 

# Remove zero-stress or invalid days
zero_stress_days <- combined_hourly_stress %>%
  mutate(day = as.Date(hour)) %>%
  group_by(participant_id, day) %>%
  summarise(avg_stress = mean(mean_stress, na.rm = TRUE), .groups = "drop") %>%
  filter(is.na(avg_stress) | avg_stress == 0)

cleaned_stress_data <- combined_hourly_stress %>%
  mutate(day = as.Date(hour)) %>%
  anti_join(zero_stress_days, by = c("participant_id", "day")) %>%
  group_by(participant_id) %>%
  arrange(participant_id, day, hour) %>%
  mutate(
    day_number = dense_rank(day),
    hour_number = row_number()
  ) %>%
  ungroup()

cleaned_stress_data$id <- gsub("AIREADI-", "", cleaned_stress_data$user_id)

# Merge with matched_ppts for modeling, if desired
df <- merge(matched_ppts, cleaned_stress_data, by = "id")
#write.csv(df, "stress-data_hourly.csv")

length(unique(df$id)) #100


##filter missing data
df_filtered <- df[df$mean_stress >= 0, ]
df_filtered <- df_filtered %>%
  group_by(id) %>%
  mutate(
    date_only = as.Date(hour),
    day_number_updated = as.numeric(date_only - min(date_only)) + 1
  ) %>%
  ungroup() %>%
  filter(day_number <= 16) ## remove ppts who have more than 16 days. 


# Step 1: Aggregate stress per day per person
daily_stress <- df %>%
  group_by(id, day_number, glycemic_control) %>%
  summarise(daily_mean_stress = mean(mean_stress, na.rm = TRUE), .groups = "drop")

# Step 2: Aggregate to person-level mean stress
person_stress <- daily_stress %>%
  group_by(id, glycemic_control) %>%
  summarise(mean_stress = mean(daily_mean_stress, na.rm = TRUE),
            n_days = n(),
            .groups = "drop")

##descriptives
round(sd(person_stress$mean_stress[person_stress$glycemic_control == 'poor']), digits=2)
round(sd(person_stress$mean_stress[person_stress$glycemic_control == 'good']), digits=2)

###more than 5 days of data
complete_ids <- daily_stress %>%
  group_by(id, glycemic_control) %>%
  summarise(n_days = n_distinct(day_number), .groups = "drop") %>%
  filter(n_days > 5)

###descriptives 
# # Step 2: Filter full data for participants
plot_data <- daily_stress %>%
  semi_join(complete_ids, by = c("id", "glycemic_control")) %>%
  mutate(id = as.factor(id))

plot_data <- plot_data %>% mutate(id = as.factor(id))

# Just to check:
print(head(plot_data))
print(names(plot_data))


plot_data_10 <- plot_data %>%
  filter(day_number <= 10)

group_means_stress <- plot_data %>%
  group_by(glycemic_control, day_number) %>%
  summarise(mean_stress = mean(daily_mean_stress), .groups = "drop")

group_means_10 <- group_means_stress %>%
  filter(day_number <= 10)



#######

#### plot one main paper

set.seed(42)  # for reproducibility
# Sample 30 unique IDs per group
sampled_ids <- plot_data_10 %>%
  distinct(id, glycemic_control) %>%
  group_by(glycemic_control) %>%
  slice_sample(n = 30) %>%
  ungroup()

# Filter plot_data_10 to keep only those sampled participants
plot_data_10_subset <- plot_data_10 %>%
  semi_join(sampled_ids, by = c("id", "glycemic_control"))

# Recalculate group means for the subset
group_means_10_subset <- plot_data_10_subset %>%
  group_by(glycemic_control, day_number) %>%
  summarise(mean_stress = mean(daily_mean_stress), .groups = "drop")

# Now plot with the subsetted data

glycemic_labels <- c(
  "good" = "Controlled",
  "poor" = "Uncontrolled"
)

####main paper plot 1 

ggplot(plot_data_10_subset, aes(x = day_number, y =  daily_mean_stress, group = id, color = glycemic_control)) +
  geom_line(alpha = 0.3, linewidth = 0.7) +  # individual participant lines thicker
  geom_line(data = group_means_10_subset,
            aes(x = day_number, y = mean_stress, color = glycemic_control, group = glycemic_control),
            linewidth = 1.8) +  # bold group mean line
  facet_wrap(~glycemic_control, ncol = 2, labeller = labeller(glycemic_control = glycemic_labels)) +
  scale_color_manual(values = okabe_ito[c("poor", "good")]) +
  labs(
    title = "Daily Step Counts Over First 10 Days by Glycemic Control Group (30 Random Participants per Group)",
    x = "Day",
    y = "Daily stress",
    color = "Glycemic Control"
  ) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

###


####supplement plot 1

# Subset to first 10 days
plot_data_subset <- plot_data %>% filter(day_number <= 10)
group_means_stress_subset <- group_means_stress %>% filter(day_number <= 10)

# Plot
ggplot(plot_data_subset, aes(x = day_number, y = daily_mean_stress, group = id, color = glycemic_control)) +
  geom_line(alpha = 0.3, linewidth = 0.7) +
  geom_line(data = group_means_stress_subset,
            aes(x = day_number, y = mean_stress, color = glycemic_control, group = glycemic_control),
            linewidth = 1.8) +
  facet_wrap(~glycemic_control, ncol = 2, labeller = labeller(glycemic_control = glycemic_labels)) +
  scale_color_manual(values = okabe_ito[c("poor", "good")]) +
  labs(
    title = "Daily Stress Over First 10 Days by Glycemic Control Group (30 Random Participants per Group)",
    x = "Day",
    y = "Daily Stress",
    color = "Glycemic Control"
  ) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")


##likelihood ratio tests 

##align dfs 
df_sub <- df_filtered %>%
  filter(!is.na(mean_stress), !is.na(glycemic_control), !is.na(id), !is.na(day_number))

# Fit models on the same data
model_full <- lmer(
  mean_stress ~ glycemic_control +
    (1 | id) +
    (1 | id:day_number),
  data = df_sub,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

model_reduced <- lmer(
  mean_stress ~ glycemic_control + 
    (1 | id), 
  data = df_sub,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

#  likelihood ratio test
lrt <- anova(model_reduced, model_full)


###modeling 

model_stress <- lmer(
  mean_stress ~ 1 + glycemic_control +
    (1 | id) +
    (1 | id:day_number),
  data = df,
  REML = TRUE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)


# Extract variance components
vc <- as.data.frame(VarCorr(model_stress))
var_id <- vc$vcov[vc$grp == "id"]           # between-person variance
var_day <- vc$vcov[vc$grp == "id:day_number"]  # day-to-day variance within person
var_res <- vc$vcov[vc$grp == "Residual"]    # residual variance

var_total <- var_id + var_day + var_res

# Calculate ICC (proportion of variance due to between-person differences)
icc <- var_id / var_total

# Estimate between-group variance explained by fixed effect glycemic_control using marginal R2
r2_full <- r2(model_stress)
r2_marginal <- r2_full$R2_marginal
var_between_group <- r2_marginal * var_total

# Average within-group variance (between-person + day-to-day)
var_within_group_avg <- var_id + var_day

# Calculate ratio of within-group variance to between-group variance
ratio_within_to_between <- var_within_group_avg / var_between_group

# Print results
cat("Between-group variance (fixed effect glycemic_control):", round(var_between_group, 2), "\n")
cat("Between-person variance:", round(var_id, 2), "\n")
cat("Day-to-day variance:", round(var_day, 2), "\n")
cat("Residual variance:", round(var_res, 2), "\n")
cat("Total variance:", round(var_total, 2), "\n")
cat("ICC (between-person variance / total variance):", round(icc, 3), "\n")
cat("Average within-group variance (between-person + day-to-day):", round(var_within_group_avg, 2), "\n")
cat("Ratio (within-group / between-group variance):", round(ratio_within_to_between, 2), "\n")


### main
r2_full <- r2(model_stress)  # your full model with glycemic_control fixed effect
r2_marginal <- r2_full$R2_marginal
var_id <- vc$vcov[vc$grp == "id"]
var_day <- vc$vcov[vc$grp == "id:day_number"]
var_res <- vc$vcov[vc$grp == "Residual"]
var_total <- var_id + var_day + var_res
var_between_group <- r2_marginal * var_total


