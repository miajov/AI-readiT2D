library(dplyr)
library(lubridate)
library(jsonlite)
library(tidyr)
library(stringr)
library(purrr)  
library(lme4)
library(performance)


###included ppts 
matched_ppts <- read.csv("~/matched107_w_wearable_data_ppts.csv", header = TRUE)

# Base directory
base_dir <- "~/Downloads/garmin_vivosmart5 3/"
participant_folders <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

all_participants_data <- list()

for (folder in participant_folders) {
  sleep_files <- list.files(folder, pattern = "_sleep.json$", full.names = TRUE)
  
  for (json_file in sleep_files) {
    participant_id <- basename(json_file) %>% str_remove("_sleep.json")
    json_data <- fromJSON(json_file, flatten = TRUE)
    
    # Extract and prepare sleep data
    sleep_data <- json_data$body$sleep %>%
      bind_rows() %>%
      mutate(
        user_id = json_data$header$user_id,
        start_time = ymd_hms(sleep_stage_time_frame.time_interval.start_date_time),
        end_time = ymd_hms(sleep_stage_time_frame.time_interval.end_date_time),
        date = as.Date(start_time),
        sleep_stage = sleep_stage_state
      ) %>%
      arrange(user_id, date, start_time)
    
    # Remove duplicate or overlapping intervals
    sleep_data <- sleep_data %>%
      group_by(user_id, date) %>%
      arrange(start_time, .by_group = TRUE) %>%
      mutate(
        next_start = lead(start_time),
        curr_end = end_time,
        overlaps_with_next = !is.na(next_start) & curr_end > next_start
      ) %>%
      filter(!overlaps_with_next | is.na(overlaps_with_next)) %>%
      ungroup() %>%
      select(-next_start, -curr_end, -overlaps_with_next)
    
    # Recalculate sleep durations
    sleep_data <- sleep_data %>%
      mutate(
        sleep_duration = as.numeric(difftime(end_time, start_time, units = "mins"))
      )
    
    # Summarize per night
    daily_sleep <- sleep_data %>%
      group_by(user_id, date) %>%
      summarise(
        bed_time = min(start_time),
        wake_time = max(end_time),
        total_sleep_time_mins = sum(sleep_duration, na.rm = TRUE),
        total_sleep_time_hours = total_sleep_time_mins / 60,
        awake_duration = sum(sleep_duration[sleep_stage == "awake"], na.rm = TRUE) / 60,
        deep_duration = sum(sleep_duration[sleep_stage == "deep"], na.rm = TRUE) / 60,
        rem_duration = sum(sleep_duration[sleep_stage == "rem"], na.rm = TRUE) / 60,
        light_duration = sum(sleep_duration[sleep_stage == "light"], na.rm = TRUE) / 60,
        stage_sleep_time = sum(sleep_duration[sleep_stage %in% c("light", "deep", "rem")], na.rm = TRUE) / 60,
        .groups = "drop"
      )
    
    # Summarize per participant
    avg_sleep_summary <- daily_sleep %>%
      group_by(user_id) %>%
      summarise(
        avg_total_sleep_hours = mean(total_sleep_time_hours, na.rm = TRUE),
        avg_stage_sleep_hours = mean(stage_sleep_time, na.rm = TRUE),
        sd_total_sleep_hours = sd(total_sleep_time_hours, na.rm = TRUE),
        sd_stage_sleep_hours = sd(stage_sleep_time, na.rm = TRUE),
        avg_bed_time_sd = sd(as.numeric(difftime(bed_time, as.Date(bed_time), units = "mins")), na.rm = TRUE),
        avg_wake_time_sd = sd(as.numeric(difftime(wake_time, as.Date(wake_time), units = "mins")), na.rm = TRUE),
        avg_awake_hours = mean(awake_duration, na.rm = TRUE),
        avg_deep_hours = mean(deep_duration, na.rm = TRUE),
        avg_rem_hours = mean(rem_duration, na.rm = TRUE),
        avg_light_hours = mean(light_duration, na.rm = TRUE),
        total_days_with_data = n(),
        .groups = "drop"
      )
    
    # Combine with daily data
    daily_sleep <- left_join(daily_sleep, avg_sleep_summary, by = "user_id")
    
    all_participants_data[[participant_id]] <- daily_sleep
    cat("Processed:", participant_id, "\n")
  }
}

# Combine all participants
final_sleep_data <- bind_rows(all_participants_data)

# Summary of corrected durations
summary(final_sleep_data$total_sleep_time_hours)
# Report overlapping count before removal
cat("Overlap removal complete.\nTotal participants:", length(unique(final_sleep_data$user_id)), "\n")
# Save the final combined data
#write.csv(final_sleep_data, "all_participants_sleep_data.csv")

final_sleep_data$id <- str_remove(sleep$user_id, "^AIREADI-")

df_sleep = merge(matched_ppts, final_sleep_data, by= 'id' )


# Total number of unique participants before filtering
total_participants <- length(unique(df_sleep$id))

# Total number of observations before filtering
n_before <- nrow(df_sleep)

# Store original for comparison
df_sleep_before <- df_sleep

# Apply filtering: remove sleep time outliers > 13h
df_sleep <- df_sleep[df_sleep$total_sleep_time_hours <= 13, ]
# Number of observations after filtering
n_after <- nrow(df_sleep)
# Number of observations removed
n_removed <- n_before - n_after

# Percent of observations removed
percent_removed <- round(100 * n_removed / n_before, 2)
# Number of unique participants with at least one removed observation
affected_participants <- length(unique(
  df_sleep_before$id[df_sleep_before$total_sleep_time_hours > 15]
))


cat("Number of observations removed:", n_removed, "\n")
cat("Percentage of total observations removed:", percent_removed, "%\n")
cat("Number of participants affected:", affected_participants, "out of", total_participants, "\n")

summary(df_sleep$total_sleep_time_hours)
summary(df_sleep$total_days_with_data)

##group differences 
df_sleep_p = df_sleep %>% 
  distinct(id, .keep_all = T)
chisq.test(df_sleep_p$glycemic_control, df_sleep_p$total_sleep_time_hours)
chisq.test(df_sleep_p$glycemic_control, df_sleep_p$total_days_with_data)

round(mean(df_sleep_p$total_days_with_data), digits=2)
round(sd(df_sleep_p$total_days_with_data), digits=2)
round(median(df_sleep_p$total_days_with_data), digits=2)
round(sd(df_sleep_p$total_days_with_data[df_sleep_p$glycemic_control == 'poor']), digits=2)
round(sd(df_sleep_p$total_days_with_data[df_sleep_p$glycemic_control == 'good']), digits=2)




###descriptive 1 plot


# # Step 2: Filter full data for those participants
plot_data_sleep <- df_sleep %>%
  semi_join(complete_ids, by = c("id", "glycemic_control")) %>%
  mutate(id = as.factor(id))

plot_data_sleep <- plot_data_sleep %>% mutate(id = as.factor(id))

# Just to check:
print(head(plot_data_sleep))
print(names(plot_data_sleep))

plot_data_sleep <- plot_data_sleep %>%
  group_by(id) %>%
  arrange(date) %>%
  mutate(day_number = as.integer(difftime(date, min(date), units = "days")) + 1) %>%
  ungroup()
# 

plot_data_10 <- plot_data_sleep %>%
  filter(day_number <= 10)

plot_data_sleep$date

group_means_sleep <- plot_data_sleep %>%
  group_by(glycemic_control, day_number) %>%
  summarise(mean_sleep = mean(total_sleep_time_hours), .groups = "drop")

group_means_10 <- group_means_sleep  %>%
  filter(day_number <= 10)


#### plot 1 main paper

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


okabe_ito <- c(
  "poor" = "#F0E442",
  "good" = "#00796B",
  "Between-group" = "#4D4D4D"
)

# Recalculate group means for the subset
group_means_10_subset <- plot_data_10_subset %>%
  group_by(glycemic_control, day_number) %>%
  summarise(mean_sleep = mean(total_sleep_time_hours), .groups = "drop")

# Now plot with the subsetted data

glycemic_labels <- c(
  "good" = "Controlled",
  "poor" = "Uncontrolled"
)


### main paper

ggplot(plot_data_10_subset, aes(x = day_number, y = total_sleep_time_hours, group = id, color = glycemic_control)) +
  geom_line(alpha = 0.3, linewidth = 0.7) +  # individual participant lines thicker
  geom_line(data = group_means_10_subset,
            aes(x = day_number, y = mean_sleep, color = glycemic_control, group = glycemic_control),
            linewidth = 1.8) +  # bold group mean line
  facet_wrap(~glycemic_control, ncol = 2, labeller = labeller(glycemic_control = glycemic_labels)) +
  scale_color_manual(values = okabe_ito[c("poor", "good")]) +
  labs(
    title = "Daily Step Counts Over First 10 Days by Glycemic Control Group (30 Random Participants per Group)",
    x = "Day",
    y = "Sleep duration",
    color = "Glycemic Control"
  ) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

## supplement

ggplot(plot_data_10, aes(x = day_number, y = total_sleep_time_hours, group = id, color = glycemic_control)) +
  geom_line(alpha = 0.3, linewidth = 0.7) +  # individual participant lines thicker
  geom_line(data = group_means_10_subset,
            aes(x = day_number, y = mean_sleep, color = glycemic_control, group = glycemic_control),
            linewidth = 1.8) +  # bold group mean line
  facet_wrap(~glycemic_control, ncol = 2, labeller = labeller(glycemic_control = glycemic_labels)) +
  scale_color_manual(values = okabe_ito[c("poor", "good")]) +
  labs(
    title = "Daily Step Counts Over First 10 Days by Glycemic Control Group (30 Random Participants per Group)",
    x = "Day",
    y = "Sleep duration",
    color = "Glycemic Control"
  ) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")



###

# ggplot(plot_data_10, aes(x = day_number, y = daily_mean_stress, group = id, color = glycemic_control)) +
#   geom_line(alpha = 0.3, linewidth = 0.7) +  # individual participant lines thicker
#   geom_line(data = group_means_10_subset,
#             aes(x = day_number, y = mean_stress, color = glycemic_control, group = glycemic_control),
#             linewidth = 1.8) +  # bold group mean line
#   facet_wrap(~glycemic_control, ncol = 2, labeller = labeller(glycemic_control = glycemic_labels)) +
#   scale_color_manual(values = okabe_ito[c("poor", "good")]) +
#   labs(
#     title = "Daily Step Counts Over First 10 Days by Glycemic Control Group (30 Random Participants per Group)",
#     x = "Day",
#     y = "Steps per Day",
#     color = "Glycemic Control"
#   ) +
#   theme_classic(base_size = 20) +
#   theme(legend.position = "none")
# 


#### plot 1 supplement 

ggplot(plot_data_10, aes(x = day_number, y = total_sleep_time_hours, group = id, color = glycemic_control)) +
  geom_line(alpha = 0.3) +  # individual participant lines
  geom_line(data = group_means_10,
            aes(x = day_number, y = mean_sleep, color = glycemic_control, group = glycemic_control),
            linewidth = 1.5) +  # bold group mean line
  facet_wrap(~glycemic_control, ncol = 2) +
  scale_color_manual(values = okabe_ito[c("poor", "good")]) +
  labs(
    title = "Daily Step Counts Over First 10 Days by Glycemic Control Group",
    x = "Day",
    y = "Steps per Day",
    color = "Glycemic Control"
  ) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")


### modeling 


# Fit model without temporal component (random intercept only)
model <- lmer(
  total_sleep_time_hours ~ glycemic_control + day_number + (1 | id),
  data = plot_data_sleep,
  REML = TRUE
)


# Extract variance components
vc <- as.data.frame(VarCorr(model))
var_id <- vc$vcov[vc$grp == "id"]        # Between-person variance
var_res <- vc$vcov[vc$grp == "Residual"] # Residual variance (within-person/day-to-day)
var_total <- var_id + var_res            # Total variance

# Intra-class correlation (ICC)
icc <- var_id / var_total
ratio_within_to_between <- var_res / var_id

# R² values
r2_vals <- r2(model)
r2_marginal <- r2_vals$R2_marginal       # Fixed effects only
r2_conditional <- r2_vals$R2_conditional # Fixed + random effects

# Approximate variance explained by fixed effect (glycemic_control)
var_between_group <- r2_marginal * var_total
var_within_group <- var_res
ratio_within_to_between_fixed <- var_within_group / var_between_group

# Print results
cat("Between-person variance (random intercept 'id'):", round(var_id, 2), "\n")
cat("Residual variance (within-person, day-to-day):", round(var_res, 2), "\n")
cat("Total variance:", round(var_total, 2), "\n")
cat("ICC (between-person variance / total variance):", round(icc, 3), "\n")
cat("Marginal R² (fixed effects):", round(r2_marginal, 3), "\n")
cat("Conditional R² (fixed + random effects):", round(r2_conditional, 3), "\n")
cat("Approximate variance explained by glycemic_control (fixed effect):", round(var_between_group, 2), "\n")
cat("Ratio of within-person (day-to-day) variance to fixed effect variance:", round(ratio_within_to_between_fixed, 2), "\n")



