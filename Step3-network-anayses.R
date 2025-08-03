
library(dplyr)
library(lubridate)
library(tidyr)
library(mlVAR)
library(qgraph)
library(ggplot2)


####  1. create merged df for network modeling 

##load dfs (ids, pa, stress, cgm)
matched_ppts <- read.csv("~/matched107_w_wearable_data_ppts.csv", header = TRUE)
pa <-  read.csv("~/physical_activity_clean-scireports.csv", header = TRUE)
stress <-  read.csv("~/stress-data_hourly-scirep.csv", header = TRUE) %>% 
  select(id, mean_stress, hour)
# Ensure stress hour is POSIXct
stress <- stress %>%
  mutate(hour = ymd_hms(hour))
cgm <-  read.csv("~/CGM_data.csv", header = TRUE)
cgm$id =cgm$folder_id
cgm= cgm %>% 
  select(body_cgm_value, id, time)


matched_ppts_ids <- read.csv("~/matched107_w_wearable_data_ppts.csv", header = TRUE) %>% 
  select(id)
cgm_matched <- cgm %>%
  filter(id %in% matched_ppts_ids$id)
# Convert time columns to POSIXct
cgm_matched <- cgm_matched %>%
  mutate(time = ymd_hms(time),
         hour = floor_date(time, unit = "hour"))
pa <- pa %>%
  mutate(hour = ymd_hms(hour))
# Merge on id and floored hour
merged_data <- cgm_matched %>%
  inner_join(pa, by = c("id", "hour"))
# Join stress to the existing merged_data
merged_all <- merged_data %>%
  left_join(stress, by = c("id", "hour"))

length(merged_all$id) 
length(unique(merged_all$id))

n_total <- nrow(merged_all)
n_unique <- nrow(distinct(merged_all))

cat("Total rows:", n_total, "\nDistinct rows:", n_unique, "\nDuplicates:", n_total - n_unique, "\n")

# Extract date from time
merged_all <- merged_all %>%
  mutate(date = as.Date(time))  # Use `time`, not `hour`, to preserve actual CGM timestamps
# Count number of unique days per participant
days_per_participant <- merged_all %>%
  group_by(id) %>%
  summarise(n_days = n_distinct(date)) %>%
  arrange(desc(n_days))

print(days_per_participant)

summary(days_per_participant$n_days)
hist(days_per_participant$n_days, main = "Distribution of Days per Participant", xlab = "Days")

merged_all$mean_stress
merged_all$total_steps
merged_all$body_cgm_value

merged_all <- merged_all %>%
  mutate(body_cgm_value = case_when(
    body_cgm_value == "High" ~ 401,  # or whatever your device upper limit is
    body_cgm_value == "Low" ~ 39,    # or your device lower limit
    TRUE ~ as.numeric(body_cgm_value)
  ))


 ####  2. network modeling 

# Make sure data is sorted by participant and time
merged_all <- merged_all %>%
  arrange(id, time) %>%
  select(id, time, glycemic_control_binary, mean_stress, total_steps, body_cgm_value) %>%
  group_by(id) %>%
  mutate(across(c(mean_stress, total_steps, body_cgm_value), scale)) %>%  # Optional: z-score per person
  ungroup()

data_list <- merged_all %>%
  select(id, glycemic_control_binary, time, mean_stress, total_steps, body_cgm_value) %>%
  group_by(id) %>%
  group_split() %>%
  lapply(function(df) df %>% select(mean_stress, total_steps, body_cgm_value) %>% as.matrix())

# Create group vector (0 = poor, 1 = good)
group_vector <- merged_all %>%
  select(id, glycemic_control_binary) %>%
  distinct() %>%
  pull(glycemic_control_binary)

merged_all$glycemic_control <- factor(merged_all$glycemic_control, levels = c("good", "poor"))

# Check unique values to confirm coding
unique(merged_all$glycemic_control)

# Recode numeric binary to factor with descriptive labels
merged_all$glycemic_control_f <- factor(
  ifelse(merged_all$glycemic_control_binary == 0, "good", "poor"),
  levels = c("good", "poor")
)

# Confirm
table(merged_all$glycemic_control_f)

# Variables for mlVAR
vars <- c("mean_stress", "total_steps", "body_cgm_value")

# Subset data by group
data_good <- subset(merged_all, glycemic_control_f == "good")
data_poor <- subset(merged_all, glycemic_control_f == "poor")

# Load mlVAR if not loaded
library(mlVAR)

data_good <- data_good %>%
  arrange(id, time) 
data_poor <- data_poor %>%
  arrange(id, time) 

vars <- c("mean_stress[,1]", "total_steps[,1]", "body_cgm_value[,1]")
colnames(data_good) <- gsub("\\[.*\\]", "", colnames(data_good))
vars <- c("mean_stress", "total_steps", "body_cgm_value")
colnames(data_poor) <- gsub("\\[.*\\]", "", colnames(data_poor))
vars <- c("mean_stress", "total_steps", "body_cgm_value")


# Fit mlVAR for good group
model_good <- mlVAR(
  data = data_good,
  vars = vars,
  idvar = "id",
  # lags = 12,
  lags = 12,
  temporal = "orthogonal",
  estimator = "lmer",
  verbose = TRUE
)


# Check how many subjects
num_subjects <- length(model_good$results$Beta$subject)
print(num_subjects)

# Get unique subject ids from your data in the same order (should match model ordering)
subject_ids <- unique(data_good$id)
print(subject_ids)

### main manuscript figure 2: controlled 

# Select 15 random partiicpants
set.seed(123)
plot_indices <- sample(seq_len(num_subjects), size = min(15, num_subjects))

# Set up plot layout for 2 rows x 5 columns (for 10 plots)

# Set layout: 2 rows x 2 columns (you had 3 rows x 2 columns comment but code uses 2x2)
par(mfrow = c(3, 5))

par(oma = c(2, 2, 4, 2))  # outer margins
par(mar = c(2, 2, 4, 2))  # inner margins

custom_labels <- c("Stress", "Steps", "Glucose")

for (i in plot_indices) {
  net_i <- model_good$results$Beta$subject[[i]][, , 1]  # lag 1 matrix
  
  # Remove self-loops by setting diagonal to zero
  diag(net_i) <- 0
  
  pid <- subject_ids[i]
  
  qgraph(
    net_i,
    layout = "circle",
    theme = "colorblind",
    labels = custom_labels,
    title = paste("", pid),
    edge.labels = TRUE,
    edge.label.cex = 3.9,
    maximum = max(abs(net_i)),
    minimum = 0.000000001,
    fade = FALSE,
    asize = 2.5,        # arrow size
    vsize = 33,         # node size
    label.cex = 1.45,   
    label.scale.equal= TRUE, # node label size
    margin = 0.4,       # adds more space inside plot
    loopAngle = pi/4,
    loopSize = 0.05     # smaller, tighter loops (appear more inside)
  )
}

# Reset plotting layout
par(mfrow = c(1, 1))

# Reset
par(mfrow = c(1, 1))

# Extract average network at lag 1
mean_net <- model_good$results$Beta$mean[, , 1]

custom_labels <- c("Stress", "Steps", "Glucose")

# Increase plot margins outside of qgraph
par(mar = c(6, 4, 3, 2))  # bottom, left, top, 
#right margins (default is c(5,4,4,2)+0.1)

# Remove self-loops by setting diagonal to zero
diag(mean_net) <- 0

qgraph(
  mean_net,
  layout = "circle",
  theme = "colorblind",
  labels = custom_labels,
  # title = paste("", pid),
  edge.labels = TRUE,
  edge.label.cex = 2.4,
  maximum = max(abs(net_i)),
  minimum = 0.000000001,
  fade = FALSE,
  asize = 2.5,        # arrow size
  vsize = 24,         # node size
  label.cex = 1.2,   
  label.scale.equal= TRUE, # node label size
  margin = 0.1,       # adds more space inside plot
  loopAngle = pi/4,
  loopSize = 0.05     # smaller, tighter loops (appear more inside)
)



### supplement network figure: controlled 

### plotting remaining individuals

# Total number of subjects
num_subjects <- length(model_good$results$Beta$subject)
subject_ids <- unique(data_good$id)

# Indices of already plotted subjects
set.seed(123)
plot_indices <- sample(seq_len(num_subjects), size = min(15, num_subjects))

# Get indices of remaining subjects
remaining_indices <- setdiff(seq_len(num_subjects), plot_indices)

# Optional: if no remaining subjects, stop
if (length(remaining_indices) == 0) {
  message("All subjects were already plotted.")
} else {
  # Set up plotting layout (adjust rows/columns to fit number of subjects)
  n <- length(remaining_indices)
  n_cols <- 5
  n_rows <- ceiling(n / n_cols)
  par(mfrow = c(n_rows, n_cols))
  par(oma = c(2, 2, 4, 2))  # outer margins
  par(mar = c(2, 2, 4, 2))  # inner margins
  
  custom_labels <- c("Stress", "Steps", "Glucose")
  
  for (i in remaining_indices) {
    net_i <- model_good$results$Beta$subject[[i]][, , 1]
    diag(net_i) <- 0
    pid <- subject_ids[i]
    
    qgraph(
      net_i,
      layout = "circle",
      theme = "colorblind",
      labels = custom_labels,
      title = paste("", pid),
      edge.labels = TRUE,
      edge.label.cex = 3.9,
      maximum = max(abs(net_i)),
      minimum = 0.000000001,
      fade = FALSE,
      asize = 2.5,
      vsize = 33,
      label.cex = 1.45,
      label.scale.equal = TRUE,
      margin = 0.4,
      loopAngle = pi/4,
      loopSize = 0.05
    )
  }
  
  # Reset plotting layout
  par(mfrow = c(1, 1))
}


######## plotting uncontrolled network 

# Fit mlVAR for good group
model_poor <- mlVAR(
  data = data_poor,
  vars = vars,
  idvar = "id",
  # lags = 12,
  lags = 12,
  temporal = "orthogonal",
  estimator = "lmer",
  verbose = TRUE
)



### main manuscript figure 2: controlled 


# Check how many subjects
num_subjects <- length(model_poor$results$Beta$subject)
print(num_subjects) #49 total

# Get unique subject ids from your data in the same order (should match model ordering)
subject_ids <- unique(data_bad$id)
print(subject_ids)

# Select 15 random indices 
set.seed(123)
plot_indices <- sample(seq_len(num_subjects), size = min(15, num_subjects))


par(mfrow = c(3, 5))

par(oma = c(2, 2, 4, 2))  # outer margins
par(mar = c(2, 2, 4, 2))  # inner margins

custom_labels <- c("Stress", "Steps", "Glucose")

for (i in plot_indices) {
  net_i <- model_poor$results$Beta$subject[[i]][, , 1]  # lag 1 matrix
  
  # Remove self-loops by setting diagonal to zero
  diag(net_i) <- 0
  
  pid <- subject_ids[i]
  
  qgraph(
    net_i,
    layout = "circle",
    theme = "colorblind",
    labels = custom_labels,
    title = paste("", pid),
    edge.labels = TRUE,
    edge.label.cex = 3.9,
    maximum = max(abs(net_i)),
    minimum = 0.000000001,
    fade = FALSE,
    asize = 2.5,        # arrow size
    vsize = 33,         # node size
    label.cex = 1.45,   
    label.scale.equal= TRUE, # node label size
    margin = 0.4,       # adds more space inside plot
    loopAngle = pi/4,
    loopSize = 0.05     # smaller, tighter loops (appear more inside)
  )
}

# Reset plotting layout
par(mfrow = c(1, 1))

# Reset
par(mfrow = c(1, 1))

##### average network: uncontrolled


# Extract average network at lag 1
mean_net_un <- model_poor$results$Beta$mean[, , 1]

custom_labels <- c("Stress", "Steps", "Glucose")

# Increase plot margins outside of qgraph
par(mar = c(6, 4, 3, 2))  # bottom, left, top, 
#right margins (default is c(5,4,4,2)+0.1)


# Remove self-loops by setting diagonal to zero
diag(mean_net) <- 0

qgraph(
  mean_net_un,
  layout = "circle",
  theme = "colorblind",
  labels = custom_labels,
  # title = paste("", pid),
  edge.labels = TRUE,
  edge.label.cex = 2.4,
  maximum = max(abs(mean_net_un)),
  minimum = 0.0000001,
  fade = FALSE,
  asize = 2.5,        # arrow size
  vsize = 24,         # node size
  label.cex = 1.2,   
  label.scale.equal= TRUE, # node label size
  margin = 0.1,       # adds more space inside plot
  loopAngle = pi/4,
  loopSize = 0.05     # smaller, tighter loops (appear more inside)
)


################# supplement figure network: uncontrolled (remaining ppts)

# Total number of subjects
num_subjects <- length(model_poor$results$Beta$subject)
subject_ids <- unique(data_poor$id)

# Indices of already plotted subjects
set.seed(123)
plot_indices <- sample(seq_len(num_subjects), size = min(15, num_subjects))

# Get indices of remaining subjects
remaining_indices <- setdiff(seq_len(num_subjects), plot_indices)

# Optional: if no remaining subjects, stop
if (length(remaining_indices) == 0) {
  message("All subjects were already plotted.")
} else {
  # Set up plotting layout (adjust rows/columns to fit number of subjects)
  n <- length(remaining_indices)
  n_cols <- 5
  n_rows <- ceiling(n / n_cols)
  par(mfrow = c(n_rows, n_cols))
  par(oma = c(2, 2, 4, 2))  # outer margins
  par(mar = c(2, 2, 4, 2))  # inner margins
  
  custom_labels <- c("Stress", "Steps", "Glucose")
  
  for (i in remaining_indices) {
    net_i <- model_poor$results$Beta$subject[[i]][, , 1]
    diag(net_i) <- 0
    pid <- subject_ids[i]
    
    qgraph(
      net_i,
      layout = "circle",
      theme = "colorblind",
      labels = custom_labels,
      title = paste("", pid),
      edge.labels = TRUE,
      edge.label.cex = 3.9,
      maximum = max(abs(net_i)),
      minimum = 0.000000001,
      fade = FALSE,
      asize = 2.5,
      vsize = 33,
      label.cex = 1.45,
      label.scale.equal = TRUE,
      margin = 0.4,
      loopAngle = pi/4,
      loopSize = 0.05
    )
  }
  
  # Reset plotting layout
  par(mfrow = c(1, 1))
}


################# testing difference in SD of random effects within vs. between

# --- Extract and compute SD of random effects for "good" group ---
random_matrices_good <- model_good$results$Beta$random  # List of 3x3 matrices
random_array_good <- simplify2array(random_matrices_good)  # Convert to 3x3xN array
sd_random_good <- apply(random_array_good, c(1, 2), sd, na.rm = TRUE)  # SD across participants
sd_vec_good <- as.vector(sd_random_good)  # Flatten to vector

# --- Extract and compute SD of random effects for "poor" group ---
random_matrices_poor <- model_poor$results$Beta$random
random_array_poor <- simplify2array(random_matrices_poor)
sd_random_poor <- apply(random_array_poor, c(1, 2), sd, na.rm = TRUE)
sd_vec_poor <- as.vector(sd_random_poor)

# --- Combine into one dataframe for comparison ---
sd_df <- data.frame(
  group = rep(c("good", "poor"), each = length(sd_vec_good)),
  sd = c(sd_vec_good, sd_vec_poor)
)

# --- Summary statistics ---
mean_good <- mean(sd_vec_good, na.rm = TRUE)
mean_poor <- mean(sd_vec_poor, na.rm = TRUE)

cat("Mean SD (good):", mean_good, "\n")
cat("Mean SD (poor):", mean_poor, "\n")

# --- Statistical comparison of heterogeneity ---
test_result <- wilcox.test(sd ~ group, data = sd_df)
print(test_result)

#custom_labels <- c("Stress \n score", "Step \n count", "Glucose \n levels")
custom_labels <- c("Stress", "Steps", "Glucose")


# Inspect summaries
summary(model_good)
summary(model_poor)

sd_good <- model_good$results$Beta$SD[, , 1]
sd_poor <- model_poor$results$Beta$SD[, , 1]

# Flatten SD matrices
vec_good <- as.vector(sd_good)
vec_poor <- as.vector(sd_poor)

# Combine into a data frame for comparison
sd_df <- data.frame(
  edge = rep(1:length(vec_good), 2),
  sd_value = c(vec_good, vec_poor),
  group = rep(c("good", "poor"), each = length(vec_good))
)


aggregate(sd_value ~ group, data = sd_df, mean)

ggplot(sd_df, aes(x = group, y = sd_value)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.7) +
  labs(title = "Random Effect SDs by Group", y = "SD of Random Effect", x = "Glycemic Control Group")


# === Extract SDs of random effects ===
sd_good <- model_good$results$beta$sd
sd_poor <- model_poor$results$beta$sd

model_good$results$Beta$sd

# Paired t-test (assumes normality of differences)
t.test(vec_poor, vec_good, paired = TRUE)

# OR non-parametric alternative
wilcox.test(vec_poor, vec_good, paired = TRUE)

mean(vec_good)  # Average SD (good group)
mean(vec_poor)  # Average SD (poor group)

vec_sd_good <- as.vector(sd_good)
vec_sd_poor <- as.vector(sd_poor)

fixed_good <- model_good$results$Beta$mean[, , 1]
fixed_poor <- model_poor$results$Beta$mean[, , 1]

vec_fixed_diff <- abs(as.vector(fixed_good - fixed_poor))


within_variability <- mean(c(vec_sd_good, vec_sd_poor))
between_variability <- mean(vec_fixed_diff)

cat("Average within-group SD:", within_variability, "\n")
cat("Average between-group difference in fixed effects:", between_variability, "\n")

#within-group heterogeneity appears ~2.5× larger than between-group differences. 

#### start here

# Step 1: Extract within-group SD matrices for GOOD and POOR at lag 1
sd_good <- model_good$results$Beta$SD[, , 1]
sd_poor <- model_poor$results$Beta$SD[, , 1]

# Step 2: Average the within-group SD matrices element-wise
sd_within_avg <- (sd_good + sd_poor) / 2

# Step 3: Extract mean networks for GOOD and POOR at lag 1
mean_good <- model_good$results$Beta$mean[, , 1]
mean_poor <- model_poor$results$Beta$mean[, , 1]

# Step 4: Calculate absolute between-group difference matrix
between_diff <- abs(mean_good - mean_poor)

# Step 5: Vectorize matrices to vectors
vec_sd_within_avg <- as.vector(sd_within_avg)
vec_between_diff <- as.vector(between_diff)

# Step 6: Run Wilcoxon signed-rank test (paired, one-sided)
wilcox_test <- wilcox.test(
  vec_sd_within_avg,
  vec_between_diff,
  paired = TRUE,
  alternative = "greater"
)

# Step 7: Print test result and effect size (ratio)
print(wilcox_test)
ratio <- mean(vec_sd_within_avg) / mean(vec_between_diff)
cat("Ratio of average within-group SD to between-group difference:", ratio, "\n")



#extract all individual random effects

library(tidyr)
library(dplyr)

# Assume:
# - vars: character vector of variable names, e.g. c("mean_stress", "total_steps", "body_cgm_value")
# - model_good$results$Beta$subject is a list of arrays [variables x variables x lags] or matrices [variables x variables]

# Number of variables
n_vars <- length(vars)

# Number of subjects
n_subjects_good <- length(model_good$results$Beta$subject)

# Initialize an empty list to store data frames
list_df <- vector("list", n_subjects_good)

for (i in seq_len(n_subjects_good)) {
  # Extract subject's network matrix at lag 1 (assuming it's 3D: var x var x lag)
  subj_net <- model_good$results$Beta$subject[[i]][, , 1]
  
  # Convert matrix to dataframe in long format
  df_subj <- as.data.frame(subj_net)
  colnames(df_subj) <- vars
  df_subj$from <- vars
  
  # Gather into long format: from, to, value
  df_long <- df_subj %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "random_effect") %>%
    mutate(
      subject_id = paste0("good_", i),
      lag = 1
    ) %>%
    select(subject_id, from, to, lag, random_effect)
  
  list_df[[i]] <- df_long
}

# Bind all rows into a single dataframe
df_good_random_effects <- bind_rows(list_df)

# View result
head(df_good_random_effects)



# For POOR group
n_subjects_poor <- length(model_poor$results$Beta$subject)

list_df_poor <- vector("list", n_subjects_poor)

for (i in seq_len(n_subjects_poor)) {
  subj_net <- model_poor$results$Beta$subject[[i]][, , 1]
  
  df_subj <- as.data.frame(subj_net)
  colnames(df_subj) <- vars
  df_subj$from <- vars
  
  df_long <- df_subj %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "random_effect") %>%
    mutate(
      subject_id = paste0("poor_", i),
      lag = 1
    ) %>%
    select(subject_id, from, to, lag, random_effect)
  
  list_df_poor[[i]] <- df_long
}

df_poor_random_effects <- bind_rows(list_df_poor)

head(df_poor_random_effects)

df_all_random_effects <- bind_rows(df_good_random_effects, df_poor_random_effects)

# length(unique(df_all_random_effects$subject_id))
# [1] 96


library(tidyr)
library(dplyr)

# Create a new column for edge names (from_to)
df_wide <- df_all_random_effects %>%
  mutate(edge = paste(from, to, sep = "_")) %>%
  select(subject_id, edge, random_effect) %>%
  pivot_wider(names_from = edge, values_from = random_effect)

# View the wide dataframe
head(df_wide)


##### right ids 

library(dplyr)
library(tidyr)

# Variables list (example)
# vars <- c("mean_stress", "total_steps", "body_cgm_value")

# --- GOOD group ---
n_subjects_good <- length(model_good$results$Beta$subject)

list_df_good <- vector("list", n_subjects_good)

for (i in seq_len(n_subjects_good)) {
  subj_net <- model_good$results$Beta$subject[[i]][, , 1]
  
  df_subj <- as.data.frame(subj_net)
  colnames(df_subj) <- vars
  df_subj$from <- vars
  
  df_long <- df_subj %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "random_effect") %>%
    mutate(
      subject_id = model_good$IDs[i],   # Use correct IDs here
      lag = 1
    ) %>%
    select(subject_id, from, to, lag, random_effect)
  
  list_df_good[[i]] <- df_long
}

df_good_random_effects <- bind_rows(list_df_good)

# --- POOR group ---
n_subjects_poor <- length(model_poor$results$Beta$subject)

list_df_poor <- vector("list", n_subjects_poor)

for (i in seq_len(n_subjects_poor)) {
  subj_net <- model_poor$results$Beta$subject[[i]][, , 1]
  
  df_subj <- as.data.frame(subj_net)
  colnames(df_subj) <- vars
  df_subj$from <- vars
  
  df_long <- df_subj %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "random_effect") %>%
    mutate(
      subject_id = model_poor$IDs[i],  # Use correct IDs here
      lag = 1
    ) %>%
    select(subject_id, from, to, lag, random_effect)
  
  list_df_poor[[i]] <- df_long
}

df_poor_random_effects <- bind_rows(list_df_poor)

# Combine GOOD and POOR dataframes
df_all_random_effects <- bind_rows(df_good_random_effects, df_poor_random_effects)

# Reshape to wide format: one row per subject, one column per edge (from_to)
df_wide <- df_all_random_effects %>%
  mutate(edge = paste(from, to, sep = "_")) %>%
  select(subject_id, edge, random_effect) %>%
  pivot_wider(names_from = edge, values_from = random_effect)

# View first few rows of wide dataframe
head(df_wide)
#write.csv(df_wide, "df_edgeweights_clustering.csv")