library(dplyr)
library(tidyr)
library(purrr)
library(forcats)

load_subject_performance_data <- function () {
  
  # Load and prepare all data
  var_data_path <- "../../Park-MOVE_fnirs_dataset_v2/IMU_data/imu_variability_parameters.csv"
  gait_path <- "../../Park-MOVE_fnirs_dataset_v2/IMU_data/imu_gait_parameters.csv"
  turns_file <- "../data/mobility_lab_turn_parameters.csv"
  
  time_file <- paste("../../Park-MOVE_fnirs_dataset_v2/Task_data/auditory_stroop_answer_time.csv", sep="")
  nav_file <- paste("../../Park-MOVE_fnirs_dataset_v2/Task_data/navigation_mistakes.csv", sep="")
  acc_file <- paste("../../Park-MOVE_fnirs_dataset_v2/Task_data/auditory_stroop_accuracy.csv", sep="")
  
  # Do we filter any subjects?
  filter_subjects <- c("NA")
  
  # Identifiers
  csv_path <- paste("../../Park-MOVE_fnirs_dataset_v1/identifiers_YA.csv", sep="")
  identifiers_ya <- read.csv(csv_path)
  csv_path <- paste("../../Park-MOVE_fnirs_dataset_v1/identifiers_OA.csv", sep="")
  identifiers_oa <- read.csv(csv_path)
  csv_path <- paste("../../Park-MOVE_fnirs_dataset_v1/identifiers_PD.csv", sep="")
  identifiers_pd <- read.csv(csv_path)
  
  # Assign group function
  assign_group <- function(df, identifiers_ya, identifiers_oa, identifiers_pd){
    df$group <- case_when(
      df$subject %in% identifiers_ya$id_nummer ~ "YA",
      df$subject %in% identifiers_oa$id_nummer ~ "OA",
      df$subject %in% identifiers_pd$id_nummer ~ "PD",
      TRUE ~ NA_character_
    )
    df$group <- factor(df$group, levels=c('YA', 'OA', 'PD'))
    return(df)
  }
  
  # Helper function to calculate DT cost
  calculate_dt_cost <- function(dt, st) {
    -(dt - st) / st * 100
  }

  # Process gait data
  gait_data <- read.csv(gait_path)
  
  # Gait variability data
  gait_variability_data <- read.csv(var_data_path)
  
  # Add averages of left and right feet
  gait_data <- gait_data %>%
    mutate(Cadence.LR = (Cadence.L + Cadence.R) / 2) %>%
    mutate(Single.Support.LR = (Single.Support.L + Single.Support.R) / 2) %>%
    mutate(Step.Count.LR = (Step.Count.L + Step.Count.R) / 2) %>%
    mutate(Step.Time.LR = (Step.Time.L + Step.Time.R) / 2) %>%
    mutate(Stride.Length.LR = (Stride.Length.L + Stride.Length.R) / 2) %>%
    mutate(Walking.Speed.LR = (Walking.Speed.L + Walking.Speed.R) / 2) %>%
    mutate(Stance.Time.LR = (Stance.Time.L + Stance.Time.R) / 2)
  
  # We're not interested in ST_stand
  gait_data <- gait_data %>% filter(trial_type != "Stand_still_and_Aud_Stroop")
  
  # Rename condition
  gait_data$trial_type <- gsub('Navigation', 'Navigated_walking', gait_data$trial_type)
  gait_data$trial_type <- gsub("^Straight_walking$", 'ST_walk', gait_data$trial_type)
  gait_data$trial_type <- gsub("^Straight_walking_and_Aud_Stroop$", 'DT_walk', gait_data$trial_type)
  gait_data$trial_type <- gsub("^Navigated_walking$", 'ST_navigation', gait_data$trial_type)
  gait_data$trial_type <- gsub("^Navigated_walking_and_Aud_Stroop$", 'DT_navigation', gait_data$trial_type)
  
  # Get averages of all blocks per subject
  gait_data_per_subj <- gait_data %>%
    group_by(subject, session, trial_type) %>%
    summarise(mean_single_support = mean(Single.Support.LR, na.rm = TRUE), 
              mean_cadence = mean(Cadence.LR, na.rm = TRUE), 
              mean_step_count = mean(Step.Count.LR, na.rm = TRUE),
              mean_step_time = mean(Step.Time.LR, na.rm = TRUE),
              mean_stride_length = mean(Stride.Length.LR, na.rm = TRUE),
              mean_walking_speed = mean(Walking.Speed.LR, na.rm = TRUE),
              mean_stance_time = mean(Stance.Time.LR, na.rm = TRUE),
              .groups='drop')
  
  # Make wider so we have one column per trial type
  gait_data_per_subj_wide <- gait_data_per_subj %>%
    pivot_wider(names_from = c(session, trial_type),
                values_from = c(mean_single_support,
                                mean_cadence,
                                mean_step_count,
                                mean_step_time,
                                mean_stride_length,
                                mean_walking_speed,
                                mean_stance_time))
  
  
  # Also do variability data
  gait_variability_data$trial_type <- gsub('Navigation', 'Navigated_walking', gait_variability_data$trial_type)
  gait_variability_data$trial_type <- gsub("^Straight_walking$", 'ST_walk', gait_variability_data$trial_type)
  gait_variability_data$trial_type <- gsub("^Straight_walking_and_Aud_Stroop$", 'DT_walk', gait_variability_data$trial_type)
  gait_variability_data$trial_type <- gsub("^Navigated_walking$", 'ST_navigation', gait_variability_data$trial_type)
  gait_variability_data$trial_type <- gsub("^Navigated_walking_and_Aud_Stroop$", 'DT_navigation', gait_variability_data$trial_type)
  
  # Get step time variability in ms instead of s
  gait_variability_data$Step.Time.Variability = gait_variability_data$Step.Time.Variability * 1000
  
  # Get averages of all blocks per subject
  gait_variability_data_per_subj <- gait_variability_data %>%
    group_by(subject, session, trial_type) %>%
    summarise(step_time_variability = mean(Step.Time.Variability, na.rm = TRUE),
              stride_length_variability = mean(Stride.Length.Variability, na.rm = TRUE), 
              step_time_asymmetry_percent = mean(Step.Time.Asymmetry.Percent, na.rm = TRUE),
              stride_length_asymmetry_percent = mean(Stride.Length.Asymmetry.Percent, na.rm = TRUE),
              .groups='drop')
  
  # Make wider so we have one column per trial type
  gait_variability_data_per_subj_wide <- gait_variability_data_per_subj %>%
    pivot_wider(names_from = c(session, trial_type),
                values_from = c(step_time_variability, 
                                stride_length_variability,
                                step_time_asymmetry_percent,
                                stride_length_asymmetry_percent))
  

  turn_data <- read.csv(turns_file)
  turn_data <- subset(turn_data, select = -c(turns_velocity_std, turns_steps_std,
                                             turns_duration_std, turns_angle_std, gait_cycles))
  turn_data_wide <- turn_data %>%
    pivot_wider(names_from = c(protocol),
                values_from = c(turns_velocity, 
                                turns_steps,
                                turns_duration,
                                turns_angle))

  # Merge
  all_gait_data <- merge(gait_data_per_subj_wide, gait_variability_data_per_subj_wide, by = "subject", all = TRUE)
  all_gait_data <- merge(all_gait_data, turn_data_wide, by = "subject", all = TRUE)
  
  # Add DT costs
  # Note: a lower number is better = add (-1) in front
  all_gait_data <- all_gait_data %>%
    mutate(diff_walk_speed_protocol1 = mean_walking_speed_protocol1_DT_walk - mean_walking_speed_protocol1_ST_walk,
           diff_walk_speed_protocol3 = mean_walking_speed_protocol3_DT_navigation - mean_walking_speed_protocol3_ST_navigation,
           dt_cost_walk_speed_protocol1 = calculate_dt_cost(mean_walking_speed_protocol1_DT_walk, mean_walking_speed_protocol1_ST_walk),
           dt_cost_walk_speed_protocol3 = calculate_dt_cost(mean_walking_speed_protocol3_DT_navigation, mean_walking_speed_protocol3_ST_navigation),
           dt_cost_stride_length_protocol1 = calculate_dt_cost(mean_stride_length_protocol1_DT_walk, mean_stride_length_protocol1_ST_walk),
           dt_cost_stride_length_protocol3 = calculate_dt_cost(mean_stride_length_protocol3_DT_navigation, mean_stride_length_protocol3_ST_navigation),
           dt_cost_cadence_protocol1 = calculate_dt_cost(mean_cadence_protocol1_DT_walk, mean_cadence_protocol1_ST_walk),
           dt_cost_cadence_protocol3 = calculate_dt_cost(mean_cadence_protocol3_DT_navigation, mean_cadence_protocol3_ST_navigation),
           dt_cost_single_support_protocol1 = calculate_dt_cost(mean_single_support_protocol1_DT_walk, mean_single_support_protocol1_ST_walk),
           dt_cost_single_support_protocol3 = calculate_dt_cost(mean_single_support_protocol3_DT_navigation, mean_single_support_protocol3_ST_navigation),
           dt_cost_step_time_variability_protocol1 = -calculate_dt_cost(step_time_variability_protocol1_DT_walk, step_time_variability_protocol1_ST_walk),
           dt_cost_step_time_variability_protocol3 = -calculate_dt_cost(step_time_variability_protocol3_DT_navigation, step_time_variability_protocol3_ST_navigation),)

  # Process auditory stroop data
  acc_data <- read.csv(acc_file)
  acc_data_long <- pivot_longer(acc_data, 
                                cols = starts_with("accuracy_"), 
                                names_to = "accuracy_variable", 
                                values_to = "accuracy_value")
  
  # Only take protocol 3
  acc_data <- acc_data[acc_data$protocol == 'protocol_3', ]
  
  # Make wider so we have one column per trial type
  acc_data_wide <- acc_data[c('subject', 'block_type', 'protocol', 'accuracy_total')]
  acc_data_wide <- acc_data_wide %>%
    pivot_wider(names_from = c(block_type, protocol),
                values_from = c(accuracy_total))
  
  # Stroop times
  time_data <- read.csv(time_file)
  
  # Only take protocol 3
  time_data <- time_data[time_data$protocol == 'protocol_3', ]
  
  # Get avg answer time per subject
  stroop_time_subject_means <- time_data %>%
    group_by(subject, protocol, block_type) %>%
    summarise(stroop_time_mean_value = mean(answer_time, na.rm = TRUE), .groups='drop')
  
  # Make wider so we have one column per trial type
  stroop_time_subject_means_wide <- stroop_time_subject_means %>%
    pivot_wider(names_from = c(block_type, protocol),
                values_from = c(stroop_time_mean_value))
  
  # Merge with acc
  as_data_wide <- merge(acc_data_wide, stroop_time_subject_means_wide, by = "subject", all = TRUE, suffixes = c("_stroop_acc", "_stroop_time"))
  
  # Navigation data
  nav_data <- read.csv(nav_file)
  names(nav_data)[names(nav_data) == 'id_nummer'] <- 'subject'
  nav_data <- nav_data %>% select(-nav_s2_case_interpret, -nav_s3_case_interpret, -comment, -navigation_score_complete)
  
  # Merge all performance data
  performance_data <- merge(all_gait_data, as_data_wide, by = "subject", all = TRUE)
  performance_data <- merge(performance_data, nav_data, by = "subject", all = TRUE)
  
  # Assign group
  performance_data <- assign_group(performance_data, identifiers_ya, identifiers_oa, identifiers_pd)
  
  # Filter out YA
  performance_data <- performance_data[!performance_data$group == 'YA', ]
  
  # Cleaner names
  labels(performance_data)  <- c(
    dt_cost_walk_speed_protocol1 = "DT cost walking speed, %",
    dt_cost_walk_speed_protocol3 = "DT cost walking speed, %",
    dt_cost_stride_length_protocol1 = "DT cost stride length, %",
    dt_cost_stride_length_protocol3 = "DT cost stride length, %",
    dt_cost_step_time_variability_protocol1 = "DT cost step time variability, %",
    dt_cost_step_time_variability_protocol3 = "DT cost step time variability, %",
    mean_walking_speed_protocol1_ST_walk = "Walking speed (ST), m/s",
    mean_walking_speed_protocol1_DT_walk = "Walking speed (DT), m/s",
    mean_walking_speed_protocol2_ST_walk = "Walking speed (straight), m/s",
    mean_walking_speed_protocol2_ST_navigation = "Walking speed (navigation), m/s",
    mean_walking_speed_protocol3_ST_navigation = "Walking speed (ST), m/s",
    mean_walking_speed_protocol3_DT_navigation = "Walking speed (DT), m/s",
    mean_stride_length_protocol1_ST_walk = "Stride length (ST), m",
    mean_stride_length_protocol1_DT_walk = "Stride length (DT), m",
    mean_stride_length_protocol2_ST_walk = "Stride length (straight), m",
    mean_stride_length_protocol2_ST_navigation = "Stride length (navigation), m",
    mean_stride_length_protocol3_ST_navigation = "Stride length (ST), m",
    mean_stride_length_protocol3_DT_navigation = "Stride length (DT), m",
    mean_step_time_protocol1_ST_walk = "Step time (ST), s",
    mean_step_time_protocol1_DT_walk = "Step time (DT), s",
    mean_step_time_protocol2_ST_walk = "Step time (straight), s",
    mean_step_time_protocol2_ST_navigation = "Step time (navigation), s",
    mean_step_time_protocol3_ST_navigation = "Step time (ST), s",
    mean_step_time_protocol3_DT_navigation = "Step time (DT), s",
    step_time_variability_protocol2_ST_walk = "Step time variability (straight), ms",
    step_time_variability_protocol1_ST_walk = "Step time variability (ST), ms",
    step_time_variability_protocol1_DT_walk = "Step time variability (DT), ms",
    step_time_variability_protocol2_ST_navigation = "Step time variability (navigation), ms",
    step_time_variability_protocol3_ST_navigation = "Step time variability (ST), ms",
    step_time_variability_protocol3_DT_navigation = "Step time variability (DT), ms",
    turns_velocity_protocol2 = "Turn velocity (peak), deg/s",
    turns_velocity_protocol3 = "Turn velocity (peak), deg/s",
    step_time_variability_protocol3_DT_navigation = "Step time variability (DT), ms",
    DT_protocol_3_stroop_acc = "Accuracy (DT), %",
    DT_protocol_3_stroop_time = "Answer time (DT), s")

  return(performance_data)
  
}
