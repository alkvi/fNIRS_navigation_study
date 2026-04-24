library(dplyr)
library(tidyr)
library(purrr)
library(forcats)

load_subject_data <- function () {
  
  # Load and prepare all data
  demo_file <- paste("../../Park-MOVE_fnirs_dataset_v2/basic_demographics.csv", sep="")
  redcap_file <- paste("../../Park-MOVE_fnirs_dataset_v2/REDcap_data/All_REDcap_data.csv", sep="")

  # Do we filter any subjects?
  filter_subjects <- c("NA")
  
  # Identifiers
  csv_path <- paste("../../Park-MOVE_fnirs_dataset_v2/identifiers_YA.csv", sep="")
  identifiers_ya <- read.csv(csv_path)
  csv_path <- paste("../../Park-MOVE_fnirs_dataset_v2/identifiers_OA.csv", sep="")
  identifiers_oa <- read.csv(csv_path)
  csv_path <- paste("../../Park-MOVE_fnirs_dataset_v2/identifiers_PD.csv", sep="")
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
  
  # Demographic data
  demo_data <- read.csv(demo_file)
  demo_data['sex'][demo_data['sex'] == 0] <- 'Male'
  demo_data['sex'][demo_data['sex'] == 1] <- 'Female'
  demo_data$sex <- factor(demo_data$sex, levels=c('Male', 'Female'))
  
  # REDcap data
  redcap_data <- read.csv(redcap_file) 
  names(redcap_data)[names(redcap_data) == 'id_nummer'] <- 'subject'
  redcap_data['ramlat_12_man'][redcap_data['ramlat_12_man'] == 0] <- 'No'
  redcap_data['ramlat_12_man'][redcap_data['ramlat_12_man'] == 1] <- 'Yes'
  redcap_data$ramlat_12_man <- factor(redcap_data$ramlat_12_man, levels=c('No', 'Yes'))
  
  # Get TMT contrast
  redcap_data <- redcap_data %>%
    mutate(tmt_4_tmt_2_contrast = tmt_4 - tmt_2)
  
  # Merge
  all_subject_data <- merge(demo_data, redcap_data, by = "subject", all = TRUE)
  all_subject_data <- assign_group(all_subject_data, identifiers_ya, identifiers_oa, identifiers_pd)
  
  # Cleaner names
  labels(all_subject_data)  <- c(
    age = "Age, yrs", 
    sex = "Gender, female",
    crf_utbildning_ar = "Education, yrs",
    weight = "Weight, kg",
    height = "Height, cm",
    frandin_grimby.x = "Frändin-Grimby",
    ramlat_12_man = "Falls in last 12 months, yes",
    mb_total = "Mini-BESTest score",
    g12_sum = "Walk-12 sum",
    tmt_2 = "TMT2, s",
    tmt_4 = "TMT4, s",
    tmt_4_tmt_2_contrast = "TMT4 - TMT2, s",
    cwit_3 = "CWIT3, s",
    ravlt_ret = "RAVLT retention, words",
    vf_sum = "Verbal fluency, total words")

  return(all_subject_data)
  
}
