## Define all the support functions to use in the main script

## Load libraries
library(readxl)
library(ggplot2)
library(SNFtool)
library(tidyverse)
library(survival)
library(survminer)
library(gplots)
library(forcats)
library(TCGAbiolinks)

## Load utils
### Colors
steel_blue  <- "steelblue"
steel_red   <- "#B22222"
steel_green <- "#99FF66"
pink        <- "pink"
gray        <- "gray"

## Custom Functions
## 01
## Load different path where data are located
relative_path <- function(subdir) {
  "Input: subdirectory with .txt files
   Output: path where files are located"
  ## Get the current data directory
  current_path <- getwd()
  ## Find the index of "Personal" in the path
  personal_index <- regexpr("Personal", current_path)
  ## Check if "Personal" is found in the path
  if (personal_index == -1) {
    stop("The 'Personal' directory was not found in the path.")
  }
  ## Extract the path up to "Personal" including "Personal"
  path_up_to_personal <- substr(current_path, 1, personal_index + attr(personal_index, "match.length") - 1)
  ## Specify the relative path
  relative_path <- "Data"
  ## Construct the new path using file.path()
  new_path <- gsub("/+", "/", file.path(path_up_to_personal, relative_path, subdir, "/"))
  return(new_path)
}

## 02
## Load data from multiple paths into different data frames with dynamic names
read_data <- function(prefix, folder_path) {
  "Input: prefix for folder data
   Output: dataframes combined into a named list"
  ## Construct the full file paths
  file_suffixes    <- c("Gene_Expression.txt", "Methy_Expression.txt", "Mirna_Expression.txt", "Survival.txt", "Clinical.csv")
  file_paths       <- paste0(folder_path, "/", prefix, file_suffixes)
  ## Read the files into a list of data frames
  data_frames <- lapply(file_paths, function(file) {
    ## Condition
    if (endsWith(file, "COLON_Methy_Expression.txt"))  {
      t(read.table(file, header = TRUE, row.names = NULL))
    } 
    else if (endsWith(file, "KIDNEY_Methy_Expression.txt"))  {
      t(read.delim(file, header = TRUE))
    } 
    else if (endsWith(file, "COLON_Methy_Expression.txt"))  {
      t(read.table(file, header = TRUE, row.names = NULL))
    } 
    else if (endsWith(file, "GLIO_Clinical.csv"))  {
      read.csv(file, header = TRUE)
    } 
    else{
      t(read.table(file, header = TRUE, fill = TRUE, fileEncoding = "UTF-8"))}
  })
  ## Combine the data frames into a named list
  names(data_frames) <- c("df_mRNA", "df_methy", "df_miRNA", "df_survival", "df_clinical")
  return(data_frames)
}

## 03
## Function to preprocess data
preprocess <- function(df) {
  "Input: dataframe
   Output: modified dataframe"
  row.names(df) <- sub("^(\\w+\\.\\w+\\.\\w+).*", "\\1", row.names(df))
  row.names(df) <- gsub("\\.", "-", row.names(df))
  return(df)
}

## 04
## Function to check the number of patients, number of genes and to see the data
check_df <- function(df, tumor_type) {
  "Input: dataframe, tumor_type
   Output: number of patients, and number of genes"
  ## Tumor Type
  cat("Tumor type:", tumor_type, "\n")
  ## Cycle for each data type
  for (data_type in c("mRNA", "methy", "miRNA")) {
    data_df <- df[[paste0("df_", data_type)]]
    cat("Patients", paste0(data_type, ":"), nrow(data_df), "\n")
    cat("Genes", paste0(data_type, ":"), ncol(data_df), "\n")
    cat("\n")
  }
}

## 05
## Define a function to generate data frames for each subtype
generate_subtype_df <- function(subtype_id, patient_ids, df_clinical) {
  "Input: subtype id, patient ids, clinical df
   Output: dataframe with merged informations"
  df <- df_clinical %>%
    filter(bcr_patient_barcode %in% patient_ids) %>%
    select(Patient_Code = bcr_patient_barcode,
           Age_At_Diagnosis = age_at_diagnosis,
           Age_At_Index = age_at_index,
           Gender = gender,
           Sync_Malignancy = synchronous_malignancy,
           Stage = ajcc_pathologic_stage,
           Origin_Point = tissue_or_organ_of_origin,
           Primary_Diagnosis = primary_diagnosis,
           Prior_Malignancy = prior_malignancy,
           State = state,
           Prior_Treatment = prior_treatment,
           Tumor_Classification = classification_of_tumor,
           Site_Biopsy = site_of_resection_or_biopsy,
           Tumor_Grade = tumor_grade,
           Race = race,
           Ethnicity = ethnicity,
           Vital_Status = vital_status,
           Radiation_Treatment = treatments_radiation_treatment_type,
           ) %>%
    mutate(Subtype = subtype_id)
  return(df)
}

## 06
## Define a function that computes kruskal test
run_kruskal_test <- function(variable_name, data){
  "Input: variable name, data
   Output: kruskal test"
  ## Create the formula dynamically
  formula <- as.formula(paste(variable_name, "~ Subtype"))
  ## Perform the Kruskal-Wallis test
  result <- kruskal.test(formula, data = data)
  ## Print the result
  print(result)
}
