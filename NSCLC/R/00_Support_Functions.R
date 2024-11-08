## Script useful to export data in .txt files that will be used for MATLAB

rm(list = ls())

## Utils
library(readxl)
library(survival)
library(survminer)
library(ggsurvfit)
library(dplyr)
library(aricode)

## Support function
## 01
## Create a function that returns the string in correct form
clean_patient <- function(patients_list){
  ## Input: list of patients code to format
  ## Output: code formatted
  return(gsub("\\.", "-", patients_list))
}

## 02
## Define a function to generate data frames for each subtype
generate_subtype_df <- function(subtype_id, patient_ids, df_clinical) {
  "Input: subtype id, patient ids, clinical df
   Output: dataframe with merged informations"
  df <- df_clinical %>%
    filter(Harmonized_SU2C_Participant_ID_v2 %in% patient_ids) %>%
    select(Patient_Code = Harmonized_SU2C_Participant_ID_v2,
           Age_At_Diagnosis = Patient_Age_at_Diagnosis,
           Gender = Patient_Sex,
           Smoking_Status = Patient_Smoking_Status,
           Stage = Initial_Stage,
           Histology = Histology_Harmonized,
           Line_Theraphy = Line_of_Therapy,
           Prior_Platinum = Prior_Platinum,
           Prior_TKI = Prior_TKI,
           Response = Response) %>%
    mutate(Subtype = subtype_id,
           Age_At_Diagnosis = as.numeric(Age_At_Diagnosis))
  return(df)
}

## 5 types of data
## Rna Gene Counts
## Mutations
## Clinical Signatures
## Immune Signatures
## Curated Signatures

## Path
data_path <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/NSCLC/Data/"

message("Starting Data Import")
message("--------------------")
## Load Data
df_rna_counts <- read.csv(paste0(data_path, "rna_counts.csv"))
df_wes        <- read.csv(paste0(data_path, "Exome_Genes.csv"))
df_immune     <- read_excel(paste0(data_path, "SU2C-MARK_Supplementary_Tables_Combined_v5_Filtered.xlsx"),
                            sheet = "Table_S18_Immune_Signatures")
df_myeloid    <- read_excel(paste0(data_path, "SU2C-MARK_Supplementary_Tables_Combined_v5_Filtered.xlsx"),
                            sheet = "Table_S19_Myeloid_Signatures")
df_curated    <- read_excel(paste0(data_path, "SU2C-MARK_Supplementary_Tables_Combined_v5_Filtered.xlsx"),
                            sheet = "Table_S20_Curated_Signatures")
df_clinical   <- read_excel(paste0(data_path, "SU2C-MARK_Supplementary_Tables_Combined_v5_Filtered.xlsx"),
                            sheet = "Table_S1_Clinical_Annotations",
                            skip = 2)
## Load 10% Differential Co-Expressed Genes
df_10_hubs <- read.delim(paste0(data_path, "10_genes.txt"),
                         header = T)
message("Data Import Done")