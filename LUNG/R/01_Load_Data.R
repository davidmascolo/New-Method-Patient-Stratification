## Load all data types
rm(list = ls())
source("00_Support_Functions.R")
cat("Starting Data Import\n")
## Compute Path Data
path_lung      <- relative_path("LUNG")

## Load and preprocess data data frames from the specified paths
cat("LUNG\n")
data_lung    <- read_data("LUNG_", path_lung)

## LUNG 
df_mRNA_lung     <- preprocess(data_lung$df_mRNA)
df_methy_lung    <- preprocess(data_lung$df_methy)
df_miRNA_lung    <- preprocess(data_lung$df_miRNA)
df_survival_lung <- t(data_lung$df_survival)
df_clinical_lung <- data_lung$df_clinical[,-1]

## Check
check_df(data_lung, "LUNG")
cat("Data Import Done\n")