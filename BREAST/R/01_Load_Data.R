## Load all data types
rm(list = ls())
source("00_Support_Functions.R")
cat("Starting Data Import\n")
## Compute Path Data
path_breast      <- relative_path("BREAST")

## Load and preprocess data data frames from the specified paths
cat("BREAST\n")
data_breast    <- read_data("BREAST_", path_breast)

## BREAST 
df_mRNA_breast     <- preprocess(data_breast$df_mRNA)
df_methy_breast    <- preprocess(data_breast$df_methy)
df_miRNA_breast    <- preprocess(data_breast$df_miRNA)
df_survival_breast <- t(data_breast$df_survival)
df_clinical_breast <- data_breast$df_clinical[,-1]

## Check
check_df(data_breast, "BREAST")
cat("Data Import Done\n")