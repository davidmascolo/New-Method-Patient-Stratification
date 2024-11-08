## Load all data types
rm(list = ls())
source("00_Support_Functions.R")
cat("Starting Data Import\n")
## Compute Path Data
path_colon     <- relative_path("COLON")

## Load and preprocess data data frames from the specified paths
cat("COLON\n")
data_colon    <- read_data("COLON_", path_colon)

## COLON 
df_mRNA_colon     <- preprocess(data_colon$df_mRNA)
df_methy_colon    <- preprocess(data_colon$df_methy)
df_miRNA_colon    <- preprocess(data_colon$df_miRNA)
df_survival_colon <- t(data_colon$df_survival)
df_clinical_colon <- data_colon$df_clinical[,-1]
## Other wrangling operations
colnames(df_methy_colon) <- df_methy_colon[1,]
df_methy_colon <- df_methy_colon[-1,]

## Check
check_df(data_colon, "COLON")
cat("Data Import Done\n")