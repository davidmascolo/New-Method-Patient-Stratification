## Load all data types
rm(list = ls())
source("00_Support_Functions.R")
cat("Starting Data Import\n")
## Compute Path Data
path_kidney     <- relative_path("KRCCC")

## Load and preprocess data data frames from the specified paths
cat("KRCCC\n")
data_krccc    <- read_data("KIDNEY_", path_kidney)

## KRCCC 
df_mRNA_krccc     <- preprocess(data_krccc$df_mRNA)
df_methy_krccc    <- preprocess(data_krccc$df_methy)
df_miRNA_krccc    <- preprocess(data_krccc$df_miRNA)
df_survival_krccc <- t(data_krccc$df_survival)
df_clinical_krccc <- data_krccc$df_clinical[,-1]
## Other wrangling operations
colnames(df_methy_krccc) <- df_methy_krccc[1,]
df_methy_krccc <- df_methy_krccc[-1,]

## Check
check_df(data_krccc, "KRCCC")
cat("Data Import Done\n")