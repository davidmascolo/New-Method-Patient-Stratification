## Load all data types
source("00_Support_Functions.R")

cat("Starting Data Import\n")
## Compute Path Data
path_gbm      <- relative_path("GBM")

## Load and preprocess data data frames from the specified paths
cat("GBM\n")
data_gbm    <- read_data("GLIO_", path_gbm)

## GBM 
df_mRNA_gbm     <- preprocess(data_gbm$df_mRNA)
df_methy_gbm    <- preprocess(data_gbm$df_methy)
df_miRNA_gbm    <- preprocess(data_gbm$df_miRNA)
df_survival_gbm <- t(data_gbm$df_survival)
df_clinical_gbm <- data_gbm$df_clinical[,-1]

## Check
check_df(data_gbm, "GBM")
cat("Data Import Done\n")