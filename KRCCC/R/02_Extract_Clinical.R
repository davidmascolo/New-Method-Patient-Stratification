## Load useful
rm(list = ls())
source("00_Support_Functions.R")
path_kidney     <- relative_path("KRCCC")

## Load raw id patients
id_patients <- data.frame(t(read.table(paste0(path_kidney,
                                              "id_Krccc_patients.txt"),
                                       sep = "")))
colnames(id_patients) <- "Patient_ID"
id_patients_clean <- sub("^([^-]+-[^-]+-[^-]+).*", "\\1",
                         gsub("\\.", "-", id_patients$Patient_ID))

## Download clincal data from TCGA
df_clinical <- GDCquery_clinic(project = "TCGA-KIRC",
                               type = "clinical")
## Check
dim(df_clinical)

## Format barcode patients
df_clinical$bcr_patient_barcode <- sub("^([^-]+-[^-]+-[^-]+).*", "\\1",
                                       gsub("\\.", "-", df_clinical$bcr_patient_barcode))

## Check intersection id patients
cat("N. subset id patients:", length(id_patients_clean))
cat("N. complete id patients", length(df_clinical$bcr_patient_barcode))
cat("N. common patients: ", length(intersect(id_patients_clean,
                                             df_clinical$bcr_patient_barcode)))
## Save output
write.csv(df_clinical, "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/KRCCC/Data/KIDNEY_Clinical.csv")
write.csv(id_patients_clean, "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/KRCCC/Data/Id_KIDNEY_patients.csv",
          row.names = F)
