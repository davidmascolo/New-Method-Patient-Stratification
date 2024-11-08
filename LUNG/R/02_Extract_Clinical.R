## Load useful
rm(list = ls())

## Load raw id patients
id_patients <- data.frame(t(read.table(paste0(path_lung,
                                              "id_Lung_patients.txt"),
                                       sep = "")))
colnames(id_patients) <- "Patient_ID"
id_patients_clean <- sub("^([^-]+-[^-]+-[^-]+).*", "\\1",
                         gsub("\\.", "-", id_patients$Patient_ID))

## Download clincal data from TCGA
df_clinical <- GDCquery_clinic(project = "TCGA-LUSC",
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
write.csv(df_clinical, "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/LUNG/Data/LUNG_Clinical.csv")
write.csv(id_patients_clean, "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/LUNG/Data/Id_LUNG_patients.csv",
          row.names = F)