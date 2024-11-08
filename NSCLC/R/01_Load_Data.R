## Script useful to export data in .txt files that will be used for MATLAB

## Utils
source("00_Support_Functions.R")

## 5 types of data
## Rna Gene Counts
## Mutations
## Clinical Signatures
## Immune Signatures
## Curated Signatures


# Pre-Processing ----------------------------------------------------------
## Rna Gene Counts
df_rna_counts$Name        <- NULL
gene_names                <- df_rna_counts$Description
df_rna_counts$Description <- NULL
colnames(df_rna_counts)   <- clean_patient(colnames(df_rna_counts))
df_rna_counts             <- data.frame(t(df_rna_counts))
colnames(df_rna_counts)   <- gene_names
## Exome
wes_patients     <- df_wes$Tumor_Sample_Barcode
rownames(df_wes) <- clean_patient(wes_patients)
df_wes           <- df_wes[,-2]
df_wes$X         <- NULL
## Immune
immune_patients <- df_immune$Harmonized_SU2C_RNA_Tumor_Sample_ID_v2
df_immune$Harmonized_SU2C_RNA_Tumor_Sample_ID_v2 <- NULL
rownames(df_immune) <- immune_patients
## Myeloid
myeloid_patients <- df_myeloid$Harmonized_SU2C_RNA_Tumor_Sample_ID_v2
df_myeloid$Harmonized_SU2C_RNA_Tumor_Sample_ID_v2 <- NULL
rownames(df_myeloid) <- myeloid_patients
## Curated
curated_patients <- df_curated$Harmonized_SU2C_RNA_Tumor_Sample_ID_v2
df_curated$Harmonized_SU2C_RNA_Tumor_Sample_ID_v2 <- NULL
rownames(df_curated) <- curated_patients
## Clinical data
df_clinical$Harmonized_SU2C_Participant_ID_v2 <- paste(
  df_clinical$Harmonized_SU2C_Participant_ID_v2, "T1", sep = "-")

## Check
dim(df_rna_counts);dim(df_wes);dim(df_immune);dim(df_myeloid);dim(df_curated)
## Dimensions
## RNA Counts: 57523 genes
## Mutations : 16270 genes


# Filtering ---------------------------------------------------------------
## Filter the patients of WES data using the patients that are in RNA-Seq data
## Number of patients in common: 65
common_patients <- Reduce(intersect, list(rownames(df_rna_counts),
                                          rownames(df_wes),
                                          rownames(df_curated),
                                          rownames(df_myeloid),
                                          rownames(df_curated)))
## Filter the dataframe with two different sets of patients
df_rna_counts_final <- df_rna_counts[rownames(df_rna_counts) %in% common_patients,
                                     colnames(df_rna_counts) %in% df_10_hubs$Genes]
df_wes_final        <- df_wes[rownames(df_wes) %in% common_patients,
                              colnames(df_wes) %in% df_10_hubs$Genes]
df_immune_final     <- df_immune[rownames(df_immune) %in% common_patients, ]
df_myeloid_final    <- df_myeloid[rownames(df_myeloid) %in% common_patients,]
df_curated_final    <- df_curated[rownames(df_curated) %in% common_patients,]
df_clinical_final   <- df_clinical[df_clinical$Harmonized_SU2C_Participant_ID_v2 %in% common_patients,]
## Check
dim(df_rna_counts_final);dim(df_wes_final);dim(df_immune_final);dim(df_myeloid_final);dim(df_curated_final);dim(df_clinical_final)


# Export ------------------------------------------------------------------
## Export data
file_path_01 <- "/filtered_gene.txt"
file_path_02 <- "/filtered_mutation.txt"
file_path_03 <- "/filtered_immune.txt"
file_path_04 <- "/filtered_myeloid.txt"
file_path_05 <- "/filtered_curated.txt"
file_path_06 <- "/filtered_id_clinical.txt"
## Save the dataframe to a text file
write.table(df_rna_counts_final, file = paste0(data_path, "/", file_path_01),
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(df_wes_final, file = paste0(data_path, "/", file_path_02),
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(df_immune_final, file = paste0(data_path, "/", file_path_03),
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(df_myeloid_final, file = paste0(data_path, "/", file_path_04),
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(df_curated_final, file = paste0(data_path, "/", file_path_05),
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(df_clinical_final, file = paste0(data_path, "/", file_path_06),
            sep = "\t", row.names = FALSE, col.names = FALSE)