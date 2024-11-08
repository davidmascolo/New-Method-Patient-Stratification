## Load useful
rm(list = ls())
source("00_Support_Functions.R")

## EDA considering clinical data
## Set data path
data_path <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/KRCCC/Data/"

## Load clinical variables
df_clinical <- read.csv(paste0(data_path, "KIDNEY_Clinical.csv"))
id_patients <- read.csv(paste0(data_path, "id_KIDNEY_patients.csv"))
id_patients <- id_patients$x

## Select only useful features
df_clinical_subset <- df_clinical %>% 
  select(synchronous_malignancy, ajcc_pathologic_stage, 
         days_to_diagnosis, last_known_disease_status,
         tissue_or_organ_of_origin, days_to_last_follow_up,
         age_at_diagnosis, primary_diagnosis, prior_malignancy,
         year_of_diagnosis, state, prior_treatment, ajcc_staging_system_edition,
         ajcc_pathologic_t, ajcc_pathologic_n, ajcc_pathologic_m, 
         classification_of_tumor, site_of_resection_or_biopsy,
         tumor_grade, progression_or_recurrence, alcohol_history,
         race, gender, ethnicity, vital_status, age_at_index,
         days_to_birth, year_of_birth, treatments_radiation_treatment_type,
         treatments_radiation_treatment_or_therapy, bcr_patient_barcode) %>% 
  filter(bcr_patient_barcode %in% id_patients)
## Check
dim(df_clinical_subset)

## Formatting
df_clinical_subset$ajcc_pathologic_stage <-
  fct_infreq(df_clinical_subset$ajcc_pathologic_stage)
df_clinical_subset$ajcc_pathologic_stage <- factor(df_clinical_subset$ajcc_pathologic_stage,
                                                   levels = c("Stage I", "Stage IA", "Stage IIA", "Stage IIIA",
                                                              "Stage IB", "Stage IIB", "Stage IIIC", "Stage IV", "Stage X"),
                                                   ordered = TRUE)

## EDA
## Gender
(p1 <- ggplot(df_clinical_subset, aes(x = gender, fill = gender)) +
    geom_bar(stat = "count") +
    geom_text(aes(label = paste0(after_stat(count), ", ", "(",
                                 round(after_stat(count)/sum(after_stat(count)) * 100, 2), "%)")),
              stat = "count", vjust = -0.3, hjust = 0.3, size = 4) +
    scale_fill_manual(values = c("female" = "pink",
                                 "male" = "steelblue"),
                      name = "Gender",
                      labels = c("female" = "F", "male" = "M")) +
    labs(title = "Gender Distribution", x = "", y = "Count") +
    theme(axis.text.x = element_blank()))

## Race
(p2 <- ggplot(df_clinical_subset, aes(x = race, fill = race)) +
    geom_bar(stat = "count") +
    geom_text(aes(label = paste0(after_stat(count), ", ", "(",
                                 round(after_stat(count)/sum(after_stat(count)) * 100, 2), "%)")),
              stat = "count", vjust = -0.3, hjust = 0.2, size = 4) +
    scale_fill_manual(values = c("asian" = steel_red,
                                 "black or african american" = steel_blue,
                                 "not reported" = gray,
                                 "white" = pink),
                      name = "Race") +
    labs(title = "Race Distribution", x = "", y = "Count") +
    theme(axis.text.x = element_blank()))

## Vital Status
(p3 <- ggplot(df_clinical_subset, aes(x = vital_status, fill = vital_status)) +
    geom_bar(stat = "count") +
    geom_text(aes(label = paste0(after_stat(count), ", ", "(",
                                 round(after_stat(count)/sum(after_stat(count)) * 100, 2), "%)")),
              stat = "count", vjust = -0.3, hjust = 0.2, size = 4) +
    scale_fill_manual(values = c("Alive" = steel_red,
                                 "Dead" = steel_blue),
                      name = "Vital Status") +
    labs(title = "Vital Status Distribution", x = "", y = "Value") +
    theme(axis.text.x = element_blank()))

## Age Distribution
(p4 <- ggplot(df_clinical_subset, aes(x = gender, y = age_at_index,
                                      fill = gender)) +
    geom_boxplot() +
    scale_fill_manual(values = c("female" = pink,
                                 "male" = steel_blue),
                      name = "Gender") +
    geom_jitter(size = 1, alpha = 5) +
    labs(title = "Age Distribution",
         x = "", y = "Age") +
    theme(axis.text.x = element_blank()))

## Radiation Treatment
(p5 <- ggplot(df_clinical_subset,
              aes(x = treatments_radiation_treatment_or_therapy,
                  fill = treatments_radiation_treatment_or_therapy)) +
    geom_bar(stat = "count") +
    geom_text(aes(label = paste0(after_stat(count), ", ", "(",
                                 round(after_stat(count)/sum(after_stat(count)) * 100, 2), "%)")),
              stat = "count", vjust = -0.3, hjust = 0.2, size = 4) +
    scale_fill_manual(values = c("no" = steel_blue,
                                 "not reported" = gray,
                                 "yes" = steel_red),
                      name = "Radiation Treatment") +
    labs(title = "Radiation Treatment Distribution",
         x = "", y = "Value") +
    theme(axis.text.x = element_blank()))

## Tumor Stage
(p6 <- ggplot(df_clinical_subset,
              aes(x = factor(ajcc_pathologic_stage, levels = c("Stage I", "Stage IA", "Stage IIA", "Stage IIIA", "Stage IB", "Stage IIB", "Stage IIIC", "Stage IV", "Stage X"), ordered = TRUE),
                  fill = ajcc_pathologic_stage)) +
    geom_bar(stat = "count") +
    geom_text(aes(label = paste0(after_stat(count), ", ", "(",
                                 round(after_stat(count)/sum(after_stat(count)) * 100, 2), "%)")),
              stat = "count", vjust = -0.3, hjust = 0.3, size = 4) +
    scale_fill_manual(values = c("Stage I"    = "#C6DBEF",
                                 "Stage IV"   = "#08519C"),
                      name = "Stage") +
    labs(title = "Tumor Stage Distribution",
         x = "", y = "Count"))

## Tissue or Organ of Origin
(p7 <- ggplot(df_clinical_subset,
              aes(x = factor(tissue_or_organ_of_origin),
                  fill = tissue_or_organ_of_origin)) +
    geom_bar(stat = "count") +
    geom_text(aes(label = paste0(round(after_stat(count)/sum(after_stat(count)) * 100, 2), "%")),
              stat = "count", vjust = -0.3, hjust = 0.3, size = 4) +
    scale_fill_manual(values = c("Kidney, NOS"      = "#08519C"),
                      name = "Origin Point") +
    labs(title = "Tissue Organ Origin Point Distribution",
         x = "", y = "Count"))

## Primary Diagnosis
(p8 <- ggplot(df_clinical_subset,
              aes(x = factor(primary_diagnosis),
                  fill = primary_diagnosis)) +
    geom_bar(stat = "count") +
    geom_text(aes(label = paste0(round(after_stat(count)/sum(after_stat(count)) * 100, 2), "%")),
              stat = "count", vjust = -0.3, hjust = 0.1, size = 4) +
    scale_fill_manual(values = c("Clear cell adenocarcinoma, NOS"= "#4292C6"),
                      name = "Diagnosis") +
    labs(title = "Primary Diagnosis Distribution",
         x = "", y = "Count") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)))

## Prior Malignacy
(p9 <- ggplot(df_clinical_subset,
              aes(x = factor(prior_malignancy),
                  fill = prior_malignancy)) +
    geom_bar(stat = "count") +
    geom_text(aes(label = paste0(after_stat(count), ", ", "(",
                                 round(after_stat(count)/sum(after_stat(count)) * 100, 2), "%)")),
              stat = "count", vjust = -0.3, hjust = 0.2, size = 4)  +
    scale_fill_manual(values = c("no" ="#2171B5",
                                 "yes"= "#08306B"),
                      name = "Prior Malignacy") +
    labs(title = "Prior Malignacy Distribution",
         x = "", y = "Count") +
    theme(axis.text.x = element_blank()))

## Median Age
(median_age <- group_by(df_clinical_subset, gender) %>% 
    summarise(count = n(),
              mean = mean(age_at_index, na.rm = T),
              sd = sd(age_at_index, na.rm = T),
              median = median(age_at_index, na.rm = T),
              IQR = IQR(age_at_index, na.rm = T)))

## Print
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
print(median_age)