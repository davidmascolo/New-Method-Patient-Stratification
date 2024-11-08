## This code is useful to analyze the clusters' significance for the
## partitions obtained using MATLAB
rm(list = ls())
source("00_Support_Functions.R")
cat("Starting Significance Generalized Louvain\n")

## Path
data_path <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/LUNG/Data/"
res_path  <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/LUNG/Res/"

## Load data
LUNG_Survival <- read.delim(paste0(data_path, "LUNG_Survival.txt"), header = T)
LUNG_Survival$PatientIDModified <- gsub("TCGA-([0-9]{2})-([0-9]{4}).*",
                                          "\\1.\\2", LUNG_Survival$PatientID)
id <- LUNG_Survival
## Extract data
## Load MATLAB partitions (for GenLouvain)
labels <- read.table(paste0(res_path, "/", "L.txt"), sep = ",")

## K-means
#set.seed()
k <- 4
kmeans_res <- kmeans(t(labels), k, algorithm = "Lloyd")
clusters   <- kmeans_res$cluster
(t <- table(clusters))

## Extract survival profiles
type_01 <- which(clusters == 1)
type_02 <- which(clusters == 2)
type_03 <- which(clusters == 3)
type_04 <- which(clusters == 4)

## Extract survival IDs and death IDs
surv_01 <- id$Survival[type_01]
surv_02 <- id$Survival[type_02]
surv_03 <- id$Survival[type_03]
surv_04 <- id$Survival[type_04]

death_01 <- id$Death[type_01]
death_02 <- id$Death[type_02]
death_03 <- id$Death[type_03]
death_04 <- id$Death[type_04]
## Collect dataframe
d1 <- data.frame(surv_01, death_01); colnames(d1) <- c("surv","death")
d2 <- data.frame(surv_02, death_02); colnames(d2) <- c("surv","death")
d3 <- data.frame(surv_03, death_03); colnames(d3) <- c("surv","death")
d4 <- data.frame(surv_04, death_04); colnames(d4) <- c("surv", "death")

## Divide datafram for Surv and Death
d1_ord <- data.frame(d1$surv[order(d1$surv)]/30, d1$death[order(d1$surv)]); colnames(d1_ord) <- c("surv","death")
d2_ord <- data.frame(d2$surv[order(d2$surv)]/30, d2$death[order(d2$surv)]); colnames(d2_ord) <- c("surv","death")
d3_ord <- data.frame(d3$surv[order(d3$surv)]/30, d3$death[order(d3$surv)]); colnames(d3_ord) <- c("surv","death")
d4_ord <- data.frame(d4$surv[order(d4$surv)]/30, d4$death[order(d4$surv)]); colnames(d4_ord) <- c("surv","death")

## Collect
d <- cbind(rbind(d1_ord, d2_ord, d3_ord, d4_ord),
           c(rep(1,t[1]), rep(2,t[2]),
             c(rep(3,t[3]), c(rep(4,t[4])))));colnames(d)[3] = "Subtype"

## Survival Analysis
fit_survival <- survfit(Surv(surv, death) ~ Subtype , data = d)

## Visualization with GG Plot
(p1 <- ggsurvplot(fit_survival, data = d, xlim = c(0, 120),
                  ggtheme = theme_bw(), palette = c("green", "blue", "red", "yellow"),
                  legend.title = "Subtype",
                  legend.labs = c("01", "02", "03", "04"),
                  xlab = "Survival Time (Months)",
                  ylab = "Survival Probability",
                  title = "LUNG Survival Analysis",
                  subtitle = "3 layers: mRNA Expression, DNA Methylation, miRNA",
                  legend = "right"))
## Significance
print(surv_pvalue(fit_survival)[,2])

## ***********************************************************
## Extract clinical fetures
df_clinical <- read.delim(paste0(data_path, "/",
                                 "LUNG_Clinical.csv"),
                          sep = ",", header = T)

## Extract for each subtype the patient ID
patient_id_type_01 <- id$PatientID[type_01]
patient_id_type_02 <- id$PatientID[type_02]
patient_id_type_03 <- id$PatientID[type_03]
patient_id_type_04 <- id$PatientID[type_04]

## Extract only common string
patient_id_type_01 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_01)
patient_id_type_02 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_02)
patient_id_type_03 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_03) 
patient_id_type_04 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_04)

## Define the list of subtype IDs and corresponding patient IDs
subtype_patient_ids <- list(
  "01" = patient_id_type_01,
  "02" = patient_id_type_02,
  "03" = patient_id_type_03,
  "04" = patient_id_type_04)

## Generate data frames for each subtype
subtype_dfs <- lapply(names(subtype_patient_ids), function(subtype_id) {
  generate_subtype_df(subtype_id,
                      subtype_patient_ids[[subtype_id]],
                      df_clinical)
})

## Median age
print(median(subtype_dfs[[1]]$Age_At_Index))
print(median(subtype_dfs[[2]]$Age_At_Index))
print(median(subtype_dfs[[3]]$Age_At_Index))
print(median(subtype_dfs[[4]]$Age_At_Index))

## Combine all data frames into a single data frame
df_patients <- bind_rows(subtype_dfs)
## Kruksal test
run_kruskal_test("Gender", df_patients) ## OK
run_kruskal_test("Age_At_Index", df_patients)
run_kruskal_test("Sync_Malignancy", df_patients)
run_kruskal_test("Stage", df_patients)
run_kruskal_test("Origin_Point", df_patients)
run_kruskal_test("Primary_Diagnosis", df_patients)
run_kruskal_test("Prior_Treatment", df_patients) ## OK
run_kruskal_test("Site_Biopsy", df_patients)
run_kruskal_test("Race", df_patients) ## OK
run_kruskal_test("Ethnicity", df_patients) ## OK
run_kruskal_test("Vital_Status", df_patients)
run_kruskal_test("Ajcc_Staging_System_Edition", df_patients)
run_kruskal_test("Ajcc_Pathologic_T", df_patients) ## OK
run_kruskal_test("Ajcc_pathologic_N", df_patients)
run_kruskal_test("Ajcc_Pathologic_M", df_patients)

## Analytics
## Gender
(p1 <- ggplot(df_patients, aes(x = Subtype,
                               fill = Gender)) +
    geom_bar(sta = "count", position = "dodge") +
    scale_fill_manual(values = c("female" = "#F8766D",
                                 "male" = "#D39200"),
                      name = "Gender") +
    labs(title = "Gender Distribution",
         x = "Subtype", y = "Count"))
## Percentages
df_patients %>%
  group_by(Subtype, Gender) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>%
  arrange(Subtype, desc(count))


## Prior Treatment Distribution
(p2 <- ggplot(df_patients, aes(x = Subtype,
                               fill = Prior_Treatment)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(values = c("No" ="#F8766D",
                                 "Yes"= "#D39200"),
                      name = "Prior Treatment") +
    labs(title = "Prior Treatment Distribution",
         x = "Subtype", y = "Count"))
## Percentages
df_patients %>%
  group_by(Subtype, Prior_Treatment) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>%
  arrange(Subtype, desc(count))

## Race Distribution
(p3 <- ggplot(df_patients, aes(x = Subtype,
                               fill = Race)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(values = c("asian" ="#F8766D",
                                 "black or african american"= "#D39200",
                                 "not reported" = "#1f77b4",
                                 "white" = "#2ca02c"),
                      name = "Race") +
    labs(title = "Race Distribution",
         x = "Subtype", y = "Count"))
## Percentages
df_patients %>%
  group_by(Subtype, Race) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>%
  arrange(Subtype, desc(count))

## Ethnicity Distribution
(p4 <- ggplot(df_patients, aes(x = Subtype,
                               fill = Ethnicity)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(values = c("hispanic or latino"  = "#1f77b4",
                                 "not hispanic or latino" = "#ff7f0e",
                                 "not reported" = "#2ca02c"),
                      name = "Ethnicity") +
    labs(title = "Ethnicity Distribution",
         x = "Subtype", y = "Count"))
## Percentages
df_patients %>%
  group_by(Subtype, Ethnicity) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>%
  arrange(Subtype, desc(count))

## Ajcc Pathologic M Distribution
(p5 <- ggplot(df_patients, aes(x = Subtype,
                               fill = fct_explicit_na(Ajcc_Pathologic_M, na_level = "NA"))) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(values = c("M0" = "#1f77b4",
                                 "M1" = "#ff7f0e",
                                 "NA" = "gray"),
                      name = "Ajcc Pathologic M Distribution") +
    labs(title = "Ajcc Pathologic M Distribution",
         x = "Subtype", y = "Count"))
## Percentages
df_patients %>% 
  group_by(Subtype, Ajcc_Pathologic_M) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>% 
  arrange(Subtype, desc(count))

cat("End Significance Generalized Louvain\n")