## This code is useful to analyze the clusters' significance for the
## partitions obtained using MATLAB
rm(list = ls())
source("00_Support_Functions.R")
cat("Starting Significance Generalized Louvain\n")

## Path
data_path <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/KRCCC/Data/"
res_path  <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/KRCCC/Res/"

## Load data
KIDNEY_Survival <- read.delim(paste0(data_path, "KIDNEY_Survival.txt"), header = T)
KIDNEY_Survival$PatientIDModified <- gsub("TCGA-([0-9]{2})-([0-9]{4}).*",
                                         "\\1.\\2", KIDNEY_Survival$PatientID)
id <- KIDNEY_Survival
## Extract data
## Load MATLAB partitions (for GenLouvain)
labels <- read.table(paste0(res_path, "/", "L.txt"), sep = ",")

## K-means
#set.seed()
k <- 3
kmeans_res <- kmeans(t(labels), k, algorithm = "Lloyd")
clusters   <- kmeans_res$cluster
(t <- table(clusters))

## Extract survival profiles
type_01 <- which(clusters == 1)
type_02 <- which(clusters == 2)
type_03 <- which(clusters == 3)

## Extract survival IDs and death IDs
surv_01  <- id$Survival[type_01]
surv_02  <- id$Survival[type_02]
surv_03  <- id$Survival[type_03]

death_01 <- id$Death[type_01]
death_02 <- id$Death[type_02]
death_03 <- id$Death[type_03]

## Collect dataframe
d1 <- data.frame(surv_01, death_01); colnames(d1) <- c("surv","death")
d2 <- data.frame(surv_02, death_02); colnames(d2) <- c("surv","death")
d3 <- data.frame(surv_03, death_03); colnames(d3) <- c("surv","death")

## Divide datafram for Surv and Death
d1_ord <- data.frame(d1$surv[order(d1$surv)]/30, d1$death[order(d1$surv)]); colnames(d1_ord) <- c("surv","death")
d2_ord <- data.frame(d2$surv[order(d2$surv)]/30, d2$death[order(d2$surv)]); colnames(d2_ord) <- c("surv","death")
d3_ord <- data.frame(d3$surv[order(d3$surv)]/30, d3$death[order(d3$surv)]); colnames(d3_ord) <- c("surv","death")

## Collect
d <- cbind(rbind(d1_ord, d2_ord, d3_ord),
           c(rep(1,t[1]), rep(2,t[2]), c(rep(3,t[3]))));colnames(d)[3] = "Subtype"

## Survival Analysis
fit_survival <- survfit(Surv(surv, death) ~ Subtype , data = d)

## Visualization with GG Plot
(p1 <- ggsurvplot(fit_survival, data = d, xlim = c(0, 120),
                  ggtheme = theme_bw(), palette = c("green", "blue", "red"),
                  legend.title = "Subtype",
                  legend.labs = c("01", "02", "03"),
                  xlab = "Survival Time (Months)",
                  ylab = "Survival Probability",
                  title = "KIRC Survival Analysis",
                  subtitle = "3 layers: mRNA Expression, DNA Methylation, miRNA",
                  legend = "right"))
## Significance
print(surv_pvalue(fit_survival)[,2])

## ***********************************************************
## Extract clinical fetures
df_clinical <- read.delim(paste0(data_path, "/",
                                 "KIDNEY_Clinical.csv"),
                          sep = ",", header = T)

## Extract for each subtype the patient ID
patient_id_type_01 <- id$PatientID[type_01]
patient_id_type_02 <- id$PatientID[type_02]
patient_id_type_03 <- id$PatientID[type_03]

## Extract only common string
patient_id_type_01 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_01)
patient_id_type_02 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_02)
patient_id_type_03 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_03) 

## Define the list of subtype IDs and corresponding patient IDs
subtype_patient_ids <- list(
  "01" = patient_id_type_01,
  "02" = patient_id_type_02,
  "03" = patient_id_type_03)

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

## Combine all data frames into a single data frame
df_patients <- bind_rows(subtype_dfs)
## Kruksal test
run_kruskal_test("Gender", df_patients)
run_kruskal_test("Age_At_Index", df_patients)
run_kruskal_test("Sync_Malignancy", df_patients) ## OK
run_kruskal_test("Stage", df_patients) ## OK
run_kruskal_test("Origin_Point", df_patients)
run_kruskal_test("Primary_Diagnosis", df_patients)
run_kruskal_test("Prior_Treatment", df_patients)
run_kruskal_test("Site_Biopsy", df_patients)
run_kruskal_test("Race", df_patients)
run_kruskal_test("Ethnicity", df_patients)
run_kruskal_test("Vital_Status", df_patients) ## OK
run_kruskal_test("Ajcc_Stating_System_Edition", df_patients)
run_kruskal_test("Ajcc_Pathologic_T", df_patients) ## OK
run_kruskal_test("Ajcc_pathologic_N", df_patients)
run_kruskal_test("Ajcc_Pathologic_M", df_patients)

## Analytics
## Sync Malignancy
(p1 <- ggplot(df_patients, aes(x = Subtype,
                               fill = Sync_Malignancy)) +
    geom_bar(sta = "count", position = "dodge") +
    scale_fill_manual(values = c("No" = "#F8766D",
                                 "Not Reported" = "#D39200"),
                      name = "Sync Malignancy") +
    labs(title = "Sync Malignancy Distribution",
         x = "Subtype", y = "Count"))
## Percentages
df_patients %>%
  group_by(Subtype, Sync_Malignancy) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>%
  arrange(Subtype, desc(count))


## Stage Distribution
(p2 <- ggplot(df_patients, aes(x = Subtype,
                               fill = Stage)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(values = c("Stage I" ="#F8766D",
                                 "Stage II"= "#D39200",
                                 "Stage III" = "#93AA00",
                                 "Stage IV" = "#00BA38"),
                      name = "Stage") +
    labs(title = "Stage Distribution",
         x = "Subtype", y = "Count"))
## Percentages
df_patients %>%
  group_by(Subtype, Stage) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>%
  arrange(Subtype, desc(count))

## Vital_Status Distribution
(p3 <- ggplot(df_patients, aes(x = Subtype,
                               fill = Vital_Status)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(values = c("Alive" ="#F8766D",
                                 "Dead"= "#D39200"),
                      name = "Vital Status") +
    labs(title = "Vital Status Distribution",
         x = "Subtype", y = "Count"))
## Percentages
df_patients %>%
  group_by(Subtype, Vital_Status) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>%
  arrange(Subtype, desc(count))

## Ajcc Pathologic T Distribution
(p4 <- ggplot(df_patients, aes(x = Subtype,
                               fill = Ajcc_Pathologic_T)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(values = c("T1"  = "#1f77b4",
                                 "T1a" = "#ff7f0e",
                                 "T1b" = "#2ca02c",
                                 "T2"  = "#d62728",
                                 "T2a" = "#9467bd",
                                 "T2b" = "#8c564b",
                                 "T3"  = "#e377c2",
                                 "T3a" = "#7f7f7f",
                                 "T3b" = "#bcbd22",
                                 "T4"  = "#17becf"),
                      name = "Ajcc Pathologic T") +
    labs(title = "Ajcc Pathologic T Distribution",
         x = "Subtype", y = "Count"))
## Percentages
df_patients %>%
  group_by(Subtype, Ajcc_Pathologic_T) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>%
  arrange(Subtype, desc(count))

cat("End Significance Generalized Louvain\n")