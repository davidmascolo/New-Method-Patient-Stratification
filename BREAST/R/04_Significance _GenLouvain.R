## This code is useful to analyze the clusters' significance for the
## partitions obtained using MATLAB
cat("Starting Significance Generalized Louvain\n")

## Path
data_path <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/BREAST/Data/"
res_path  <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/BREAST/Res/"

## Load data
BREAST_Survival          <- read.delim(paste0(data_path, "BREAST_Survival.txt"), header = T)
BREAST_Survival$PatientIDModified <- gsub("TCGA-([0-9]{2})-([0-9]{4}).*",
                                       "\\1.\\2", BREAST_Survival$PatientID)
id <- BREAST_Survival
## Extract data
## Load MATLAB partitions (for GenLouvain)
labels <- read.table(paste0(res_path, "/", "L.txt"), sep = ",")

## K-means
#set.seed()
k <- 5
kmeans_res <- kmeans(t(labels), k, algorithm = "Lloyd")
clusters   <- kmeans_res$cluster
(t <- table(clusters))

## Extract survival profiles
type_01 <- which(clusters == 1)
type_02 <- which(clusters == 2)
type_03 <- which(clusters == 3)
type_04 <- which(clusters == 4)
type_05 <- which(clusters == 5)
## Extract survival IDs and death IDs
surv_01  <- id$Survival[type_01]
surv_02  <- id$Survival[type_02]
surv_03  <- id$Survival[type_03]
surv_04  <- id$Survival[type_04]
surv_05  <- id$Survival[type_05]
death_01 <- id$Death[type_01]
death_02 <- id$Death[type_02]
death_03 <- id$Death[type_03]
death_04 <- id$Death[type_04]
death_05 <- id$Death[type_05]

## Collect dataframe
d1 <- data.frame(surv_01, death_01); colnames(d1) <- c("surv","death")
d2 <- data.frame(surv_02, death_02); colnames(d2) <- c("surv","death")
d3 <- data.frame(surv_03, death_03); colnames(d3) <- c("surv","death")
d4 <- data.frame(surv_04, death_04); colnames(d4) <- c("surv","death")
d5 <- data.frame(surv_05, death_05); colnames(d5) <- c("surv","death")

## Divide datafram for Surv and Death
d1_ord <- data.frame(d1$surv[order(d1$surv)]/30, d1$death[order(d1$surv)]); colnames(d1_ord) <- c("surv","death")
d2_ord <- data.frame(d2$surv[order(d2$surv)]/30, d2$death[order(d2$surv)]); colnames(d2_ord) <- c("surv","death")
d3_ord <- data.frame(d3$surv[order(d3$surv)]/30, d3$death[order(d3$surv)]); colnames(d3_ord) <- c("surv","death")
d4_ord <- data.frame(d4$surv[order(d4$surv)]/30, d4$death[order(d4$surv)]); colnames(d4_ord) <- c("surv","death")
d5_ord <- data.frame(d5$surv[order(d5$surv)]/30, d5$death[order(d5$surv)]); colnames(d5_ord) <- c("surv","death")

## Collect
d <- cbind(rbind(d1_ord, d2_ord, d3_ord, d4_ord, d5_ord),
           c(rep(1,t[1]), rep(2,t[2]), c(rep(3,t[3]),
             rep(4,t[4]), rep(5,t[5]))));colnames(d)[3] = "Subtype"

## Survival Analysis
fit_survival <- survfit(Surv(surv, death) ~ Subtype , data = d)

## Alternative Visualization with GG Plot
p1 <- ggsurvplot(fit_survival, data = d, xlim = c(0, 200),
                 ggtheme = theme_bw(), palette = c("green", "blue", "red",
                                                   "purple", "yellow"),
                 legend.title = "Subtype",
                 legend.labs = c("01", "02", "03", "04", "05"),
                 xlab = "Survival Time (Months)",
                 ylab = "Survival Probability",
                 title = "BREAST Survival Analysis",
                 subtitle = "3 layers: mRNA Expression, DNA Methylation, miRNA",
                 legend = "right")
print(p1)
## Significance
print(surv_pvalue(fit_survival)[,2])

## ***********************************************************
## Extract clinical fetures
df_clinical <- read.delim(paste0(data_path, "/",
                                 "BREAST_Clinical.csv"),
                          sep = ",", header = T)

## Extract for each subtype the patient ID
patient_id_type_01 <- id$PatientID[type_01]
patient_id_type_02 <- id$PatientID[type_02]
patient_id_type_03 <- id$PatientID[type_03]
patient_id_type_04 <- id$PatientID[type_04]
patient_id_type_05 <- id$PatientID[type_05]

## Extract only common string
patient_id_type_01 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_01)
patient_id_type_02 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_02)
patient_id_type_03 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_03) 
patient_id_type_04 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_04)
patient_id_type_05 <- gsub("^(([^-]*-){2}[^-]*).*", "\\1", patient_id_type_05) 

## Define the list of subtype IDs and corresponding patient IDs
subtype_patient_ids <- list(
  "01" = patient_id_type_01,
  "02" = patient_id_type_02,
  "03" = patient_id_type_03,
  "04" = patient_id_type_04,
  "05" = patient_id_type_05)

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
print(median(subtype_dfs[[5]]$Age_At_Index))

## Combine all data frames into a single data frame
df_patients <- bind_rows(subtype_dfs)
## Kruksal test
run_kruskal_test("Gender", df_patients)
run_kruskal_test("Age_At_Index", df_patients)
run_kruskal_test("Sync_Malignancy", df_patients)
run_kruskal_test("Stage", df_patients)
run_kruskal_test("Origin_Point", df_patients)
run_kruskal_test("Primary_Diagnosis", df_patients) ## OK
run_kruskal_test("Prior_Treatment", df_patients)
run_kruskal_test("Site_Biopsy", df_patients)
run_kruskal_test("Race", df_patients)
run_kruskal_test("Ethnicity", df_patients)
run_kruskal_test("Vital_Status", df_patients)

## Analytics
## Primary Diagnosis Distribution
(p2 <- ggplot(df_patients, aes(x = Subtype,
                               fill = Primary_Diagnosis)) +
    geom_bar(stat = "count", position = "dodge") +
    geom_text(aes(label = paste0(round(after_stat(count)/sum(after_stat(count)) * 100, 2), "%")),
              stat = "count", vjust = -0.3, hjust = 0.3, size = 4) +
    scale_fill_manual(values = c("Apocrine adenocarcinoma" ="#F8766D",
                                 "Cribriform carcinoma, NOS"= "#D39200",
                                 "Infiltrating duct and lobular carcinoma" = "#93AA00",
                                 "Infiltrating duct carcinoma, NOS" = "#00BA38",
            "Infiltrating duct mixed with other types of carcinoma" = "#00C19F",
                                 "Lobular carcinoma, NOS" = "#00B9E3",
                                 "Metaplastic carcinoma, NOS" = "#619CFF"),
                      name = "Primary Diagnosis") +
    labs(title = "Primary Diagnosis Distribution",
         x = "Subtype", y = "Count") +
    coord_flip())

## See percentages
df_patients %>%
  group_by(Subtype, Primary_Diagnosis) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>%
  arrange(Subtype, desc(count))

cat("End Significance Generalized Louvain\n")
