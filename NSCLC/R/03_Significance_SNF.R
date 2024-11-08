## This code is useful to analyze the clusters' significance for the
## partitions obtained using MATLAB
rm(list = ls())

## Load Useful
source("01_Load_Data.R")

## Path
res_path  <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/NSCLC/Res"
data_path <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/code/NSCLC/Data"

## Extract data
id <- df_clinical_final

## Load MATLAB partitions (for SNF)
labels <- read.table(paste0(res_path, "/", "labg0.txt"))
## Count
(t <- table(labels))

## Extract survival profiles
type_01 <- which(labels == 1)
type_02 <- which(labels == 2)
# type_03 <- which(labels == 3)
## Extract survival IDs and death IDs
surv_01   <- id$Harmonized_OS_Days[type_01]
surv_02   <- id$Harmonized_OS_Days[type_02]
# surv_03   <- id$Harmonized_OS_Days[type_03]
status_01 <- id$Harmonized_OS_Event[type_01]
status_02 <- id$Harmonized_OS_Event[type_02]
# status_03 <- id$Harmonized_OS_Event[type_03]

## Collect dataframe
d1 <- data.frame(surv_01, status_01); colnames(d1) <- c("surv","status")
d2 <- data.frame(surv_02, status_02); colnames(d2) <- c("surv","status")
# d3 <- data.frame(surv_03, status_03); colnames(d3) <- c("surv","status")
## Divide datafram for Surv and Death
d1_ord <- data.frame(d1$surv[order(d1$surv)], d1$status[order(d1$surv)]); colnames(d1_ord) <- c("surv","status")
d2_ord <- data.frame(d2$surv[order(d2$surv)], d2$status[order(d2$surv)]); colnames(d2_ord) <- c("surv","status")
# d3_ord <- data.frame(d3$surv[order(d3$surv)], d3$status[order(d3$surv)]); colnames(d3_ord) <- c("surv","status")

## Collect
d <- cbind(rbind(d1_ord, d2_ord),
           c(rep(1,t[1]), rep(2,t[2]))) ; colnames(d)[3] = "Subtype"

## Survival Analysis
fit_survival <- survfit(Surv(surv, status) ~ Subtype , data = d)

## Alternative Visualization with GG Plot
ggsurvplot(fit_survival,
           data = d, xlim = c(min(d$surv, na.rm = T),
                              max(d$surv, na.rm = T)),
           ylim = c(0, 1),
           ggtheme = theme_bw(), palette = c("red", "green"),
           legend.title = "Subtype",
           legend.labs = c("01", "02"),
           xlab = "Survival Time (Days)",
           ylab = "Survival Probability",
           title = "NSCLC Survival Analysis - SNF ",
           subtitle = "5 layers: Transcriptomic, Mutations, {Immune, Myeloid, Curated} Signatures",
           legend = "right")

## Significance
surv_pvalue(fit_survival)[,2]

## ***********************************************************
## Extract clinical fetures
df_clinical <- read.delim(paste0(data_path, "/",
                                 "GLIO_Clinical.csv"),
                          sep = ",", header = T)

## Extract for each subtype the patient ID
patient_id_type_01 <- id$PatientID[type_01]
patient_id_type_02 <- id$PatientID[type_02]
patient_id_type_03 <- id$PatientID[type_03]

## Extract only common string
patient_id_type_01 <- sub("^(TCGA-[0-9]+-[0-9]+).*", "\\1",
                          patient_id_type_01)  
patient_id_type_02 <- sub("^(TCGA-[0-9]+-[0-9]+).*", "\\1",
                          patient_id_type_02)  
patient_id_type_03 <- sub("^(TCGA-[0-9]+-[0-9]+).*", "\\1",
                          patient_id_type_03)  

## Define the list of subtype IDs and corresponding patient IDs
subtype_patient_ids <- list(
  "Type 01" = patient_id_type_01,
  "Type 02" = patient_id_type_02,
  "Type 03" = patient_id_type_03)

## Generate data frames for each subtype
subtype_dfs <- lapply(names(subtype_patient_ids), function(subtype_id) {
  generate_subtype_df(subtype_id,
                      subtype_patient_ids[[subtype_id]],
                      df_clinical)
})

## Median age
median(subtype_dfs[[1]]$Age_At_Index)
median(subtype_dfs[[2]]$Age_At_Index)
median(subtype_dfs[[3]]$Age_At_Index)

## Distributions
percentages_df <- data.frame(Gender = c("Female", "Male"),
                             Type_01 = c(round(prop.table(table(subtype_dfs[[1]]$Gender))[[1]],3)*100,
                                         round(prop.table(table(subtype_dfs[[1]]$Gender))[[2]],3)*100),
                             Type_02 = c(round(prop.table(table(subtype_dfs[[2]]$Gender))[[1]],3)*100,
                                         round(prop.table(table(subtype_dfs[[2]]$Gender))[[2]],3)*100),
                             Type_03 = c(round(prop.table(table(subtype_dfs[[3]]$Gender))[[1]],3)*100,
                                         round(prop.table(table(subtype_dfs[[3]]$Gender))[[2]],3)*100))
## Check
percentages_df

## Combine all data frames into a single data frame
df_patients <- bind_rows(subtype_dfs)
## Kruksal test
kruskal.test(Gender ~ Subtype, data = df_patients)
kruskal.test(Age_At_Index ~ Subtype, data = df_patients)

## Analytics
## Age Distribution
ggplot(df_patients, aes(x = Subtype, y = Age_At_Index,
                        fill = Subtype)) +
  geom_boxplot() +
  geom_jitter(size = 1, alpha = 0.5) +
  labs(title = "Age Distribution", subtitle = "Median Age: 56",
       x = "Subtype", y = "Age")

## Gender Distribution
ggplot(df_patients, aes(x = Subtype,
                        fill = Gender)) +
  geom_bar(position = "dodge") +
  labs(title = "Gender Distribution", x = "Gender", y = "Count")
