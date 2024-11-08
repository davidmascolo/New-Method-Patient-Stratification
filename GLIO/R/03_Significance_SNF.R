## This code is useful to analyze the clusters' significance for the
## partitions obtained using MATLAB
cat("Starting Significance SNF\n")

## Load useful scripts
source("00_support_functions.R")
source("01_load_data.R")
source("02_EDA.R")

## Path
data_path <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/GLIO/Data/"
res_path  <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/GLIO/Res/"

## Load survival data
GBM_Survival          <- read.delim(paste0(data_path, "GLIO_Survival.txt"), header = T)
id <- GBM_Survival

## Load MATLAB partitions (for SNF)
labels <- read.table(paste0(res_path, "labg01.txt"), sep = ",")
## Count
(t <- table(labels))

## Extract survival profiles
type_01 <- which(labels == 1)
type_02 <- which(labels == 2)
type_03 <- which(labels == 3)
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
           c(rep(1,t[1]), rep(2,t[2]), rep(3,t[3]))) ; colnames(d)[3] = "Subtype"

## Survival Analysis
fit_survival <- survfit(Surv(surv, death) ~ Subtype , data = d)

## Alternative Visualization with GG Plot
p1 <- ggsurvplot(fit_survival,
                 data = d, xlim = c(0, 125),
                 ggtheme = theme_bw(), palette = c("red", "blue", "green"),
                 legend.title = "Subtype",
                 legend.labs = c("01", "02", "03"),
                 xlab = "Survival Time (Months)",
                 ylab = "Survival Probability",
                 title = "GBM Survival Analysis - SNF ",
                 subtitle = "3 layers: mRNA Expression, DNA Methylation, miRNA",
                 legend = "right")
print(p1)
        
## Significance
print(surv_pvalue(fit_survival)[,2])

## ***********************************************************
## Extract clinical fetures
df_clinical <- read.delim(paste0(data_path, "GLIO_Clinical.csv"),
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
print(median(subtype_dfs[[1]]$Age_At_Index))
print(median(subtype_dfs[[2]]$Age_At_Index))
print(median(subtype_dfs[[3]]$Age_At_Index))

## Distributions
percentages_df <- data.frame(Gender = c("Female", "Male"),
                             Type_01 = c(round(prop.table(table(subtype_dfs[[1]]$Gender))[[1]],3)*100,
                                         round(prop.table(table(subtype_dfs[[1]]$Gender))[[2]],3)*100),
                             Type_02 = c(round(prop.table(table(subtype_dfs[[2]]$Gender))[[1]],3)*100,
                                         round(prop.table(table(subtype_dfs[[2]]$Gender))[[2]],3)*100),
                             Type_03 = c(round(prop.table(table(subtype_dfs[[3]]$Gender))[[1]],3)*100,
                                         round(prop.table(table(subtype_dfs[[3]]$Gender))[[2]],3)*100))
## Check
print(percentages_df)

## Combine all data frames into a single data frame
df_patients <- bind_rows(subtype_dfs)
## Kruksal test
print(kruskal.test(Gender ~ Subtype, data = df_patients))
print(kruskal.test(Age_At_Index ~ Subtype, data = df_patients))

## Analytics
## Age Distribution
p2 <- ggplot(df_patients, aes(x = Subtype, y = Age_At_Index,
                        fill = Subtype)) +
  geom_boxplot() +
  geom_jitter(size = 1, alpha = 0.5) +
  labs(title = "Age Distribution", subtitle = "Median Age: 56",
       x = "Subtype", y = "Age")

## Gender Distribution
p3 <- ggplot(df_patients, aes(x = Subtype,
                        fill = Gender)) +
  geom_bar(position = "dodge") +
  labs(title = "Gender Distribution", x = "Gender", y = "Count")

print(p2)
print(p3)

cat("End Significance Generalized Louvain\n")
