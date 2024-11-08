## This code is useful to analyze the clusters' significance for the
## partitions obtained using MATLAB

## Load Useful
source("01_Load_Data.R")

## Path
res_path  <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/NSCLC/Res"
data_path <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/NSCLC/Data"

## Extract data
id <- df_clinical_final
id$Response <- ifelse(id$Harmonized_Confirmed_BOR %in% c("CR", "PR"), 1, 0)

## Load MATLAB partitions 
labels <- read.table(paste0(res_path, "/", "L.txt"), sep = ",")

## Clustering
k <- 2
kmeans_res <- kmeans(t(labels), k, algorithm = "Lloyd")
clusters   <- kmeans_res$cluster
t          <- table(clusters)
## Extract survival profiles
type_01 <- which(clusters == 1)
type_02 <- which(clusters == 2)
## Extract survival IDs and death IDs
surv_01   <- id$Harmonized_OS_Days[type_01]
surv_02   <- id$Harmonized_OS_Days[type_02]
status_01 <- id$Harmonized_OS_Event[type_01]
status_02 <- id$Harmonized_OS_Event[type_02]
  
## Collect dataframe
d1 <- data.frame(surv_01, status_01); colnames(d1) <- c("surv","status")
d2 <- data.frame(surv_02, status_02); colnames(d2) <- c("surv","status")
## Divide datafram for Surv and Death
d1_ord <- data.frame(d1$surv[order(d1$surv)], d1$status[order(d1$surv)]); colnames(d1_ord) <- c("surv","status")
d2_ord <- data.frame(d2$surv[order(d2$surv)], d2$status[order(d2$surv)]); colnames(d2_ord) <- c("surv","status")
  
## Collect
d <- cbind(rbind(d1_ord, d2_ord),
           c(rep(1,t[1]), rep(2,t[2]))) ; colnames(d)[3] = "Subtype"
  
## Survival Analysis
fit_survival <- survfit(Surv(surv, status) ~ Subtype , data = d)

## Alternative Visualization with GG Plot
ggsurvplot(fit_survival,
           data = d, xlim = c(0,
                              max(d$surv, na.rm = T)),
           ylim = c(0, 1),
           ggtheme = theme_bw(), palette = c("red", "green"),
           legend.title = "Subtype",
           legend.labs = c("01", "02"),
           xlab = "Survival Time (Months)",
           ylab = "Survival Probability",
           title = "NSCLC Survival Analysis - Algorithm Name",
           subtitle = "5 layers: Transcriptomic, Mutations, {Immune, Myeloid, Curated} Signatures",
           legend = "right")

## Significance
surv_pvalue(fit_survival)[,2]



## ***********************************************************
## Extract for each subtype the patient ID
patient_id_type_01 <- id$Harmonized_SU2C_Participant_ID_v2[type_01]
patient_id_type_02 <- id$Harmonized_SU2C_Participant_ID_v2[type_02]


## Define the list of subtype IDs and corresponding patient IDs
subtype_patient_ids <- list(
  "01" = patient_id_type_01,
  "02" = patient_id_type_02)

## Generate data frames for each subtype
subtype_dfs <- lapply(names(subtype_patient_ids), function(subtype_id) {
  generate_subtype_df(subtype_id,
                      subtype_patient_ids[[subtype_id]],
                      id)
})

## Median age
median(subtype_dfs[[1]]$Age_At_Diagnosis)
median(subtype_dfs[[2]]$Age_At_Diagnosis)
surv_pvalue(fit_survival)[,2]


## Distributions
percentages_df_gender <- data.frame(Gender = c("Female", "Male"),
                                    Type_01 = c(round(prop.table(table(subtype_dfs[[1]]$Gender))[[1]],3)*100,
                                                round(prop.table(table(subtype_dfs[[1]]$Gender))[[2]],3)*100),
                                    Type_02 = c(round(prop.table(table(subtype_dfs[[2]]$Gender))[[1]],3)*100,
                                                round(prop.table(table(subtype_dfs[[2]]$Gender))[[2]],3)*100))
## Check
percentages_df_gender
## Combine all data frames into a single data frame
df_patients <- bind_rows(subtype_dfs)

## Kruksal test
kruskal.test(Age_At_Diagnosis ~ Subtype, data = df_patients) ## OK
kruskal.test(Gender ~ Subtype, data = df_patients) ## OK
kruskal.test(Smoking_Status ~ Subtype, data = df_patients)
kruskal.test(Stage ~ Subtype, data = df_patients)
kruskal.test(Histology ~ Subtype, data = df_patients)
kruskal.test(Line_Theraphy ~ Subtype, data = df_patients) ## OK
kruskal.test(Prior_Platinum ~ Subtype, data = df_patients) 
kruskal.test(Prior_TKI ~ Subtype, data = df_patients)
kruskal.test(Response ~ Subtype, data = df_patients) ## OK

## Analytics
## Age Distribution
ggplot(df_patients, aes(x = Subtype, y = Age_At_Diagnosis,
                        fill = Subtype)) +
  geom_boxplot() +
  geom_jitter(size = 1, alpha = 0.5) +
  labs(title = "Age Distribution", subtitle = "Median Age: 63",
       x = "Subtype", y = "Age")

## Gender Distribution
ggplot(df_patients, aes(x = Subtype,
                        fill = Gender)) +
  geom_bar(position = "dodge") +
  labs(title = "Gender Distribution", x = "Subtype", y = "Count")

## Response Distribution
ggplot(df_patients, aes(x = Subtype,
                        fill = factor(Response))) +
  geom_bar(position = "dodge") +
  labs(title = "Response Distribution", x = "Subtype", y = "Count") +
  scale_fill_brewer(palette = "Set1",
                    name = "Response")
## Comparison Metrics
ARI(id$Response, clusters)
Chi2(id$Response, clusters)

## Line Therapy
ggplot(df_patients, aes(x = Subtype,
                        fill = factor(Line_Theraphy)),
       col) +
  geom_bar(position = "dodge") +
  labs(title = "Line Of Therapy Distribution",
       x = "Subtype", y = "Count") +
  scale_fill_brewer(palette = "Set1",
                    name = "Line of Therapy")
## Distribution
percentages_df_theraphy <- data.frame(Theraphy = c("1", "2", "3", "4", "5"),
                                      Type_01 = rep(NA, 5),
                                      Type_02 = rep(NA, 5))
percentages_df_theraphy$Type_01      <- c(round(prop.table(table(subtype_dfs[[1]]$Line_Theraphy)),3)*100)
percentages_df_theraphy$Type_02[[1]] <- c(round(prop.table(table(subtype_dfs[[2]]$Line_Theraphy)),3)*100)[[1]]
percentages_df_theraphy$Type_02[[2]] <- c(round(prop.table(table(subtype_dfs[[2]]$Line_Theraphy)),3)*100)[[2]]
## Check
percentages_df_theraphy




## BONUS PLOT FOR GENDER DISTRIBUTION IN EDA
df_patients <- data.frame(Gender = c(rep("F", 207), rep("M", 182)))
df_patients$Gender <- factor(df_patients$Gender, levels = c("F", "M"))
ggplot(df_patients, aes(x = Gender, fill = Gender)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("F" = "pink", "M" = "steelblue"),
                    name = "Gender") +
  labs(title = "Gender Distribution", x = "", y = "Count") +
  theme(axis.text.x = element_blank())
