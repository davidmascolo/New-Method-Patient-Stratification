## Similarity Network Fusion

rm(list = ls())
## Load Useful
source("00_support_functions.R")
source("01_load_data.R")

## Path
res_path  <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Code/GenLouvain-master/Results"
data_path <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Personal/Data/GBM"
## Load data
GBM_Gene_Expression   <- read.delim(paste0(data_path, "/", "GLIO_Gene_Expression.txt"), header = T)
GBM_Mirna_Expression  <- read.delim(paste0(data_path, "/", "GLIO_Methy_Expression.txt"), header = T)
GBM_Methy_Expression  <- read.delim(paste0(data_path, "/", "GLIO_Mirna_Expression.txt"), header = T)
GBM_Survival          <- read.delim(paste0(data_path, "/", "GLIO_Survival.txt"), header = T)
GBM_Mutations         <- read_excel(paste0(data_path, "/", "GLIO_Mutation+CNA.xlsx"), col_names = T)
GBM_Mutations[, -1]   <- apply(GBM_Mutations[, -1], 2,
                       function(x) ifelse(is.na(x), 0, 1))
GBM_Mutations <- GBM_Mutations[order(GBM_Mutations$`Tumor #`), ]

## Preprocess data
## Gene  
data1 <- data.frame(t(GBM_Gene_Expression))
data1 <- data1[!(rownames(data1) == "TCGA.06.0147.01A.01R.0219.01"), ]
## Mirna  
data2 <- data.frame(t(GBM_Mirna_Expression))
data2 <- data2[!(rownames(data2) == "TCGA.06.0147.01A.01D.0218.05"), ]
## Methy
data3 <- data.frame(t(GBM_Methy_Expression))
data3 <- data3[!(rownames(data3) == "TCGA.06.0147.01A.01T.0213.07"), ]
## Survival
id                <- GBM_Survival
## Mutations
data4 <- GBM_Mutations; 


## Cut
data1$Patients <- sub("^TCGA.(\\d{2}.\\d{4}).*$", "\\1", rownames(data1))
data2$Patients <- sub("^TCGA.(\\d{2}.\\d{4}).*$", "\\1", rownames(data2))
data3$Patients <- sub("^TCGA.(\\d{2}.\\d{4}).*$", "\\1", rownames(data3))
data4$Patients <- gsub("\\-", ".", data4$`Tumor #`)
id$Patients    <- gsub("\\-", ".", id$PatientID)
id$Patients    <- sub("^TCGA.(\\d{2}.\\d{4}).*$", "\\1", id$Patients)


## Common patients
common_patients <- Reduce(intersect,
                          list(data1$Patients,
                               data2$Patients,
                               data3$Patients,
                               data4$Patients,
                               id$Patients))
## Filtering
data1 <- data1[data1$Patients %in% common_patients, ]
data2 <- data2[data2$Patients %in% common_patients, ]
data3 <- data3[data3$Patients %in% common_patients, ]
data4 <- data4[data4$Patients %in% common_patients, ]
id    <- id[id$Patients %in% common_patients, ]
id    <- unique(id)

## Remove columns
data1$Patients  <- NULL
data2$Patients  <- NULL
data3$Patients  <- NULL
data4$`Tumor #` <- NULL

## Check
dim(data1);dim(data2);dim(data3);dim(data4);dim(id)

## Export data
file_path_01 <- "/filtered_gene.txt"
file_path_02 <- "/filtered_mirna.txt"
file_path_03 <- "/filtered_methy.txt"
file_path_04 <- "/filtered_id_survival.txt"
file_path_05 <- "/filtered_mutations.txt"
## Save the dataframe to a text file
write.table(t(data1), file = paste0("C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Dati TUMORI/GBM Data/", file_path_01),
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(data2, file = paste0("C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Dati TUMORI/GBM Data/", file_path_02),
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(data3, file = paste0("C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Dati TUMORI/GBM Data/", file_path_03),
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(id, file = paste0("C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Dati TUMORI/GBM Data/", file_path_04),
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(data4, file = paste0("C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Dati TUMORI/GBM Data/", file_path_05),
            sep = "\t", row.names = FALSE, col.names = FALSE)

## *********************************************************
## SNF
## Set parameters
K     <- 20		
alpha <- 0.5
t     <- 20
## Standardization
data1_scaled <- standardNormalization(data1)
data2_scaled <- standardNormalization(data2)
data3_scaled <- standardNormalization(data3)
data4_scaled <- standardNormalization(data4)
## Compute Distances
Dist1 <- dist2(as.matrix(data1_scaled),as.matrix(data1_scaled))
Dist2 <- dist2(as.matrix(data2_scaled),as.matrix(data2_scaled))
Dist3 <- dist2(as.matrix(data3_scaled),as.matrix(data3_scaled))
Dist4 <- dist2(as.matrix(data4_scaled),as.matrix(data4_scaled))
## Weight Matrices
W1 <- affinityMatrix(Dist1, K, alpha)
W2 <- affinityMatrix(Dist2, K, alpha)
W3 <- affinityMatrix(Dist3, K, alpha)
W4 <- affinityMatrix(Dist4, K, alpha)
## Affinity Matrix
W <- SNF(list(W1, W2, W3, W4), K, t)
## Clustering
C     <- 3			
group <- spectralClustering(W, C)
(t     <- table(group))

## Extract survival profiles
type_01 <- which(group == 1)
type_02 <- which(group == 2)
type_03 <- which(group == 3)
## Extract survival IDs
surv_01 <- id$Survival[type_01]
surv_02 <- id$Survival[type_02]
surv_03 <- id$Survival[type_03]
## Data
df_01 <- data.frame(surv_01, id$Death[type_01]); colnames(df_01) <- c("surv", "death") 
df_02 <- data.frame(surv_02, id$Death[type_02]); colnames(df_02) <- c("surv", "death") 
df_03 <- data.frame(surv_03, id$Death[type_03]); colnames(df_03) <- c("surv", "death") 
## Data divided for Survive and Death
df_01_ord <- data.frame(df_01$surv[order(df_01$surv)]/30, df_01$death[order(df_01$surv)]); colnames(df_01_ord) = c("surv","death")
df_02_ord <- data.frame(df_02$surv[order(df_02$surv)]/30, df_02$death[order(df_02$surv)]); colnames(df_02_ord) = c("surv","death")
df_03_ord <- data.frame(df_03$surv[order(df_03$surv)]/30, df_03$death[order(df_03$surv)]); colnames(df_03_ord) = c("surv","death")
## Merge
d <- cbind(rbind(df_01_ord, df_02_ord, df_03_ord),
           c(rep(1, t[1]), rep(2, t[2]), rep(3, t[3]))); colnames(d)[3] <- "Subtype"

## Survival Analysis
fit_survival <- survfit(Surv(surv, death) ~ Subtype , data = d)
## Plot
plot(fit_survival, 
     col = c("blue", "red", "green"),
     mark.time = T, xlim = c(0,102),
     main = "GBM Survival Analysis - SNF",
     ylab = "Survival Probability",
     xlab = "Survival Time (Months)",
     lwd = 2)
legend("topright",
       legend = c("Subtype = 1", "Subtype = 2", "Subtype = 3"),
       col = c("blue", "red", "green"),
       box.lty = 0, cex = 0.7, lty = 1)

## Significance
surv_pvalue(fit_survival)[,2]

## The log rank p-value from the Cox model can be accessed by:
coxfit <- coxph(Surv(surv, death) ~ Subtype, data = d); summary(coxfit)$sctest[3] 
