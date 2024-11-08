## This code is useful to analyze the heatmap of Generalized Louvain
## repeated for 1000 times
rm(list = ls())

## Load Useful
source("00_support_functions.R")

## Path
res_path <- "C:/Users/david/Documents/Conferences/To_Do/EMBC_2024/Patient Similarity Network/Code/Res/"

## File to load
txt_file <- "L.txt"

## Import
df <- t(read.table(paste0(res_path, txt_file), sep = ","))
## Convert
df_numeric <- as.matrix(as.data.frame(df), row.names = NULL)
## Set rownames and colnames
rownames(df_numeric) <- paste("Patient", 1:215)
colnames(df_numeric) <- paste("Iteration", 1:1000)
## Heatmap
color_palette <- c("blue1", "cadetblue2", "chartreuse",
                   "yellow1", "brown1")
x11()
heatmap(df_numeric, Colv = NA, Rowv = NULL,
        labRow = NA, labCol = NA,
        xlab = "Iterations",
        ylab = "Patients",
        main = "Generalized Louvain x 1000 times",
        col = color_palette)
legend(x = "topright", title = "Cluster",
       legend = c("1", "2", "3", "4", "5"), fill = color_palette)

## Alternative plot
heatmap(df_numeric, Rowv = NULL, Colv = NA,
        labRow = NA, labCol = NA,
        main = "Generalized Louvain 1000 times",
        ylab = "Patients", xlab = "Iterations",
        scale = "col")
