library(S4Vectors)
# library(SummarizedExperiment)
library(pmp)
library(ggplot2)
library(reshape2)
library(gridExtra)

# Reading in the data and feature engineering ----
DIMSpy_df <- read.csv("C:/Users/adamg/OneDrive - The University of Liverpool/Project Work/DIMSpy/Runs/DIMSpy_SNR3.5_ppm555_NC_Feature_Matrix.csv")
MALDIquant_df <- read.csv("C:/Users/adamg/OneDrive - The University of Liverpool/Project Work/MALDIquant/Runs/MALDIquant_Tol5e-6_keep400_C_Filtered_Feature_Matrix.csv")
MZmine_df <- read.csv("C:/Users/adamg/OneDrive - The University of Liverpool/Project Work/MZmine/Runs/MZmine_MinF7.5E3_NT3.5_Tol555_NC_Feature_Matrix.csv")

DIMSpy_df <-data.frame(DIMSpy_df, row.names = 1)
DIMSpy_df <- as.matrix(DIMSpy_df, mode = "numeric")
MALDIquant_df <-data.frame(MALDIquant_df, row.names = 1)
MALDIquant_df <- as.matrix(MALDIquant_df, mode = "numeric")
MZmine_df <-data.frame(MZmine_df, row.names = 1)
MZmine_df <- as.matrix(MZmine_df, mode = "numeric")

# View(my_data)
# Function to extract information from column names
extract_info <- function(colname) {
  # Split the column name by underscores
  split_name <- strsplit(colname, "_")[[1]]
  
  # Extract batch number (get last character)
  batch_number <- substr(split_name[1], nchar(split_name[1]), nchar(split_name[1]))
  
  # Extract animal code (remove numbers)
  animal_code <- gsub("\\d", "", split_name[2])
  
  # Return a list with extracted information in desired format
  return(list(batch = batch_number, animal = animal_code))
}

column_names <- colnames(DIMSpy_df)
extracted_info <- lapply(column_names, extract_info)

# Access the extracted information (separate vectors for clarity)
DIMSpy_batch <- sapply(extracted_info, function(x) x$batch)
DIMSpy_class <- sapply(extracted_info, function(x) x$animal)
DIMSpysample_order <- seq_along(colnames(DIMSpy_df))

column_names <- colnames(MALDIquant_df)
extracted_info <- lapply(column_names, extract_info)

# Access the extracted information (separate vectors for clarity)
MALDIquant_batch <- sapply(extracted_info, function(x) x$batch)
MALDIquant_class <- sapply(extracted_info, function(x) x$animal)
MALDIquantsample_order <- seq_along(colnames(MALDIquant_df))

column_names <- colnames(MZmine_df)
extracted_info <- lapply(column_names, extract_info)

# Access the extracted information (separate vectors for clarity)
MZmine_batch <- sapply(extracted_info, function(x) x$batch)
MZmine_class <- sapply(extracted_info, function(x) x$animal)
MZminesample_order <- seq_along(colnames(MZmine_df))

total_missing_percentage <- sum(is.na(DIMSpy_df)) / (nrow(DIMSpy_df) * ncol(DIMSpy_df)) * 100
print(total_missing_percentage)

# Normalisation ----



# DIMSpy
print(dim(DIMSpy_df))
filtered_DIMSpy_df <- filter_peaks_by_fraction(df=DIMSpy_df, classes=DIMSpy_class, method="QC",
                                 qc_label="QC", min_frac=0.8)
print(dim(filtered_DIMSpy_df))

total_missing_percentage <- sum(is.na(filtered_DIMSpy_df) | filtered_DIMSpy_df <7500) / (nrow(filtered_DIMSpy_df) * ncol(filtered_DIMSpy_df)) * 100

Norm_DIMSpy_df <- pqn_normalisation(filtered_DIMSpy_df, classes=DIMSpy_class, qc_label="QC")
Norm_KNN_DIMSpy_df <- mv_imputation(Norm_DIMSpy_df, method="KNN", k=5, rowmax=0.5,
                          colmax=1.0, maxp=NULL, check_df=FALSE)
total_missing_percentage <- sum(is.na(Norm_KNN_DIMSpy_df)) / (nrow(Norm_KNN_DIMSpy_df) * ncol(Norm_KNN_DIMSpy_df)) * 100
total_missing_percentage

DIMSpy_pca_data <- glog_transformation(Norm_KNN_DIMSpy_df, classes=DIMSpy_class, qc_label="QC")

# MALDIquant
Norm_MALDIquant_df <- pqn_normalisation(MALDIquant_df, classes=MALDIquant_class, qc_label="QC")
Norm_KNN_MALDIquant_df <- mv_imputation(Norm_MALDIquant_df, method="KNN", k=5, rowmax=0.5,
                                        colmax=1.0, maxp=NULL, check_df=FALSE)
MALDIquant_pca_data <- glog_transformation(Norm_KNN_MALDIquant_df, classes=MALDIquant_class, qc_label="QC")

# M
Norm_MZmine_df <- pqn_normalisation(MZmine_df, classes=MZmine_class, qc_label="QC")
Norm_KNN_MZmine_df <- mv_imputation(Norm_MZmine_df, method="KNN", k=5, rowmax=0.5,
                                    colmax=1.0, maxp=NULL, check_df=FALSE)
MZmine_pca_data <- glog_transformation(Norm_KNN_MZmine_df, classes=MZmine_class, qc_label="QC")

library(dplyr)
library(tidyr)
library(ggplot2)

# Boxplots -----

# Assuming DIMSpy_class contains batch information
plots <- list()

# DIMSpy
DIMSpy_df_t <- as.data.frame(t(DIMSpy_df))# Add batch information as a column
DIMSpy_df_t$batch <- DIMSpy_batch

data_long <- DIMSpy_df_t %>%
  pivot_longer(cols = -batch, names_to = "feature", values_to = "value")

plots[[1]] <- ggplot(data_long, aes(x = batch, y = value, fill = batch)) +
  geom_boxplot() +
  labs(x = "Batch", y = "Feature Value",title = "DIMSpy: Before normalisation")

# MALDIquant
MALDIquant_df_t <- as.data.frame(t(MALDIquant_df))# Add batch information as a column
MALDIquant_df_t$batch <- MALDIquant_batch

data_long <- MALDIquant_df_t %>%
  pivot_longer(cols = -batch, names_to = "feature", values_to = "value")

plots[[2]] <- ggplot(data_long, aes(x = batch, y = value, fill = batch)) +
  geom_boxplot() +
  labs(x = "Batch", y = "Feature Value",title = "MALDiquant: Before normalisation")

# MZmine
MZmine_df_t <- as.data.frame(t(MZmine_df))# Add batch information as a column
MZmine_df_t$batch <- MZmine_batch

data_long <- MZmine_df_t %>%
  pivot_longer(cols = -batch, names_to = "feature", values_to = "value")

plots[[3]] <- ggplot(data_long, aes(x = batch, y = value, fill = batch)) +
  geom_boxplot() +
  labs(x = "Batch", y = "Feature Value",title = "MZMine: Before normalisation")

# DIMSpy Data Preparation
DIMSpy_pca_data_t <- as.data.frame(t(DIMSpy_pca_data))  # Add batch information as a column
DIMSpy_pca_data_t$batch <- DIMSpy_batch

data_long <- DIMSpy_pca_data_t %>%
  pivot_longer(cols = -batch, names_to = "feature", values_to = "value")

plots[[4]] <- ggplot(data_long, aes(x = batch, y = value, fill = batch)) +
  geom_boxplot() +
  labs(x = "Batch", y = "Feature Value",,title = "DIMSpy: After normalisation") +
  ylim(9,18)

# MALDIquant Data Preparation
MALDIquant_pca_data_t <- as.data.frame(t(MALDIquant_pca_data))  # Add batch information as a column
MALDIquant_pca_data_t$batch <- MALDIquant_batch

data_long <- MALDIquant_pca_data_t %>%
  pivot_longer(cols = -batch, names_to = "feature", values_to = "value")

plots[[5]] <- ggplot(data_long, aes(x = batch, y = value, fill = batch)) +
  geom_boxplot() +
  labs(x = "Batch", y = "Feature Value",title = "MALDIquant: After normalisation") +
  ylim(9,18)

# MZmine Data Preparation
MZmine_pca_data_t <- as.data.frame(t(MZmine_pca_data)) # Add batch information as a column
MZmine_pca_data_t$batch <- MZmine_batch

data_long <- MZmine_pca_data_t %>%
  pivot_longer(cols = -batch, names_to = "feature", values_to = "value")

plots[[6]] <- ggplot(data_long, aes(x = batch, y = value, fill = batch)) +
  geom_boxplot() +
  labs(x = "Batch", y = "Feature Value",title = "MZMine: After normalisation") +
  ylim(9, 18)


grid.arrange(ncol = 3, plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]])

#grid.arrange(ncol = 3, plots[[4]],plots[[5]],plots[[6]])


plots[[1]]
plots[[2]]
plots[[3]]
