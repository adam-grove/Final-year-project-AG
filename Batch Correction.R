# Libraries
#install.packages("gridExtra")


library(S4Vectors)
# library(SummarizedExperiment)
library(pmp)
library(ggplot2)
library(reshape2)
library(gridExtra)

# Reading in the data and feature engineering ----
my_data <- read.csv("C:/Users/adamg/OneDrive - The University of Liverpool/Project Work/DIMSpy/Runs/DIMSpy_SNR2.0_ppm223_NC_Feature_Matrix.csv")
#my_data <- t(my_data)
#colnames(my_data) <- my_data[1, ]  # Set first row as column names
#my_data <- my_data[-1, ]           # Remove the first row (optional)
my_data <-data.frame(my_data, row.names = 1)
my_data <- as.matrix(my_data, mode = "numeric")



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

# Get column names excluding the first one (assuming data is your dataframe)
column_names <- colnames(my_data)

# Apply the function to all column names except the first
extracted_info <- lapply(column_names, extract_info)

# Access the extracted information (separate vectors for clarity)
batch <- sapply(extracted_info, function(x) x$batch)
class <- sapply(extracted_info, function(x) x$animal)
sample_order <- seq_along(colnames(my_data))




# Pmp Correction ----
data <- filter_peaks_by_fraction(df=my_data, classes=class, method="QC",
                                 qc_label="QC", min_frac=0.8)

corrected_data <- QCRSC(df=data, order=sample_order, batch=batch, 
                        classes=class, spar=0, minQC=4)

# Batch Correction----
data <- my_data
corrected_data <- corrected_data

manual_color = c("#386cb0", "#ef3b2c", "#7fc97f", "#fdb462", "#984ea3", 
                 "#a6cee3", "#778899", "#fb9a99", "#ffff33")

pca_data <- pqn_normalisation(data, classes=class, qc_label="QC")
pca_data <- mv_imputation(pca_data, method="KNN", k=5, rowmax=0.5,
                          colmax=1.0, maxp=NULL, check_df=FALSE)
pca_data <- glog_transformation(pca_data, classes=class, qc_label="QC")



pca_corrected_data <- pmp::pqn_normalisation(corrected_data, classes=class,
                                             qc_label="QC")

pca_corrected_data <- filter_samples_by_mv(pca_corrected_data, max_perc_mv = 0.8, classes = class, remove_samples = TRUE)


pca_corrected_data <- pmp::mv_imputation(pca_corrected_data, method="KNN", k=5,
                                         rowmax=0.5, colmax=0.5, maxp=NULL, check_df=FALSE)
#> Warning in knnimp(x, k, maxmiss = rowmax, maxp = maxp): 8 rows with more than 50 % entries missing;
#>  mean imputation used for these rows

# Get column names excluding the first one 
column_names <- colnames(pca_corrected_data)

# Apply the function to all column names except the first
extracted_info <- lapply(column_names, extract_info)

# Access the extracted information (separate vectors for clarity)
pca_batch <- sapply(extracted_info, function(x) x$batch)
pca_class <- sapply(extracted_info, function(x) x$animal)

pca_corrected_data <- pmp::glog_transformation(pca_corrected_data, 
                                               classes=pca_class, qc_label="QC")


# Outputting the data
output_data <- as.data.frame(pca_corrected_data)  # Convert to data frame
#output_folder <- "C:/Users/adamg/OneDrive - The University of Liverpool/Project Work/PCAs/Batch Correction PCA data/"
#file_name <- "MALDIquant_Tol5e-6_keep400_C_Filterd_BC_Feature_Matrix.csv"

full_file_path <- paste0(output_folder, file_name)
write.csv(output_data, full_file_path)


#View(pca_corrected_data)

pca_data_transposed <- t(pca_data)
pca_data <- prcomp(pca_data_transposed, center=TRUE, scale=FALSE)


pca_corrected_data_transposed <- t(pca_corrected_data)
pca_corrected_data <- prcomp(pca_corrected_data_transposed, center=TRUE, scale=FALSE)
                             
# Plotting PCA ----

# Calculate percentage of explained variance of the first two PC's
exp_var_pca <- round(((pca_data$sdev^2)/sum(pca_data$sdev^2)*100)[1:2],2)
exp_var_pca_corrected <- round(((pca_corrected_data$sdev^2) /
                                  sum(pca_corrected_data$sdev^2)*100)[1:2],2)

plots <- list()

plotdata <- data.frame(PC1=pca_data$x[, 1], PC2=pca_data$x[, 2], 
                       batch=as.factor(batch), class=class)

plots[[1]] <- ggplot(data=plotdata, aes(x=PC1, y=PC2, col=batch)) +
  geom_point(size=2) + theme(panel.background=element_blank()) +
  scale_color_manual(values=manual_color) +
  ggtitle("PCA scores, before correction") +
  xlab(paste0("PC1 (", exp_var_pca[1] ," %)")) +
  ylab(paste0("PC2 (", exp_var_pca[2] ," %)"))

plots[[2]] <- ggplot(data=plotdata, aes(x=PC1, y=PC2, col=class)) +
  geom_point(size=2) + theme(panel.background=element_blank()) +
  scale_color_manual(values=manual_color) +
  ggtitle("PCA scores, before correction") +
  xlab(paste0("PC1 (", exp_var_pca[1] ," %)")) +
  ylab(paste0("PC2 (", exp_var_pca[2] ," %)"))
plots[[2]]
plotdata_corr <- data.frame(PC1=pca_corrected_data$x[, 1], 
                            PC2=pca_corrected_data$x[, 2], batch=as.factor(pca_batch), class=pca_class)

plots[[3]] <- ggplot(data=plotdata_corr, aes(x=PC1, y=PC2, col=batch)) +
  geom_point(size=2) +
  theme(panel.background=element_blank()) +
  scale_color_manual(values=manual_color) +
  ggtitle("PCA scores, after correction") +
  xlab(paste0("PC1 (", exp_var_pca_corrected[1] ," %)")) +
  ylab(paste0("PC2 (", exp_var_pca_corrected[2] ," %)"))

plots[[4]] <- ggplot(data=plotdata_corr, aes(x=PC1, y=PC2, col=class)) +
  geom_point(size=2) +
  theme(panel.background=element_blank()) +
  scale_color_manual(values=manual_color) +
  ggtitle("PCA scores, after correction") +
  xlab(paste0("PC1 (", exp_var_pca_corrected[1] ," %)")) +
  ylab(paste0("PC2 (", exp_var_pca_corrected[2] ," %)"))

grid.arrange(ncol=2, plots[[1]], plots[[2]], plots[[3]], plots[[4]])




#output_data <- as.data.frame(pca_corrected_data$x)
#output_folder <- "C:/Users/adamg/OneDrive - The University of Liverpool/Project Work/PCAs/Batch Correction PCA data/"
#file_name <- "MZmine_MinF7.5E3_NT3.5_Tol555_NC_BC_PCAdata.csv"

#full_file_path <- paste0(output_folder, file_name)
#write.csv(output_data, full_file_path)
