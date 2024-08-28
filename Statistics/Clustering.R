library(igraph)
library(dplyr)

cluster_mass_spectrometry_features <- function(df, correlation_threshold) {
  # Transpose the DataFrame to have features as columns
  df_transposed <- t(df)
  
  # Compute the correlation matrix
  corr_matrix <- cor(df_transposed)
  distance_matrix <- 1 - corr_matrix
  
  # Create the adjacency matrix based on the correlation threshold
  adjacency_matrix <- as.matrix(corr_matrix >= correlation_threshold)
  diag(adjacency_matrix) <- 0  # Ensure no self-loop
  
  # Convert the adjacency matrix to a graph
  g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
  
  # Find connected components
  clusters <- components(g)$membership
  sizes <- table(Clusters)
  
  # Number of clusters
  num_clusters <- length(unique(clusters))
  cat("Number of clusters:", num_clusters, "\n")
  cat("Number of features:", ncol(df_transposed))
  # Add cluster as a column to the original DataFrame
  df$Cluster <- clusters
  return(df)
}

Cluster_dist <- function(df, correlation_threshold) {
  # Transpose the DataFrame to have features as columns
  df_transposed <- t(df)
  
  # Compute the correlation matrix
  corr_matrix <- cor(df_transposed)
  distance_matrix <- 1 - corr_matrix
  
  # Create the adjacency matrix based on the correlation threshold
  adjacency_matrix <- as.matrix(corr_matrix >= correlation_threshold)
  diag(adjacency_matrix) <- 0  # Ensure no self-loop
  
  # Convert the adjacency matrix to a graph
  g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
  
  # Find connected components
  clusters <- components(g)$membership
  
  # Number of clusters
  num_clusters <- length(unique(clusters))
  cat("Number of clusters:", num_clusters, "\n")
  cat("Number of features:", ncol(df_transposed))
  
  # Add cluster as a column to the original DataFrame
  df$Cluster <- clusters
  
  # Calculate cluster sizes
  cluster_sizes <- table(clusters)
  
  # Additional cluster analysis (optional)
  # For example, calculate average intensity per cluster
  cluster_means <- aggregate(df, list(clusters), mean)
  
  # Return a list with the DataFrame, cluster sizes, and potentially other metrics
  return(list(
    data = df,
    cluster_sizes = cluster_sizes,
    cluster_means = cluster_means  # Or other calculated metrics
  ))
}


find_max_rows <- function(df) {
  
  # Extract value columns (all but the last)
  value_cols <- names(df)[-ncol(df)]
  
  # Add a row_id column to track original indices
  df <- df %>%
    mutate(row_id = row_number(),
           feature_id = rownames(df))
  
  # Calculate maximum value for each row based on value_cols
  df <- df %>%
    rowwise() %>%
    mutate(max_value = max(across(all_of(value_cols)))) %>%
    ungroup()
  
  # Find rows with maximum value in each cluster, prioritizing first occurrence
  max_rows <- df %>%
    group_by(Cluster) %>%
    filter(max_value == max(max_value)) %>%
    arrange(Cluster) %>%  # Ensure consistent ordering within clusters
    slice_head(n = 1) %>%  # Select the first row with maximum value
    ungroup()
  
  # Return the original rows based on row_id, without unnecessary columns
  max_rows <- df[max_rows$row_id, !(names(df) %in% c("Cluster", "row_id", "max_value"))]
  max <- max_rows %>% select(feature_id, everything())
  

  return(max)
}
DIMSpy_df <- read.csv("DIMSpy_SNR3.5_ppm555_NC_BC_Feature_Matrix.csv", row.names = 1)
MALDIquant_df <- read.csv("MALDIquant_Tol5e-6_keep400_C_Filtered_BC_Feature_Matrix.csv",row.names = 1)
MZmine_df <-read.csv("MZmine_MinF7.5E3_NT3.5_Tol555_NC_BC_Feature_Matrix.csv",row.names = 1)

cluster_DIMSpy_df <- cluster_mass_spectrometry_features(DIMSpy_df,0.8)
cluster_MALDIquant_df <- cluster_mass_spectrometry_features(MALDIquant_df,0.8)
cluster_MZmine_df <- cluster_mass_spectrometry_features(MZmine_df,0.8)

DIMSpy_dist <- Cluster_dist(DIMSpy_df,0.8)
DIMSpy_dist <- DIMSpy_dist[[2]]
DIMSpy_dist <- sort(DIMSpy_dist,decreasing = TRUE)

MALDIquant_dist <- Cluster_dist(MALDIquant_df,0.8)
MALDIquant_dist <- MALDIquant_dist[[2]]
MALDIquant_dist <- sort(MALDIquant_dist,decreasing = TRUE)
View(MALDIquant_dist)

MZmine_dist <- Cluster_dist(MZmine_df, 0.8)
MZmine_dist <- MZmine_dist[[2]]
MZmine_dist <- sort(MZmine_dist, decreasing = TRUE)

DIMSpy_bp <- find_max_rows(cluster_DIMSpy_df)
MALDIquant_bp <- find_max_rows(cluster_MALDIquant_df)
MZmine_bp <- find_max_rows(cluster_MZmine_df)

write.csv(DIMSpy_bp,"DIMSpy_SNR3.5_ppm555_NC_BC_Clustered_Feature_Matrix.csv",row.names = F)
write.csv(MALDIquant_bp,"MALDIquant_Tol5e-6_keep400_C_Filtered_BC_Clustered_Feature_Matrix.csv",row.names = F)
write.csv(MZmine_bp,"MZmine_MinF7.5E3_NT3.5_Tol555_NC_BC_Clustered_Feature_Matrix.csv",row.names = F)
