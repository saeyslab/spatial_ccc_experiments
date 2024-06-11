# Load necessary libraries
library(dplyr)

# Function to downsample the dataframe
df <- read.csv("/data/gent/vo/000/gvo00070/vsc43831/ABA_data/cell_metadata_with_cluster_annotation.csv")
target_n_cells <- 150000

# Calculate the size of the original dataframe
original_n_cells <- nrow(df)

# Calculate the proportions of each feature_matrix_label and subclass combination
proportions <- df %>%
    group_by(feature_matrix_label, subclass) %>%
    tally() %>%
    mutate(proportion = n / original_n_cells)

# Calculate the target number of cells for each combination based on the desired total number of samples
proportions <- proportions %>%
mutate(target_count = round(proportion * target_n_cells))

# Initialize an empty dataframe to hold the downsampled data
downsampled_df <- data.frame()

# Loop through each combination and sample the required number of cells
for (i in 1:nrow(proportions)) {
    batch <- proportions$feature_matrix_label[i]
    celltype <- proportions$subclass[i]
    count <- proportions$target_count[i]

    sampled_cells <- df %>%
        filter(feature_matrix_label == batch, subclass == celltype) %>%
        sample_n(size = count, replace = FALSE)

    downsampled_df <- bind_rows(downsampled_df, sampled_cells)
}

write.csv(downsampled_df, "/data/gent/vo/000/gvo00070/vsc43831/ABA_data/cell_metadata_with_cluster_annotation_downsampled.csv", row.names = FALSE)