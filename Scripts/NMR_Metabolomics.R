# Load required libraries
library(tidyverse)
library(pheatmap)
library(FactoMineR)
library(factoextra)

# Define the file path
file_path <- "F:/nmr/biofluids/NMR-based-metabonomics-for-biofluids/Data/Dummy"

# Import the files with encoding specified
file_names <- list.files(file_path, full.names = TRUE, pattern = "Sample_.*\\.csv")
data_list <- lapply(file_names, read.csv)

# Bind the data into one DataFrame
data <- bind_rows(data_list, .id = "Sample_ID")

# Add group information
data <- data %>%
  mutate(Group = if_else(Sample_ID %in% c("1", "2", "3"), "V", "P"))

# Write the combined data to a new CSV file
write.csv(data, file = file.path(file_path, "Combined_Data.csv"), row.names = FALSE)

# Filter out rows where Compound_Name is "DSS"
data <- data %>% 
  filter(Compound_Name != "DSS")

# Prepare data for PCA
pca_data <- data %>%
  select(Sample_ID, Compound_Name, Concentration_uM) %>% # Keep Sample_ID
  pivot_wider(names_from = Compound_Name, values_from = Concentration_uM) %>% # Create wide format
  column_to_rownames("Sample_ID") # Set Sample_ID as row names

# Check PCA data
print(head(pca_data))


# Perform PCA using base R
pca_results <- prcomp(pca_data, scale. = TRUE)

# Extract PCA results for plotting
pca_scores <- as.data.frame(pca_results$x) # PCA scores
pca_var <- summary(pca_results)$importance # Explained variance
pca_scores$Sample_ID <- rownames(pca_data) # Add Sample_ID for identification

# Merge with group information
pca_scores <- pca_scores %>%
  left_join(sample_groups, by = "Sample_ID")


library(ggplot2)

# Plot the first two principal components
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE) +
  theme_minimal() +
  labs(title = "PCA of Metabolomics Data",
       x = paste0("PC1 (", round(pca_var[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(pca_var[2, 2] * 100, 1), "%)"))



# Pivot wider to prepare data for heatmap
heatmap_data <- data %>%
  select(Sample_ID, Compound_Name, Concentration_uM) %>%
  pivot_wider(names_from = Sample_ID, values_from = Concentration_uM)

# Extract Group information for samples
sample_groups <- data %>%
  select(Sample_ID, Group) %>%
  distinct()

# Convert the wide data to a matrix (excluding Compound_Name)
heatmap_matrix <- as.matrix(heatmap_data %>% select(-Compound_Name))

# Ensure rows are named by Compound_Name
rownames(heatmap_matrix) <- heatmap_data$Compound_Name

# Create annotation for heatmap columns (Sample_ID -> Group)
column_annotation <- sample_groups$Group
names(column_annotation) <- sample_groups$Sample_ID


# Plot heatmap
pheatmap(heatmap_matrix,
         scale = "row", # Scale rows (standardize metabolites)
         cluster_rows = TRUE, # Cluster metabolites
         cluster_cols = TRUE, # Cluster samples
         annotation_col = data.frame(Group = column_annotation), # Add group annotation
         show_rownames = TRUE, # Show metabolite names
         show_colnames = TRUE, # Show sample IDs
         main = "Heatmap of Metabolomics Data")




# Identify top 5 metabolites for each group
top_metabolites <- data %>%
  group_by(Group, Compound_Name) %>%
  summarise(Average_Concentration = mean(Concentration_uM, na.rm = TRUE)) %>%
  arrange(Group, desc(Average_Concentration)) %>%
  group_by(Group) %>%
  slice_max(order_by = Average_Concentration, n = 5)

print(top_metabolites)

# Compare "V" vs "P" for top 5 metabolites
comparison <- top_metabolites %>%
  pivot_wider(names_from = Group, values_from = Average_Concentration) %>%
  mutate(Fold_Change = V / P)

print(comparison)

library(ggplot2)

# Identify top 5 metabolites for each group
top_metabolites <- data %>%
  group_by(Group, Compound_Name) %>%
  summarise(Average_Concentration = mean(Concentration_uM, na.rm = TRUE)) %>%
  arrange(Group, desc(Average_Concentration)) %>%
  group_by(Group) %>%
  slice_max(order_by = Average_Concentration, n = 5) %>%
  ungroup()

# Filter the original data to include only the top metabolites
top_metabolites_data <- data %>%
  filter(Compound_Name %in% top_metabolites$Compound_Name)

# Plot box plots for the top 5 metabolites
ggplot(top_metabolites_data, aes(x = Group, y = Concentration_uM, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Compound_Name, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Box Plot of Top 5 Metabolites by Group",
    x = "Group",
    y = "Concentration (µM)"
  ) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "none"
  )
