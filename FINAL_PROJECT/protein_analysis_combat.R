library(readxl)       # version 1.4.5
library(dplyr)        # version 1.1.4
library(sva)          # version 3.54.0
library(pheatmap)     # version 1.0.12
library(RColorBrewer) # version 1.1-3
library(ggplot2)      # version 3.5.2
library(ggrepel)      # version 0.9.6
library(umap)         # version 0.2.10.0
library(tidyr)        # version 1.3.1
library(stringr)      # version 1.5.1

# --- 0. Session Info and Setup ---
print("=== SESSION INFO ===")
sessionInfo()

# File paths
protein_file <- "Table_S1_wt.xlsx"
poland_file <- "PolandYeast_Emma.csv"
mapping_file <- "idmapping_2025_03_02.tsv"

# Create output directory for all plots
output_dir <- "output_plots_combat"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- 1. Load ID Mapping File ---
# Skip the header row and read properly
id_mapping <- read_tsv(mapping_file, skip = 1, col_names = c("From", "To"), show_col_types = FALSE)
print(paste("Loaded", nrow(id_mapping), "ID mappings"))
print("Sample mappings (id_mapping):")
print(head(id_mapping))

# Debug: Check for any issues with the mapping
print("Unique 'From' values in id_mapping (first 10):")
print(head(unique(id_mapping$From), 10))
print("Unique 'To' values in id_mapping (first 10):")
print(head(unique(id_mapping$To), 10))

# --- 2. Load and Prepare Table S1 Data ---
protein_data <- read_excel(protein_file, skip = 1)
print("Table S1 columns:")
print(colnames(protein_data))

protein_values <- protein_data %>%
  select(1, 2, 5:16)

colnames(protein_values) <- c(
  "ProteinId", "Gene_Symbol", "W303_s_1", "W303_s_2", "W303_s_3",
  "BY4742_s_1", "BY4742_s_2", "BY4742_s_3",
  "W303_e_1", "W303_e_2", "W303_e_3",
  "BY4742_e_1", "BY4742_e_2", "BY4742_e_3"
)

protein_values_clean <- protein_values %>%
  filter(!is.na(ProteinId) & ProteinId != "" & ProteinId != "Protein Id") %>%
  filter(!is.na(Gene_Symbol) & Gene_Symbol != "")

for (col in 3:ncol(protein_values_clean)) {
  protein_values_clean[[col]] <- as.numeric(as.character(protein_values_clean[[col]]))
}
print(paste("Table S1: Proteins with Gene Symbols:", nrow(protein_values_clean)))
print("Sample Table S1 Gene Symbols:")
print(head(protein_values_clean$Gene_Symbol, 10))

protein_means <- protein_values_clean %>%
  mutate(
    W303_s = rowMeans(select(., W303_s_1, W303_s_2, W303_s_3), na.rm = TRUE),
    BY4742_s = rowMeans(select(., BY4742_s_1, BY4742_s_2, BY4742_s_3), na.rm = TRUE),
    W303_e = rowMeans(select(., W303_e_1, W303_e_2, W303_e_3), na.rm = TRUE),
    BY4742_e = rowMeans(select(., BY4742_e_1, BY4742_e_2, BY4742_e_3), na.rm = TRUE)
  ) %>%
  select(ProteinId, Gene_Symbol, W303_s, BY4742_s, W303_e, BY4742_e)

# --- 3. Load and Prepare Poland Yeast Data ---
poland_data_raw <- read_csv(poland_file, show_col_types = FALSE)
print("Poland Yeast columns (raw):")
print(colnames(poland_data_raw))
print("Sample Poland Yeast Accession (raw):")
print(head(poland_data_raw$Accession, 10))

poland_means <- poland_data_raw %>%
  select(Accession, `1_1`:`8_3`) %>%
  mutate(
    WT_YPD_Glucose = rowMeans(select(., `1_1`, `1_2`, `1_3`), na.rm = TRUE),
    maf1D_YPD_Glucose = rowMeans(select(., `2_1`, `2_2`, `2_3`), na.rm = TRUE),
    rpc128_YPD_Glucose = rowMeans(select(., `3_1`, `3_2`, `3_3`), na.rm = TRUE),
    WT_YPGly_30C = rowMeans(select(., `4_1`, `4_2`, `4_3`), na.rm = TRUE),
    WT_YPGly_37C = rowMeans(select(., `5_1`, `5_2`, `5_3`), na.rm = TRUE),
    maf1D_YPGly_30C = rowMeans(select(., `6_1`, `6_2`, `6_3`), na.rm = TRUE),
    maf1D_YPGly_37C = rowMeans(select(., `7_1`, `7_2`, `7_3`), na.rm = TRUE),
    rpc128_YPGly_30C = rowMeans(select(., `8_1`, `8_2`, `8_3`), na.rm = TRUE)
  ) %>%
  select(Accession, WT_YPD_Glucose:rpc128_YPGly_30C)
print(paste("Poland Yeast: Proteins loaded using Accession:", nrow(poland_means)))

# --- 4. Map Poland Yeast Accession to Gene Symbols ---
common_ids <- intersect(poland_means$Accession, id_mapping$From)
print(paste("Common Accession IDs found in mapping:", length(common_ids)))
if (length(common_ids) > 0) {
  print("Sample common Accession IDs for mapping:")
  print(head(common_ids, 10))
} else {
  print("No common Accession IDs found for mapping. Further checks needed:")
  print("Sample Accession from Poland data (first 20 unique):")
  print(head(unique(poland_means$Accession), 20))
  print("Sample 'From' IDs from id_mapping (first 20 unique):")
  print(head(unique(id_mapping$From), 20))
}

poland_means_mapped <- poland_means %>%
  left_join(id_mapping, by = c("Accession" = "From")) %>%
  rename(Gene_Symbol = To) %>%
  filter(!is.na(Gene_Symbol) & Gene_Symbol != "") %>%
  select(-Accession)

print(paste("Poland Yeast: Proteins with mapped Gene Symbols:", nrow(poland_means_mapped)))
if(nrow(poland_means_mapped) > 0) {
  print("Sample mapped Gene Symbols (Poland Yeast):")
  print(head(poland_means_mapped$Gene_Symbol, 10))
}

poland_means_unique <- poland_means_mapped %>%
  group_by(Gene_Symbol) %>%
  summarise(across(WT_YPD_Glucose:rpc128_YPGly_30C, \(x) mean(x, na.rm = TRUE)), .groups = 'drop')
print(paste("Poland Yeast: Unique Gene Symbols after deduplication:", nrow(poland_means_unique)))

# --- 5. Join Datasets on Gene Symbol ---
if(nrow(protein_means) > 0 && nrow(poland_means_unique) > 0) {
  common_genes_before_join <- intersect(protein_means$Gene_Symbol, poland_means_unique$Gene_Symbol)
  print(paste("Common Gene Symbols between TableS1 and PolandYeast (mapped & unique) before join:", length(common_genes_before_join)))
  if (length(common_genes_before_join) > 0) {
    print("Sample common Gene Symbols for joining:")
    print(head(common_genes_before_join, 10))
  }
}

combined_data <- inner_join(protein_means, poland_means_unique, by = "Gene_Symbol")
print(paste("Combined dataset: Common Gene Symbols found after join:", nrow(combined_data)))

if(nrow(combined_data) > 0) {
  print("First few rows of combined data:")
  print(head(combined_data))
  print("Column names of combined data:")
  print(colnames(combined_data))
}

if(nrow(combined_data) < 2) {
  stop("Not enough common proteins found for analysis. Check Gene Symbol mapping and ID parsing steps.")
}

combined_data_clean <- combined_data %>%
  select(-ProteinId) %>%
  filter(complete.cases(.))

print(paste("Combined dataset after removing NAs:", nrow(combined_data_clean)))

if(nrow(combined_data_clean) < 2) {
  stop("Not enough proteins remaining after cleaning (NA removal) for analysis.")
}

# --- 6. Prepare Data for Batch Correction and Downstream Analysis ---
combined_mat_raw <- as.matrix(combined_data_clean[,-1])
rownames(combined_mat_raw) <- combined_data_clean$Gene_Symbol

combined_mat_log <- log2(combined_mat_raw + 1) # Log-transform before ComBat

# Define sample annotation (used for ComBat and plotting)
sample_annotation <- data.frame(
  Sample = colnames(combined_mat_log),
  Dataset = c(rep("Table_S1", 4), rep("Poland_Yeast", ncol(poland_means_unique)-1)),
  Condition = c("W303_stationary", "BY4742_stationary", "W303_exponential", "BY4742_exponential",
                colnames(poland_means_unique)[-1]),
  row.names = colnames(combined_mat_log)
)
sample_annotation$Condition <- make.unique(sample_annotation$Condition)
sample_annotation$Batch <- as.factor(sample_annotation$Dataset) # Batch variable for ComBat

# Add a new 'Experimental_Group' column for more general growth conditions
sample_annotation$Experimental_Group <- case_when(
  str_detect(sample_annotation$Condition, "stationary") ~ "Stationary Phase",
  str_detect(sample_annotation$Condition, "exponential") ~ "Exponential Phase",
  str_detect(sample_annotation$Condition, "YPD_Glucose") ~ "YPD Glucose",
  str_detect(sample_annotation$Condition, "YPGly_30C") ~ "YPGlycerol 30°C",
  str_detect(sample_annotation$Condition, "YPGly_37C") ~ "YPGlycerol 37°C",
  TRUE ~ "Other"
)

# Ensure it's a factor
sample_annotation$Experimental_Group <- factor(sample_annotation$Experimental_Group)

# Debug: Print sample annotation structure
print("Sample annotation data frame:")
print(sample_annotation)
print("Experimental_Group unique values:")
print(unique(sample_annotation$Experimental_Group))
print("Experimental_Group table:")
print(table(sample_annotation$Experimental_Group))
print("Sample annotation column names:")
print(colnames(sample_annotation))

# --- 6a. Apply ComBat Batch Correction ---
# ComBat requires features (genes) as rows and samples as columns
# Ensure no NA/NaN/Inf values in combined_mat_log before ComBat
if(any(!is.finite(combined_mat_log))) {
    stop("Non-finite values (NA/NaN/Inf) found in log-transformed matrix before ComBat. Please check data cleaning steps.")
}

# Check if there are at least two batches
if (length(unique(sample_annotation$Batch)) < 2) {
  print("Only one batch detected. ComBat correction will be skipped.")
  combined_mat_corrected <- combined_mat_log # Use uncorrected data if only one batch
} else {
  print("Applying ComBat for batch effect correction...")
  # If covariates are confounded with batch (e.g., conditions are unique to batches),
  # use an intercept-only model for ComBat.
  mod_combat = model.matrix(~1, data=sample_annotation) # Changed from ~Condition

  # Ensure batch vector is correctly aligned with columns of combined_mat_log
  batch_vector <- sample_annotation$Batch[match(colnames(combined_mat_log), rownames(sample_annotation))]

  if(any(is.na(batch_vector))){
      stop("Mismatch between sample names in expression matrix and batch information.")
  }

  # Check for constant rows (genes with no variance across samples within any batch)
  # ComBat can fail if there are genes with zero variance within any batch.
  # A simple check: remove genes with zero variance across all samples first.
  variances_check <- apply(combined_mat_log, 1, var, na.rm = TRUE)
  combined_mat_log_filtered <- combined_mat_log[variances_check > 1e-6, , drop = FALSE]
  
  if(nrow(combined_mat_log_filtered) < nrow(combined_mat_log)){
      print(paste("Removed", nrow(combined_mat_log) - nrow(combined_mat_log_filtered), "genes with zero variance before ComBat."))
  }
  
  if(nrow(combined_mat_log_filtered) < 2){
      stop("Not enough genes with variance to perform ComBat.")
  }

  tryCatch({
    combined_mat_corrected <- ComBat(dat = combined_mat_log_filtered, batch = batch_vector, mod = mod_combat, par.prior = TRUE, prior.plots = FALSE)
    print("ComBat correction applied successfully.")
  }, error = function(e) {
    print(paste("Error during ComBat:", e$message))
    print("Proceeding with uncorrected data for visualization.")
    combined_mat_corrected <<- combined_mat_log_filtered # Use <<- for assignment in parent environment
  })
}


# Z-score scaling for heatmap (on corrected data)
# Ensure combined_mat_corrected exists before proceeding
if (!exists("combined_mat_corrected") || is.null(combined_mat_corrected) || nrow(combined_mat_corrected) == 0) {
  stop("combined_mat_corrected is not available or empty. Cannot proceed with heatmap/PCA.")
}
combined_mat_z <- t(scale(t(combined_mat_corrected)))


dataset_colors <- c("Table_S1" = "lightblue", "Poland_Yeast" = "lightcoral")

unique_conditions <- unique(sample_annotation$Condition)
condition_colors_palette <- rainbow(length(unique_conditions))
condition_colors <- setNames(condition_colors_palette, unique_conditions)

# Generate colors for the new Experimental_Group annotation
unique_groups <- unique(sample_annotation$Experimental_Group)
print("Unique groups for coloring:")
print(unique_groups)
print("Number of unique groups:")
print(length(unique_groups))

# Use a more robust color palette
if(length(unique_groups) <= 8) {
  group_colors_palette <- RColorBrewer::brewer.pal(n = max(3, length(unique_groups)), name = "Set2")[seq_len(length(unique_groups))]
} else {
  group_colors_palette <- rainbow(length(unique_groups))
}
group_colors <- setNames(group_colors_palette, unique_groups)

print("Group colors:")
print(group_colors)

annotation_colors <- list(
  Dataset = dataset_colors,
  Condition = condition_colors,
  Experimental_Group = group_colors
)

# --- 7. Create Combined Heatmap (Using Corrected Data) ---
# Debug: Check annotation structure for pheatmap
print("Annotation for pheatmap:")
annotation_for_heatmap <- sample_annotation[, c("Dataset", "Condition", "Experimental_Group"), drop = FALSE]
print(annotation_for_heatmap)
print("Rownames of annotation:")
print(rownames(annotation_for_heatmap))
print("Colnames of matrix:")
print(colnames(combined_mat_z))

pheatmap(
  combined_mat_z, # Use Z-scored corrected data
  annotation_col = annotation_for_heatmap,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = ifelse(nrow(combined_mat_z) <= 50, TRUE, FALSE),
  show_colnames = TRUE,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = paste("Combined Heatmap: All Samples (ComBat Corrected, ", nrow(combined_mat_z), " genes)"),
  fontsize_col = 8,
  fontsize_row = 6,
  angle_col = 45,
  filename = file.path(output_dir, "combined_heatmap_combat_corrected.png"),
  width = 12,
  height = 10
)

# --- 8. Create Combined PCA (Using Corrected Data) ---
# Use the ComBat corrected matrix for PCA
# Ensure combined_mat_corrected exists before proceeding
if (!exists("combined_mat_corrected") || is.null(combined_mat_corrected) || nrow(combined_mat_corrected) == 0) {
  message("PCA skipped: combined_mat_corrected is not available or empty.")
  combined_mat_pca_input <- matrix(nrow=0, ncol=0) # Create an empty matrix to avoid further errors
} else {
  variances_corrected <- apply(combined_mat_corrected, 1, var, na.rm = TRUE)
  combined_mat_pca_input <- combined_mat_corrected[variances_corrected > 1e-6, , drop = FALSE] # Ensure variance for PCA
}

if(nrow(combined_mat_pca_input) >= 2 && ncol(combined_mat_pca_input) >= 2) {
  pca_combined <- prcomp(t(combined_mat_pca_input), scale. = TRUE)
  pca_df <- as.data.frame(pca_combined$x)
  pca_df$Sample <- rownames(pca_df)
  # Join with sample_annotation which has Batch and Condition info
  pca_df <- pca_df %>%
    left_join(sample_annotation, by = "Sample")

  var_explained <- summary(pca_combined)$importance[2, 1:2] * 100

  pca_plot_dataset <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Dataset, shape = Dataset)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = Sample), size = 3, max.overlaps = Inf) +
    labs(
      title = "PCA: Combined Data by Dataset (ComBat Corrected)",
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)")
    ) +
    scale_color_manual(values = dataset_colors) +
    scale_shape_manual(values = c("Table_S1" = 16, "Poland_Yeast" = 17)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
  print(pca_plot_dataset)
  ggsave(file.path(output_dir, "pca_by_dataset_combat.png"), plot = pca_plot_dataset, width = 10, height = 8, dpi = 300)
  ggsave(filename = file.path(output_dir, "PCA_combined_data_by_dataset.png"), plot = pca_plot_dataset, width = 10, height = 8, dpi = 300)

  pca_plot_condition <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = Sample), size = 3, max.overlaps = Inf) +
    labs(
      title = "PCA: Combined Data by Condition (ComBat Corrected)",
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)")
    ) +
    scale_color_manual(values = condition_colors) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"), legend.text = element_text(size = 8)) +
    guides(color = guide_legend(ncol = 1))
  print(pca_plot_condition)
  ggsave(file.path(output_dir, "pca_by_condition_combat.png"), plot = pca_plot_condition, width = 12, height = 8, dpi = 300)
  ggsave(filename = file.path(output_dir, "PCA_combined_data_by_condition.png"), plot = pca_plot_condition, width = 10, height = 8, dpi = 300)
  
  pca_plot_experimental_group <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Experimental_Group)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = Sample), size = 3, max.overlaps = Inf) +
    labs(
      title = "PCA: Combined Data by Experimental Group (ComBat Corrected)",
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)")
    ) +
    scale_color_manual(values = group_colors) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"), legend.text = element_text(size = 8)) +
    guides(color = guide_legend(ncol = 1))
  print(pca_plot_experimental_group)
  ggsave(file.path(output_dir, "pca_by_experimental_group_combat.png"), plot = pca_plot_experimental_group, width = 12, height = 8, dpi = 300)
  ggsave(filename = file.path(output_dir, "PCA_combined_data_by_experimental_group.png"), plot = pca_plot_experimental_group, width = 10, height = 8, dpi = 300)
} else {
  message("PCA skipped: Not enough features with non-zero variance or insufficient samples after ComBat correction and filtering.")
}

# --- 9. UMAP Analysis (Using Corrected Data) ---
if (exists("combined_mat_corrected") && !is.null(combined_mat_corrected) && nrow(combined_mat_corrected) > 0) {
  print("Performing UMAP analysis on ComBat corrected data...")

  # UMAP is typically run on samples (columns), so we transpose the matrix
  # Use the same input as PCA for consistency
  if (exists("combined_mat_pca_input") && nrow(combined_mat_pca_input) >= 2 && ncol(combined_mat_pca_input) >= 2) {
    umap_input <- t(combined_mat_pca_input)

    # UMAP requires n_neighbors to be smaller than the number of samples.
    # The default is 15, which can be too large for small datasets.
    n_samples <- nrow(umap_input)
    n_neighbors_val <- min(15, n_samples - 1)

    if (n_neighbors_val > 1) {
      print(paste("Running UMAP with n_neighbors =", n_neighbors_val))
      # Run UMAP
      set.seed(42) # for reproducibility
      umap_results <- umap(umap_input, n_neighbors = n_neighbors_val)

      # Create a data frame for plotting
      umap_df <- as.data.frame(umap_results$layout)
      colnames(umap_df) <- c("UMAP1", "UMAP2")
      umap_df$Sample <- rownames(umap_input)

      # Join with sample_annotation to get metadata
      umap_df <- umap_df %>%
        left_join(sample_annotation, by = "Sample")

      # Plot UMAP results, colored by Experimental_Group, with every gene labeled
      umap_plot_experimental_group <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Experimental_Group)) +
        geom_point(size = 4, alpha = 0.8) +
        geom_text_repel(aes(label = Sample), size = 3, max.overlaps = Inf, show.legend = FALSE) +
        labs(
          title = "UMAP: Combined Data by Experimental Group (ComBat Corrected)",
          x = "UMAP 1",
          y = "UMAP 2",
          color = "Experimental Group"
        ) +
        scale_color_manual(values = group_colors, labels = levels(sample_annotation$Experimental_Group)) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"), legend.text = element_text(size = 10)) +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))

      print(umap_plot_experimental_group)
      ggsave(file.path(output_dir, "umap_by_experimental_group_combat.png"), plot = umap_plot_experimental_group, width = 12, height = 8, dpi = 300)
      
      # Plot UMAP results, colored by Dataset, with every gene labeled
      umap_plot_dataset <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Dataset, shape = Dataset)) +
        geom_point(size = 4, alpha = 0.8) +
        geom_text_repel(aes(label = Sample), size = 3, max.overlaps = Inf, show.legend = FALSE) +
        labs(
          title = "UMAP: Combined Data by Dataset (ComBat Corrected)",
          x = "UMAP 1",
          y = "UMAP 2",
          color = "Dataset",
          shape = "Dataset"
        ) +
        scale_color_manual(values = dataset_colors, labels = names(dataset_colors)) +
        scale_shape_manual(values = c("Table_S1" = 16, "Poland_Yeast" = 17), labels = names(dataset_colors)) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"), legend.text = element_text(size = 10))

      print(umap_plot_dataset)
      ggsave(file.path(output_dir, "umap_by_dataset_combat.png"), plot = umap_plot_dataset, width = 10, height = 8, dpi = 300)

      # --- UMAP Group Description ---
      # Describe which groups are linked together in UMAP space
      # Calculate centroids for each Experimental_Group
      group_centroids <- umap_df %>%
        group_by(Experimental_Group) %>%
        summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = 'drop')
      print("UMAP group centroids (Experimental_Group):")
      print(group_centroids)
      # Calculate pairwise distances between centroids
      if(nrow(group_centroids) > 1) {
        dist_matrix <- as.matrix(dist(group_centroids[, c("UMAP1", "UMAP2")]))
        print("Pairwise distances between Experimental_Group centroids:")
        print(dist_matrix)
        # Find closest pairs
        for(i in 1:(nrow(dist_matrix)-1)) {
          closest <- which.min(dist_matrix[i, -i])
          group1 <- group_centroids$Experimental_Group[i]
          group2 <- group_centroids$Experimental_Group[-i][closest]
          dist_val <- dist_matrix[i, -i][closest]
          print(paste("Group", group1, "is closest to group", group2, "(distance:", round(dist_val, 2), ")"))
        }
      }
    } else {
      message("UMAP analysis skipped: Not enough samples to determine neighbors.")
    }
  } else {
    message("UMAP analysis skipped: Input matrix for PCA/UMAP is not valid.")
  }
} else {
  message("UMAP analysis skipped: 'combined_mat_corrected' is not available.")
}


# --- 10. Summary Statistics (Based on original cleaned data, or corrected if preferred) ---
# Summary stats are typically done on cleaned, uncorrected data or specifically on corrected data if effect of correction is being evaluated.
# Here, we'll use the original combined_data_clean for consistency with previous script,
# but you could adapt to use combined_mat_corrected if needed.
summary_stats <- combined_data_clean %>%
  pivot_longer(cols = -Gene_Symbol, names_to = "Sample", values_to = "Intensity") %>%
  left_join(sample_annotation, by = "Sample") %>% # sample_annotation now includes Batch and Experimental_Group
  group_by(Sample, Dataset, Condition, Batch, Experimental_Group) %>%
  summarise(
    Mean_Intensity = mean(Intensity, na.rm = TRUE),
    Median_Intensity = median(Intensity, na.rm = TRUE),
    SD_Intensity = sd(Intensity, na.rm = TRUE),
    N_Genes = n(),
    .groups = 'drop'
  )
print("Summary statistics by sample (based on original cleaned data):")
print(summary_stats)

# --- 11. Boxplot of Intensities by Dataset (Using Log2 Transformed Data before Correction) ---
# Boxplot of original log2 intensities to show pre-correction distributions
intensity_data_log_orig <- as.data.frame(combined_mat_log) %>% # Use combined_mat_log (pre-ComBat)
  tibble::rownames_to_column("Gene_Symbol") %>%
  pivot_longer(cols = -Gene_Symbol, names_to = "Sample", values_to = "Log2_Intensity") %>%
  left_join(sample_annotation, by = "Sample")

boxplot_intensities_before_correction <- ggplot(intensity_data_log_orig, aes(x = Sample, y = Log2_Intensity, fill = Dataset)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Distribution of Gene Expression (Log2) - Before ComBat",
    subtitle = paste("Based on", nrow(combined_mat_log), "common genes"),
    x = "Sample",
    y = "Log2(Intensity + 1)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
print(boxplot_intensities_before_correction)
ggsave(file.path(output_dir, "boxplot_intensities_before_combat.png"), plot = boxplot_intensities_before_correction, width = 12, height = 8, dpi = 300)
ggsave(filename = file.path(output_dir, "boxplot_intensities_before_correction.png"), plot = boxplot_intensities_before_correction, width = 10, height = 8, dpi = 300)

# Boxplot of ComBat corrected log2 intensities
if (exists("combined_mat_corrected") && !is.null(combined_mat_corrected) && nrow(combined_mat_corrected) > 0 && ncol(combined_mat_corrected) > 0) {
  intensity_data_log_corrected <- as.data.frame(combined_mat_corrected) %>%
    tibble::rownames_to_column("Gene_Symbol") %>%
    pivot_longer(cols = -Gene_Symbol, names_to = "Sample", values_to = "Log2_Intensity_Corrected") %>%
    left_join(sample_annotation, by = "Sample")

  boxplot_intensities_after_correction <- ggplot(intensity_data_log_corrected, aes(x = Sample, y = Log2_Intensity_Corrected, fill = Dataset)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = dataset_colors) +
    labs(
      title = "Distribution of Gene Expression (Log2) - After ComBat",
      subtitle = paste("Based on", nrow(combined_mat_corrected), "genes"),
      x = "Sample",
      y = "Log2 Adjusted Intensity"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  print(boxplot_intensities_after_correction)
  ggsave(file.path(output_dir, "boxplot_intensities_after_combat.png"), plot = boxplot_intensities_after_correction, width = 12, height = 8, dpi = 300)
  ggsave(filename = file.path(output_dir, "boxplot_intensities_after_correction.png"), plot = boxplot_intensities_after_correction, width = 10, height = 8, dpi = 300)
  
  # Additional boxplot colored by Experimental_Group
  boxplot_by_experimental_group <- ggplot(intensity_data_log_corrected, aes(x = Sample, y = Log2_Intensity_Corrected, fill = Experimental_Group)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = group_colors) +
    labs(
      title = "Distribution of Gene Expression by Experimental Group - After ComBat",
      subtitle = paste("Based on", nrow(combined_mat_corrected), "genes"),
      x = "Sample",
      y = "Log2 Adjusted Intensity"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  print(boxplot_by_experimental_group)
  ggsave(filename = file.path(output_dir, "boxplot_intensities_by_experimental_group.png"), plot = boxplot_by_experimental_group, width = 10, height = 8, dpi = 300)
}

# --- 12. Add CV Analysis for Table S1 Data ---
print("Generating CV analysis plots...")

# Calculate CV for Table S1 protein data (original replicates)
protein_cv_data <- protein_values_clean %>%
  select(Gene_Symbol, W303_s_1:BY4742_e_3) %>%
  pivot_longer(cols = -Gene_Symbol, names_to = "Sample", values_to = "Intensity") %>%
  separate(Sample, into = c("Strain", "Phase", "Replicate"), sep = "_") %>%
  group_by(Gene_Symbol, Strain, Phase) %>%
  summarise(
    Mean_Intensity = mean(Intensity, na.rm = TRUE),
    SD_Intensity = sd(Intensity, na.rm = TRUE),
    N_Replicates = n(),
    .groups = 'drop'
  ) %>%
  filter(N_Replicates >= 2 & Mean_Intensity > 0 & !is.na(SD_Intensity)) %>%
  mutate(CV_percent = (SD_Intensity / Mean_Intensity) * 100) %>%
  filter(!is.na(CV_percent) & is.finite(CV_percent))

# Create condition labels
protein_cv_data <- protein_cv_data %>%
  mutate(
    Condition_label = case_when(
      Phase == "s" ~ "stationary",
      Phase == "e" ~ "exponential",
      TRUE ~ Phase
    ),
    Group = paste(Strain, Condition_label, sep = "_")
  )

# Calculate median CV per group for plot annotation
protein_median_cv <- protein_cv_data %>%
  group_by(Strain, Condition_label, Group) %>%
  summarise(Median_CV = median(CV_percent, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Label = paste0("(", round(Median_CV, 1), ")")) %>%
  arrange(match(Group, c("W303_stationary", "BY4742_stationary", "W303_exponential", "BY4742_exponential")))

# Prepare data for plotting
protein_plot_data <- protein_cv_data %>%
  mutate(Group = factor(Group, levels = c("W303_stationary", "BY4742_stationary", "W303_exponential", "BY4742_exponential")))

# Create labels for x-axis
x_labels_protein <- setNames(
  paste0(protein_median_cv$Strain, "\n", str_to_title(protein_median_cv$Condition_label), "\n", protein_median_cv$Label),
  protein_median_cv$Group
)

# Get number of unique proteins plotted
n_proteins <- length(unique(protein_plot_data$Gene_Symbol))

# Create CV boxplot for proteins
cv_plot_proteins <- ggplot(protein_plot_data, aes(x = Group, y = CV_percent)) +
  geom_boxplot(outlier.shape = 8, outlier.size = 0.5, width = 0.6, fill = "white") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_x_discrete(labels = x_labels_protein) +
  labs(
    title = bquote(bold("C") ~ "  Proteins (n = " ~ .(format(n_proteins, big.mark = ",")) ~ ")"),
    x = NULL,
    y = "CV (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, size = 14),
    axis.text.x = element_text(size = 10, lineheight = 0.9),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(cv_plot_proteins)
ggsave(file.path(output_dir, "cv_analysis_proteins.png"), plot = cv_plot_proteins, width = 8, height = 6, dpi = 300)

# --- 13. CV Analysis for Poland Data ---
print("Performing CV analysis for Poland data...")

# Use raw Poland data and map to Gene Symbols
poland_cv_data_prep <- poland_data_raw %>%
  left_join(id_mapping, by = c("Accession" = "From")) %>%
  rename(Gene_Symbol = To) %>%
  filter(!is.na(Gene_Symbol) & Gene_Symbol != "") %>%
  select(-Accession) %>%
  # Handle multiple accessions for one gene symbol by averaging intensities per replicate
  group_by(Gene_Symbol) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')

poland_cv_data <- poland_cv_data_prep %>%
  # Pivot to long format
  pivot_longer(cols = -Gene_Symbol, names_to = "Sample", values_to = "Intensity") %>%
  # Extract condition and replicate info
  mutate(
    Group_Num = str_split_fixed(Sample, "_", 2)[,1],
    Replicate = str_split_fixed(Sample, "_", 2)[,2]
  ) %>%
  mutate(
    Condition = case_when(
      Group_Num == "1" ~ "WT_YPD_Glucose",
      Group_Num == "2" ~ "maf1D_YPD_Glucose",
      Group_Num == "3" ~ "rpc128_YPD_Glucose",
      Group_Num == "4" ~ "WT_YPGly_30C",
      Group_Num == "5" ~ "WT_YPGly_37C",
      Group_Num == "6" ~ "maf1D_YPGly_30C",
      Group_Num == "7" ~ "maf1D_YPGly_37C",
      Group_Num == "8" ~ "rpc128_YPGly_30C",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Calculate CV
  group_by(Gene_Symbol, Condition) %>%
  summarise(
    Mean_Intensity = mean(Intensity, na.rm = TRUE),
    SD_Intensity = sd(Intensity, na.rm = TRUE),
    N_Replicates = n(),
    .groups = 'drop'
  ) %>%
  filter(N_Replicates >= 2 & Mean_Intensity > 0 & !is.na(SD_Intensity)) %>%
  mutate(CV_percent = (SD_Intensity / Mean_Intensity) * 100) %>%
  filter(!is.na(CV_percent) & is.finite(CV_percent))

# Save CV results to a file
write.csv(poland_cv_data, file.path(output_dir, "poland_cv_analysis_results.csv"), row.names = FALSE)
print("Poland CV analysis results saved to poland_cv_analysis_results.csv")

# Calculate median CV per group for plot annotation
poland_median_cv <- poland_cv_data %>%
  group_by(Condition) %>%
  summarise(Median_CV = median(CV_percent, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Label = paste0("(", round(Median_CV, 1), "%)"))

# Create labels for x-axis, ensuring order is consistent
condition_order <- c(
  "WT_YPD_Glucose", "maf1D_YPD_Glucose", "rpc128_YPD_Glucose",
  "WT_YPGly_30C", "WT_YPGly_37C", "maf1D_YPGly_30C",
  "maf1D_YPGly_37C", "rpc128_YPGly_30C"
)
poland_cv_data$Condition <- factor(poland_cv_data$Condition, levels = condition_order)
poland_median_cv <- poland_median_cv %>%
  mutate(Condition = factor(Condition, levels = condition_order)) %>%
  arrange(Condition)

x_labels_poland <- setNames(
  paste0(gsub("_", " ", poland_median_cv$Condition), "\n", poland_median_cv$Label),
  poland_median_cv$Condition
)

# Get number of unique proteins plotted
n_proteins_poland <- length(unique(poland_cv_data$Gene_Symbol))

# Create CV boxplot for Poland data
cv_plot_poland <- ggplot(poland_cv_data, aes(x = Condition, y = CV_percent)) +
  geom_boxplot(outlier.shape = 8, outlier.size = 0.5, width = 0.7, fill = "lightblue") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_x_discrete(labels = x_labels_poland) +
  labs(
    title = bquote(bold("CV Analysis of Poland Dataset") ~ " (n = " ~ .(format(n_proteins_poland, big.mark = ",")) ~ " proteins)"),
    x = "Condition",
    y = "CV (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

print(cv_plot_poland)
ggsave(file.path(output_dir, "cv_analysis_poland.png"), plot = cv_plot_poland, width = 12, height = 8, dpi = 300)


# --- 14. Summary Report ---
print("\n=== ANALYSIS SUMMARY ===")
print(paste("Total common genes analyzed:", nrow(combined_data_clean)))
print(paste("Output directory:", output_dir))
print("Generated plots:")
plot_files <- list.files(output_dir, pattern = "\\.png$", full.names = FALSE)
for(file in plot_files) {
  print(paste("  -", file))
}
print("=== ANALYSIS COMPLETED ===\n")