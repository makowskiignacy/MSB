# Load libraries
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(stringr)
library(forcats)
library(ggrepel)
library(RColorBrewer)

# File paths
protein_file <- "Table_S1_wt.xlsx"
poland_file <- "PolandYeast_Emma.csv"

# --- 1. Load and Prepare Table S1 Data ---
protein_data <- read_excel(protein_file, skip = 1)
protein_values <- protein_data %>% select(1, 5:16)
colnames(protein_values) <- c(
  "ProteinId", "W303_s_1", "W303_s_2", "W303_s_3",
  "BY4742_s_1", "BY4742_s_2", "BY4742_s_3",
  "W303_e_1", "W303_e_2", "W303_e_3",
  "BY4742_e_1", "BY4742_e_2", "BY4742_e_3"
)

# Clean and convert to numeric
protein_values_clean <- protein_values %>%
  filter(!is.na(ProteinId)) %>%
  filter(!ProteinId %in% c("", "Protein Id"))

# Convert quantification columns to numeric
for (col in 2:ncol(protein_values_clean)) {
  protein_values_clean[[col]] <- as.numeric(as.character(protein_values_clean[[col]]))
}

# Calculate mean values for each condition (combining 3 replicates)
protein_means <- protein_values_clean %>%
  mutate(
    W303_s = rowMeans(select(., W303_s_1, W303_s_2, W303_s_3), na.rm = TRUE),
    BY4742_s = rowMeans(select(., BY4742_s_1, BY4742_s_2, BY4742_s_3), na.rm = TRUE),
    W303_e = rowMeans(select(., W303_e_1, W303_e_2, W303_e_3), na.rm = TRUE),
    BY4742_e = rowMeans(select(., BY4742_e_1, BY4742_e_2, BY4742_e_3), na.rm = TRUE)
  ) %>%
  select(ProteinId, W303_s, BY4742_s, W303_e, BY4742_e)

# --- 2. Heatmap for Table S1 (using means) ---
mat <- as.matrix(protein_means[,-1])
# Remove rows with any NA
good_rows <- complete.cases(mat)
mat <- mat[good_rows, ]
rownames(mat) <- protein_means$ProteinId[good_rows]
# Z-score by row
mat_z <- t(scale(t(mat)))
pheatmap(
  mat_z,
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = FALSE, show_colnames = TRUE,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "Table S1: Protein Heatmap (Mean Values)"
)

# --- 3. PCA for Table S1 (using means) ---
# Use the same matrix as heatmap (already averaged)
pca <- prcomp(t(mat), scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$Sample <- rownames(pca_df)
pca_df <- pca_df %>%
  mutate(
    Strain = str_extract(Sample, "W303|BY4742"),
    Condition = ifelse(str_detect(Sample, "_s$"), "s", "e"),
    Label = paste0(Strain, Condition)
  )
ggplot(pca_df, aes(x = PC1, y = PC2, label = Label, color = Strain)) +
  geom_point(size = 3) +
  geom_text_repel(show.legend = FALSE) +
  labs(title = "PCA: Table S1 (Mean Values)", 
       x = paste0("PC1 (", round(100*summary(pca)$importance[2,1],1), "%)"),
       y = paste0("PC2 (", round(100*summary(pca)$importance[2,2],1), "%)")) +
  theme_minimal()

# --- 4. Volcano Plot for Table S1 (Exponential phase means) ---
volcano_data <- protein_means %>%
  filter(complete.cases(W303_e, BY4742_e)) %>%
  mutate(
    log2FC = log2(W303_e / BY4742_e),
    # For means, we can't calculate proper p-values without replicates
    # So we'll use fold change thresholds only
    sig = case_when(
      log2FC > 1 ~ "up",
      log2FC < -1 ~ "down",
      TRUE ~ "ns"
    )
  ) %>%
  filter(!is.infinite(log2FC) & !is.na(log2FC))

# Count significant
n_up <- sum(volcano_data$sig == "up", na.rm = TRUE)
n_down <- sum(volcano_data$sig == "down", na.rm = TRUE)

# Create a pseudo p-value based on fold change magnitude for visualization
volcano_data <- volcano_data %>%
  mutate(
    pseudo_pval = 1 / (abs(log2FC) + 0.1),  # Higher fold change = lower pseudo p-value
    negLog10P = -log10(pseudo_pval)
  )

# Identify top differentially expressed genes for labeling
top_n_genes <- 5
volcano_data_labeled <- volcano_data %>%
  mutate(label = NA_character_) # Initialize label column

# Top N up-regulated
top_up_proteins <- volcano_data %>%
  filter(sig == "up") %>%
  arrange(desc(log2FC)) %>%
  head(top_n_genes) %>%
  pull(ProteinId)

# Top N down-regulated
top_down_proteins <- volcano_data %>%
  filter(sig == "down") %>%
  arrange(log2FC) %>%
  head(top_n_genes) %>%
  pull(ProteinId)

# Assign labels to the data frame
volcano_data_labeled <- volcano_data %>%
  mutate(
    gene_label = ifelse(ProteinId %in% c(top_up_proteins, top_down_proteins), as.character(ProteinId), NA_character_)
  )


ggplot(volcano_data_labeled, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) + # Added color based on significance
  scale_color_manual(values = c("up" = "lightblue", "down" = "khaki", "ns" = "grey")) + # Define colors
  geom_text_repel(aes(label = gene_label), na.rm = TRUE, 
                  box.padding = 0.5, max.overlaps = Inf, 
                  min.segment.length = 0, size = 3) + # Added gene labels
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  annotate("rect", xmin = -Inf, xmax = -1, ymin = -Inf, ymax = Inf, 
           fill = "khaki", alpha = 0.1) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "lightblue", alpha = 0.1) +
  annotate("text", x = -4, y = max(volcano_data_labeled$negLog10P, na.rm = TRUE) * 0.95, 
           label = paste0("N=", n_down), hjust = 0.5, fontface = 2) +
  annotate("text", x = 4, y = max(volcano_data_labeled$negLog10P, na.rm = TRUE) * 0.95, 
           label = paste0("N=", n_up), hjust = 0.5, fontface = 2) +
  labs(title = "Volcano Plot (Exponential phase - Mean Values)", 
       x = "log2 W303/BY4742", 
       y = expression(-log[10]*" pseudo p-value"),
       color = "Regulation") + # Added legend title for color
  theme_minimal() +
  theme(legend.position = "top") # Moved legend to top


# --- 5. PolandYeast_Emma.csv Analysis (using means) ---
poland_data <- read_csv(poland_file, show_col_types = FALSE)

# Calculate means for each group with descriptive names
poland_means <- poland_data %>%
  select(Uniprot_ID, `1_1`:`8_3`) %>%
  mutate(
    `WT_YPD_Glucose` = rowMeans(select(., `1_1`, `1_2`, `1_3`), na.rm = TRUE),
    `maf1D_YPD_Glucose` = rowMeans(select(., `2_1`, `2_2`, `2_3`), na.rm = TRUE),
    `rpc128_YPD_Glucose` = rowMeans(select(., `3_1`, `3_2`, `3_3`), na.rm = TRUE),
    `WT_YPGly_30C` = rowMeans(select(., `4_1`, `4_2`, `4_3`), na.rm = TRUE),
    `WT_YPGly_37C` = rowMeans(select(., `5_1`, `5_2`, `5_3`), na.rm = TRUE),
    `maf1D_YPGly_30C` = rowMeans(select(., `6_1`, `6_2`, `6_3`), na.rm = TRUE),
    `maf1D_YPGly_37C` = rowMeans(select(., `7_1`, `7_2`, `7_3`), na.rm = TRUE),
    `rpc128_YPGly_30C` = rowMeans(select(., `8_1`, `8_2`, `8_3`), na.rm = TRUE)
  ) %>%
  select(Uniprot_ID, `WT_YPD_Glucose`:`rpc128_YPGly_30C`)

# Heatmap for PolandYeast (using means)
poland_mat <- as.matrix(poland_means[,-1])
poland_good_rows <- complete.cases(poland_mat)
poland_mat <- poland_mat[poland_good_rows, ]
rownames(poland_mat) <- poland_means$Uniprot_ID[poland_good_rows]

pheatmap(
  t(scale(t(poland_mat))),
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = FALSE, show_colnames = TRUE,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "PolandYeast: Protein Heatmap (Mean Values)"
)

# PCA for PolandYeast (using means)
# Use the same matrix as heatmap (already averaged)
pca_poland <- prcomp(t(poland_mat), scale. = TRUE)
pca_poland_df <- as.data.frame(pca_poland$x)
pca_poland_df$Sample <- rownames(pca_poland_df)
ggplot(pca_poland_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(show.legend = FALSE) +
  labs(title = "PCA: PolandYeast (Mean Values)", 
       x = paste0("PC1 (", round(100*summary(pca_poland)$importance[2,1],1), "%)"),
       y = paste0("PC2 (", round(100*summary(pca_poland)$importance[2,2],1), "%)")) +
  theme_minimal()

# --- END ---