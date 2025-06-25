library(readxl)       # version 1.4.5
library(readr)        # version 2.1.5
library(dplyr)        # version 1.1.4
library(tidyr)        # version 1.3.1
library(ggplot2)      # version 3.5.2
library(pheatmap)     # version 1.0.12
library(stringr)      # version 1.5.1
library(forcats)      # version 1.0.0
library(ggrepel)      # version 0.9.6
library(RColorBrewer) # version 1.1-3
library(limma)        # version 3.62.2

# Create output directory
output_dir <- "output_plots_protein_analysis2"
dir.create(output_dir, showWarnings = FALSE)

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
  main = "Table S1: Protein Heatmap (Mean Values)",
  filename = file.path(output_dir, "s1_heatmap.png"), width = 8, height = 6
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
s1_pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Label, color = Strain)) +
  geom_point(size = 3) +
  geom_text_repel(show.legend = FALSE) +
  labs(title = "PCA: Table S1 (Mean Values)", 
       x = paste0("PC1 (", round(100*summary(pca)$importance[2,1],1), "%)"),
       y = paste0("PC2 (", round(100*summary(pca)$importance[2,2],1), "%)")) +
  theme_minimal()
print(s1_pca_plot)
ggsave(file.path(output_dir, "s1_pca.png"), plot = s1_pca_plot, width = 8, height = 6)

# --- 4. Volcano Plot for Table S1 (Exponential phase) with t-test and BH correction ---

# Function to safely perform t-test and return p-value
# Welch's t-test is used by default in R's t.test, which is appropriate
# as we do not assume equal variance between groups.
t_test_safe <- function(x, y) {
  # Check for at least 2 valid values in each group
  if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2) {
    return(NA_real_)
  }
  # Perform t-test
  return(t.test(x, y)$p.value)
}

volcano_data_exp <- protein_values_clean %>%
  # Calculate p-value and log2FC for each protein
  rowwise() %>%
  mutate(
    p_value = t_test_safe(c(W303_e_1, W303_e_2, W303_e_3), c(BY4742_e_1, BY4742_e_2, BY4742_e_3)),
    log2FC = log2(mean(c(W303_e_1, W303_e_2, W303_e_3), na.rm = TRUE) / mean(c(BY4742_e_1, BY4742_e_2, BY4742_e_3), na.rm = TRUE))
  ) %>%
  ungroup() %>%
  # Filter out rows with NA p-values or infinite fold changes
  filter(!is.na(.data$p_value) & !is.infinite(.data$log2FC) & !is.na(.data$log2FC)) %>%
  mutate(
    adj_p_value = p.adjust(.data$p_value, method = "BH"),
    negLog10_adj_P = -log10(.data$adj_p_value),
    # Determine significance based on fold change and adjusted p-value
    sig = case_when(
      .data$log2FC > 1 & .data$adj_p_value < 0.05 ~ "up",
      .data$log2FC < -1 & .data$adj_p_value < 0.05 ~ "down",
      TRUE ~ "ns"
    )
  )

# Count significant proteins
n_up <- sum(volcano_data_exp$sig == "up", na.rm = TRUE)
n_down <- sum(volcano_data_exp$sig == "down", na.rm = TRUE)

# Identify top differentially expressed genes for labeling
top_n_genes <- 5

# Top N up-regulated (sorted by adjusted p-value then fold change)
top_up_proteins <- volcano_data_exp %>%
  filter(.data$sig == "up") %>%
  arrange(.data$adj_p_value, desc(abs(.data$log2FC))) %>%
  head(top_n_genes) %>%
  pull(.data$ProteinId)

# Top N down-regulated (sorted by adjusted p-value then fold change)
top_down_proteins <- volcano_data_exp %>%
  filter(.data$sig == "down") %>%
  arrange(.data$adj_p_value, desc(abs(.data$log2FC))) %>%
  head(top_n_genes) %>%
  pull(.data$ProteinId)

# Assign labels to the data frame
volcano_data_labeled <- volcano_data_exp %>%
  mutate(
    gene_label = ifelse(.data$ProteinId %in% c(top_up_proteins, top_down_proteins), as.character(.data$ProteinId), NA_character_)
  )

# Generate Volcano Plot
s1_volcano_exp_plot <- ggplot(volcano_data_labeled, aes(x = .data$log2FC, y = .data$negLog10_adj_P)) +
  geom_point(aes(color = .data$sig), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("up" = "lightblue", "down" = "khaki", "ns" = "grey"),
                     labels = c(up = "Up-regulated", down = "Down-regulated", ns = "Not significant")) +
  geom_text_repel(aes(label = .data$gene_label), na.rm = TRUE, 
                  box.padding = 0.5, max.overlaps = Inf, 
                  min.segment.length = 0, size = 3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  annotate("rect", xmin = -Inf, xmax = -1, ymin = -log10(0.05), ymax = Inf, 
           fill = "khaki", alpha = 0.1) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -log10(0.05), ymax = Inf, 
           fill = "lightblue", alpha = 0.1) +
  annotate("text", x = -4, y = max(volcano_data_labeled$negLog10_adj_P, na.rm = TRUE) * 0.95, 
           label = paste0("N=", n_down), hjust = 0.5, fontface = 2) +
  annotate("text", x = 4, y = max(volcano_data_labeled$negLog10_adj_P, na.rm = TRUE) * 0.95, 
           label = paste0("N=", n_up), hjust = 0.5, fontface = 2) +
  labs(title = "Volcano Plot (Exponential phase - W303 vs BY4742)", 
       x = "log2 Fold Change (W303 / BY4742)", 
       y = expression(-log[10]*" adjusted p-value"),
       color = "Regulation") +
  theme_minimal() +
  theme(legend.position = "top")
print(s1_volcano_exp_plot)
ggsave(file.path(output_dir, "s1_volcano_exp.png"), plot = s1_volcano_exp_plot, width = 8, height = 7)


# --- 5. Volcano Plot for Table S1 (Stationary phase) with t-test and BH correction ---
volcano_data_stationary <- protein_values_clean %>%
  rowwise() %>%
  mutate(
    p_value = t_test_safe(c(W303_s_1, W303_s_2, W303_s_3), c(BY4742_s_1, BY4742_s_2, BY4742_s_3)),
    log2FC = log2(mean(c(W303_s_1, W303_s_2, W303_s_3), na.rm = TRUE) / mean(c(BY4742_s_1, BY4742_s_2, BY4742_s_3), na.rm = TRUE))
  ) %>%
  ungroup() %>%
  filter(!is.na(.data$p_value) & !is.infinite(.data$log2FC) & !is.na(.data$log2FC)) %>%
  mutate(
    adj_p_value = p.adjust(.data$p_value, method = "BH"),
    negLog10_adj_P = -log10(.data$adj_p_value),
    sig = case_when(
      .data$log2FC > 1 & .data$adj_p_value < 0.05 ~ "up",
      .data$log2FC < -1 & .data$adj_p_value < 0.05 ~ "down",
      TRUE ~ "ns"
    )
  )

# Count significant
n_up_stat <- sum(volcano_data_stationary$sig == "up", na.rm = TRUE)
n_down_stat <- sum(volcano_data_stationary$sig == "down", na.rm = TRUE)

# Identify top differentially expressed genes for labeling
top_up_stationary <- volcano_data_stationary %>%
  filter(.data$sig == "up") %>%
  arrange(.data$adj_p_value, desc(abs(.data$log2FC))) %>%
  head(5) %>%
  pull(.data$ProteinId)

top_down_stationary <- volcano_data_stationary %>%
  filter(.data$sig == "down") %>%
  arrange(.data$adj_p_value, desc(abs(.data$log2FC))) %>%
  head(5) %>%
  pull(.data$ProteinId)

volcano_data_stationary_labeled <- volcano_data_stationary %>%
  mutate(
    gene_label = ifelse(.data$ProteinId %in% c(top_up_stationary, top_down_stationary), as.character(.data$ProteinId), NA_character_)
  )

s1_volcano_stat_plot <- ggplot(volcano_data_stationary_labeled, aes(x = .data$log2FC, y = .data$negLog10_adj_P)) +
  geom_point(aes(color = .data$sig), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("up" = "lightblue", "down" = "khaki", "ns" = "grey"),
                     labels = c(up = "Up-regulated", down = "Down-regulated", ns = "Not significant")) +
  geom_text_repel(aes(label = .data$gene_label), na.rm = TRUE, box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0, size = 3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  annotate("rect", xmin = -Inf, xmax = -1, ymin = -log10(0.05), ymax = Inf, fill = "khaki", alpha = 0.1) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -log10(0.05), ymax = Inf, fill = "lightblue", alpha = 0.1) +
  annotate("text", x = -2.5, y = max(volcano_data_stationary_labeled$negLog10_adj_P, na.rm = TRUE) * 0.95, 
           label = paste0("N=", n_down_stat), hjust = 0.5, fontface = 2) +
  annotate("text", x = 2.5, y = max(volcano_data_stationary_labeled$negLog10_adj_P, na.rm = TRUE) * 0.95, 
           label = paste0("N=", n_up_stat), hjust = 0.5, fontface = 2) +
  labs(title = "Volcano Plot (Stationary phase - W303 vs BY4742)", 
       x = "log2 Fold Change (W303 / BY4742)", 
       y = expression(-log[10]*" adjusted p-value"), 
       color = "Regulation") +
  theme_minimal() +
  theme(legend.position = "top")
print(s1_volcano_stat_plot)
ggsave(file.path(output_dir, "s1_volcano_stat.png"), plot = s1_volcano_stat_plot, width = 8, height = 7)


# --- 6. PolandYeast_Emma.csv Analysis (using means) ---
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
  main = "PolandYeast: Protein Heatmap (Mean Values)",
  filename = file.path(output_dir, "poland_heatmap.png"), width = 8, height = 6
)

# PCA for PolandYeast (using means)
# Use the same matrix as heatmap (already averaged)
pca_poland <- prcomp(t(poland_mat), scale. = TRUE)
pca_poland_df <- as.data.frame(pca_poland$x)
pca_poland_df$Sample <- rownames(pca_poland_df)
poland_pca_plot <- ggplot(pca_poland_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(show.legend = FALSE) +
  labs(title = "PCA: PolandYeast (Mean Values)", 
       x = paste0("PC1 (", round(100*summary(pca_poland)$importance[2,1],1), "%)"),
       y = paste0("PC2 (", round(100*summary(pca_poland)$importance[2,2],1), "%)")) +
  theme_minimal()
print(poland_pca_plot)
ggsave(file.path(output_dir, "poland_pca.png"), plot = poland_pca_plot, width = 8, height = 6)

# --- 7. Volcano Plots for Polish Data with t-test and BH correction ---
poland_volcano_ttest <- function(data, group1_cols, group2_cols, label1, label2, filename) {
  
  volcano_data <- data %>%
    rowwise() %>%
    mutate(
      p_value = t_test_safe(c_across(all_of(group1_cols)), c_across(all_of(group2_cols))),
      log2FC = log2(mean(c_across(all_of(group1_cols)), na.rm = TRUE) / mean(c_across(all_of(group2_cols)), na.rm = TRUE))
    ) %>%
    ungroup() %>%
    filter(!is.na(.data$p_value) & !is.infinite(.data$log2FC) & !is.na(.data$log2FC)) %>%
    mutate(
      adj_p_value = p.adjust(.data$p_value, method = "BH"),
      negLog10_adj_P = -log10(.data$adj_p_value),
      sig = case_when(
        .data$log2FC > 1 & .data$adj_p_value < 0.05 ~ "up",
        .data$log2FC < -1 & .data$adj_p_value < 0.05 ~ "down",
        TRUE ~ "ns"
      )
    )

  top_up <- volcano_data %>%
    filter(.data$sig == "up") %>%
    arrange(.data$adj_p_value, desc(abs(.data$log2FC))) %>%
    head(5) %>%
    pull(.data$Uniprot_ID)

  top_down <- volcano_data %>%
    filter(.data$sig == "down") %>%
    arrange(.data$adj_p_value, desc(abs(.data$log2FC))) %>%
    head(5) %>%
    pull(.data$Uniprot_ID)

  volcano_data_labeled <- volcano_data %>%
    mutate(
      gene_label = ifelse(.data$Uniprot_ID %in% c(top_up, top_down), as.character(.data$Uniprot_ID), NA_character_)
    )

  n_up <- sum(volcano_data_labeled$sig == "up")
  n_down <- sum(volcano_data_labeled$sig == "down")

  p <- ggplot(volcano_data_labeled, aes(x = .data$log2FC, y = .data$negLog10_adj_P)) +
    geom_point(aes(color = .data$sig), alpha = 0.7, size = 2) +
    scale_color_manual(values = c("up" = "lightblue", "down" = "khaki", "ns" = "grey"),
                       labels = c(up = "Up-regulated", down = "Down-regulated", ns = "Not significant")) +
    geom_text_repel(aes(label = .data$gene_label), na.rm = TRUE, box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0, size = 3) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    annotate("rect", xmin = -Inf, xmax = -1, ymin = -log10(0.05), ymax = Inf, fill = "khaki", alpha = 0.1) +
    annotate("rect", xmin = 1, xmax = Inf, ymin = -log10(0.05), ymax = Inf, fill = "lightblue", alpha = 0.1) +
    annotate("text", x = min(volcano_data_labeled$log2FC, na.rm=T) * 0.8, y = max(volcano_data_labeled$negLog10_adj_P, na.rm = TRUE) * 0.95, 
           label = paste0("N=", n_down), hjust = 0.5, fontface = 2) +
    annotate("text", x = max(volcano_data_labeled$log2FC, na.rm=T) * 0.8, y = max(volcano_data_labeled$negLog10_adj_P, na.rm = TRUE) * 0.95, 
           label = paste0("N=", n_up), hjust = 0.5, fontface = 2) +
    labs(title = paste("Volcano Plot:", label1, "vs", label2), x = paste("log2 Fold Change (", label1, "/", label2, ")"), y = expression(-log[10]*" adjusted p-value"), color = "Regulation") +
    theme_minimal() +
    theme(legend.position = "top")
  print(p)
  ggsave(filename, plot = p, width = 8, height = 7)
}

# Define column groups for comparisons
g1 <- c("1_1", "1_2", "1_3") # WT_YPD_Glucose
g2 <- c("2_1", "2_2", "2_3") # maf1D_YPD_Glucose
g3 <- c("3_1", "3_2", "3_3") # rpc128_YPD_Glucose
g4 <- c("4_1", "4_2", "4_3") # WT_YPGly_30C
g5 <- c("5_1", "5_2", "5_3") # WT_YPGly_37C
g6 <- c("6_1", "6_2", "6_3") # maf1D_YPGly_30C
g7 <- c("7_1", "7_2", "7_3") # maf1D_YPGly_37C
g8 <- c("8_1", "8_2", "8_3") # rpc128_YPGly_30C

# Generate plots for specified comparisons using t-test
poland_volcano_ttest(poland_data, g1, g2, "WT YPD Glucose", "maf1-D YPD Glucose", file.path(output_dir, "poland_volcano_wt_vs_maf1_ypd.png"))
poland_volcano_ttest(poland_data, g1, g3, "WT YPD Glucose", "rpc128 YPD Glucose", file.path(output_dir, "poland_volcano_wt_vs_rpc128_ypd.png"))
poland_volcano_ttest(poland_data, g4, g6, "WT YPGly 30C", "maf1-D YPGly 30C", file.path(output_dir, "poland_volcano_wt_vs_maf1_ypgly30.png"))
poland_volcano_ttest(poland_data, g4, g8, "WT YPGly 30C", "rpc128 YPGly 30C", file.path(output_dir, "poland_volcano_wt_vs_rpc128_ypgly30.png"))
poland_volcano_ttest(poland_data, g5, g7, "WT YPGly 37C", "maf1-D YPGly 37C", file.path(output_dir, "poland_volcano_wt_vs_maf1_ypgly37.png"))

# --- END ---