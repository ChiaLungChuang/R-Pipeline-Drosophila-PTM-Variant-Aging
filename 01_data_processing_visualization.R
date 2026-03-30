# ============================================================
# Project: PTM Proteomics — Drosophila Variant × Aging Study
# Script:  01_data_processing_visualization.R
# Purpose: End-to-end pipeline for PTM-focused proteomics
#          data. Two experimental dimensions:
#            (A) Age comparison: Group 3 vs Group 1
#            (B) Genotype comparison: Variant vs Control
#          Steps:
#            1. Load raw data, remove contaminants
#            2. Organize into long format (6 comparisons)
#            3. Volcano plots per comparison
#            4. Global PTM pattern barplot and heatmap
#            5. Overview UpSet plot (all PTM types)
#            6. Per-PTM-type UpSet plots with site-level
#               gene labels
#            7. Publication-quality scatter plots:
#               Fraction A FC vs Fraction B FC (age)
#               Fraction A FC vs Fraction B FC (genotype)
#               Fraction A FC vs fraction ratio FC (age)
# Input:   dataA_ratio_comparisons.xlsx
#            (fraction ratio limma results — 2 sheets)
#          dataA_levels_comparisons.xlsx
#            (fraction-level limma results — 4 sheets)
# Output:  organized_PTM_data.csv
#          clean_data/*.csv (6 files)
#          volcano_*.pdf
#          global_PTM_patterns_*.pdf
#          upset_*.pdf / upset_*_summary.csv
#          scatter_*.pdf/.png + _data.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
#          Organism: Drosophila melanogaster.
#          Variant identity and specific PTM types have been
#          anonymized pending publication.
#          Sheet names in input Excel files must match those
#          listed in Step 1 below; update if your file uses
#          different sheet naming conventions.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggrepel)
  library(ComplexHeatmap)
  library(UpSetR)
  library(RColorBrewer)
})

# NOTE: Set your working directory to the folder containing input files
# setwd("path/to/your/data")

# ============================================================
# Modification Type Color and Label Scheme
# NOTE: Update keys to match PTM.types values in your data
# ============================================================

ptm_colors <- c(
  mod_type_1 = "#FF0000",
  mod_type_2 = "#0000FF",
  mod_type_3 = "#00FF00",
  mod_type_4 = "#FF00FF",
  mod_type_5 = "#FFD700",
  mod_type_6 = "#00FFFF",
  mod_type_7 = "#FF8C00",
  mod_type_8 = "#8B008B"
)

ptm_labels <- c(
  mod_type_1 = "Modification Type 1",
  mod_type_2 = "Modification Type 2",
  mod_type_3 = "Modification Type 3",
  mod_type_4 = "Modification Type 4",
  mod_type_5 = "Modification Type 5",
  mod_type_6 = "Modification Type 6",
  mod_type_7 = "Modification Type 7",
  mod_type_8 = "Modification Type 8"
)

# Comparison ordering for plots (update if comparison names change)
comparison_order <- c(
  "Variant_Group3_vs_Group1_Ratio",
  "Variant_vs_Control_Group3_Ratio",
  "Variant_Group3_vs_Group1_SOL",
  "Variant_vs_Control_Group3_SOL",
  "Variant_Group3_vs_Group1_INSOL",
  "Variant_vs_Control_Group3_INSOL"
)

# ============================================================
# 1. Load Data and Remove Contaminants
# ============================================================

cat("Step 1: Loading data and removing contaminants...\n")

if (dir.exists("clean_data") && file.exists("clean_data/ratio_age_clean.csv")) {
  
  cat("  Using existing clean data\n")
  ratio_age       <- read.csv("clean_data/ratio_age_clean.csv")
  ratio_geno      <- read.csv("clean_data/ratio_geno_clean.csv")
  levels_insol_age  <- read.csv("clean_data/levels_insol_age_clean.csv")
  levels_insol_geno <- read.csv("clean_data/levels_insol_geno_clean.csv")
  levels_sol_age    <- read.csv("clean_data/levels_sol_age_clean.csv")
  levels_sol_geno   <- read.csv("clean_data/levels_sol_geno_clean.csv")
  
} else {
  
  cat("  Loading raw data...\n")
  
  # NOTE: Sheet names below must match your actual Excel files.
  #       Update sheet names if your file uses different conventions.
  ratio_age       <- read_excel("dataA_ratio_comparisons.xlsx",
                                sheet = "Variant_Group3_vs_Group1")
  ratio_geno      <- read_excel("dataA_ratio_comparisons.xlsx",
                                sheet = "Variant_vs_Control_Group3")
  levels_insol_age  <- read_excel("dataA_levels_comparisons.xlsx",
                                  sheet = "INSOL_Variant_Group3_vs_Group1")
  levels_insol_geno <- read_excel("dataA_levels_comparisons.xlsx",
                                  sheet = "INSOL_Variant_vs_Control_Group3")
  levels_sol_age    <- read_excel("dataA_levels_comparisons.xlsx",
                                  sheet = "SOL_Variant_Group3_vs_Group1")
  levels_sol_geno   <- read_excel("dataA_levels_comparisons.xlsx",
                                  sheet = "SOL_Variant_vs_Control_Group3")
  
  # Legitimate Drosophila proteins to retain even if flagged
  legitimate_drosophila <- c("Act5C","Gapdh1","Gapdh2","Pgi",
                             "Aldolase","Enolase","Pgk")
  
  # Contaminant patterns (common cRAP and non-Drosophila proteins)
  contaminant_pattern <- paste(c(
    "KRT","Keratin","Krt","keratin",
    "ALB","Albumin","Serum albumin",
    "TRYP","Trypsin","PRSS",
    "IGG","IgG","IGH","IGK","IGL","Immunoglobulin",
    "ALBU_BOVIN","CASB_BOVIN","Casein",
    "Hemoglobin","Myoglobin","Fibrinogen"
  ), collapse = "|")
  
  flag_contaminants <- function(data) {
    data %>%
      mutate(
        is_contaminant =
          (grepl("^co\\|", `Protein ID`) |
             grepl(contaminant_pattern, GN, ignore.case = TRUE) |
             grepl(contaminant_pattern, `Protein ID`, ignore.case = TRUE)) &
          !GN %in% legitimate_drosophila &
          !word(GN, 1) %in% legitimate_drosophila
      )
  }
  
  datasets <- list(ratio_age, ratio_geno, levels_insol_age,
                   levels_insol_geno, levels_sol_age, levels_sol_geno)
  datasets <- lapply(datasets, flag_contaminants)
  
  n_contam <- sum(datasets[[1]]$is_contaminant)
  cat("  Contaminants found:", n_contam,
      "(", round(n_contam / nrow(datasets[[1]]) * 100, 2), "%)\n")
  
  datasets <- lapply(datasets, function(d) d %>%
                       filter(!is_contaminant) %>%
                       select(-is_contaminant))
  
  list2env(setNames(datasets,
                    c("ratio_age","ratio_geno","levels_insol_age",
                      "levels_insol_geno","levels_sol_age","levels_sol_geno")),
           envir = .GlobalEnv)
  
  dir.create("clean_data", showWarnings = FALSE)
  write.csv(ratio_age,        "clean_data/ratio_age_clean.csv",        row.names = FALSE)
  write.csv(ratio_geno,       "clean_data/ratio_geno_clean.csv",       row.names = FALSE)
  write.csv(levels_insol_age, "clean_data/levels_insol_age_clean.csv", row.names = FALSE)
  write.csv(levels_insol_geno,"clean_data/levels_insol_geno_clean.csv",row.names = FALSE)
  write.csv(levels_sol_age,   "clean_data/levels_sol_age_clean.csv",   row.names = FALSE)
  write.csv(levels_sol_geno,  "clean_data/levels_sol_geno_clean.csv",  row.names = FALSE)
  cat("  Contaminants removed; clean data cached in clean_data/\n")
}

cat("Step 1 complete.\n\n")

# ============================================================
# 2. Organize into Long Format
# ============================================================

cat("Step 2: Organizing data into long format...\n")

make_long <- function(data, comparison_name, fraction_type) {
  data %>%
    select(
      Peptides, GN, `Protein ID`, `PTM types`, ModSites,
      logFC       = starts_with("logFC"),
      adj.P.Val   = starts_with("adj.P.Val"),
      P.Value     = starts_with("P.Value"),
      AveExpr, t, B
    ) %>%
    mutate(
      Comparison = comparison_name,
      Fraction   = fraction_type,
      Identifier = paste(GN, `Protein ID`, `PTM types`, ModSites, sep = "|")
    )
}

data_long <- bind_rows(
  make_long(ratio_age,        "Variant_Group3_vs_Group1_Ratio",    "Ratio"),
  make_long(ratio_geno,       "Variant_vs_Control_Group3_Ratio",   "Ratio"),
  make_long(levels_insol_age, "Variant_Group3_vs_Group1_INSOL",    "INSOL"),
  make_long(levels_insol_geno,"Variant_vs_Control_Group3_INSOL",   "INSOL"),
  make_long(levels_sol_age,   "Variant_Group3_vs_Group1_SOL",      "SOL"),
  make_long(levels_sol_geno,  "Variant_vs_Control_Group3_SOL",     "SOL")
)

write.csv(data_long, "organized_PTM_data.csv", row.names = FALSE)

cat("  Total entries:", nrow(data_long), "\n")
cat("  Unique PTM sites:", n_distinct(data_long$Identifier), "\n")
cat("  Unique proteins:", n_distinct(data_long$GN), "\n")
cat("Step 2 complete.\n\n")

# ============================================================
# 3. Volcano Plots
# ============================================================

cat("Step 3: Creating volcano plots...\n")

create_volcano <- function(data, comparison_name, label_top_n = 20,
                           fc_thr = 1, pval_thr = 0.05,
                           labeled = TRUE) {
  
  plot_data <- data %>%
    filter(Comparison == comparison_name) %>%
    mutate(
      significant = adj.P.Val < pval_thr & abs(logFC) > fc_thr,
      PTM_types   = `PTM types`
    )
  
  if (labeled) {
    top_ids <- plot_data %>%
      filter(significant) %>%
      arrange(adj.P.Val, desc(abs(logFC))) %>%
      head(label_top_n) %>%
      pull(Identifier)
    plot_data <- plot_data %>%
      mutate(label = ifelse(Identifier %in% top_ids & significant,
                            paste0(GN, "\n", ModSites), ""))
  } else {
    plot_data <- mutate(plot_data, label = "")
  }
  
  p <- ggplot(plot_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = PTM_types, alpha = significant,
                   shape = significant), size = 2.5) +
    scale_color_manual(values = ptm_colors, labels = ptm_labels,
                       name = "Modification Type") +
    scale_alpha_manual(values = c("TRUE" = 0.6, "FALSE" = 0.3),
                       name = "Significance") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       name = "Significance") +
    geom_vline(xintercept = c(-fc_thr, fc_thr),
               linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(pval_thr),
               linetype = "dashed", color = "gray50") +
    labs(title    = gsub("_", " ", comparison_name),
         x = "log2 Fold Change",
         y = "-log10(adjusted p-value)",
         caption  = paste0("Significance: adj.p < ", pval_thr,
                           " & |log2FC| > ", fc_thr)) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          legend.title  = element_text(face = "bold"))
  
  if (labeled && any(plot_data$label != ""))
    p <- p + geom_text_repel(
      data = filter(plot_data, label != ""),
      aes(label = label), size = 2.5, max.overlaps = Inf,
      force = 10, box.padding = 0.3, segment.color = "grey50")
  
  p
}

for (comp in unique(data_long$Comparison)) {
  safe <- gsub("[^A-Za-z0-9]","_", comp)
  ggsave(paste0("volcano_", safe, "_labeled.pdf"),
         create_volcano(data_long, comp, labeled = TRUE),
         width = 10, height = 8, dpi = 600)
  ggsave(paste0("volcano_", safe, "_nolabel.pdf"),
         create_volcano(data_long, comp, labeled = FALSE),
         width = 10, height = 8, dpi = 600)
}

cat("  Volcano plots saved.\n")
cat("Step 3 complete.\n\n")

# ============================================================
# 4. Global PTM Pattern Plots
# ============================================================

cat("Step 4: Global PTM pattern plots...\n")

make_summary <- function(data, sig_thr = 0.05) {
  data %>%
    group_by(`PTM types`, Comparison) %>%
    summarise(
      n_sig       = sum(adj.P.Val < sig_thr, na.rm = TRUE),
      total_sites = n(),
      pct_changed = n_sig / total_sites * 100,
      mean_logFC  = mean(logFC[adj.P.Val < sig_thr], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      PTM_label      = ptm_labels[`PTM types`],
      Comparison     = factor(Comparison, levels = comparison_order),
      Comparison_lab = factor(gsub("_"," ", Comparison),
                              levels = gsub("_"," ", comparison_order))
    )
}

summary_data <- make_summary(data_long)

# Barplot
p_bar <- ggplot(summary_data,
                aes(x = Comparison_lab, y = pct_changed)) +
  geom_col(aes(fill = mean_logFC), width = 0.7,
           color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", pct_changed)),
            vjust = -0.5, size = 2.5) +
  geom_text(aes(label = sprintf("(n=%d)", n_sig)),
            vjust = -2, size = 2.5) +
  facet_wrap(~ PTM_label, scales = "free_y", ncol = 2) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                       midpoint = 0, name = "Mean log2FC") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "Global Modification Regulation Patterns",
       x = "", y = "% Sites Changed") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
        strip.text   = element_text(face = "bold"))

ggsave("global_PTM_patterns_barplot.pdf", p_bar,
       width = 12, height = 10, dpi = 600)

# Heatmap
p_heat <- ggplot(summary_data,
                 aes(x = Comparison_lab, y = PTM_label)) +
  geom_tile(aes(fill = pct_changed), color = "white") +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct_changed, n_sig)),
            size = 3) +
  scale_fill_gradient(low = "#4575B4", high = "#D73027",
                      name = "% Changed") +
  labs(title = "Global Modification Regulation Heatmap",
       x = "", y = "Modification Type") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank())

ggsave("global_PTM_patterns_heatmap.pdf", p_heat,
       width = 10, height = 6, dpi = 600)

cat("  Global pattern plots saved.\n")
cat("Step 4 complete.\n\n")

# ============================================================
# 5. Overview UpSet Plot (All Modification Types)
# ============================================================

cat("Step 5: Overview UpSet plot...\n")

upset_wide <- data_long %>%
  filter(adj.P.Val < 0.05) %>%
  select(Identifier, GN, ModSites, `PTM types`, Comparison) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(id_cols = c(Identifier, GN, ModSites, `PTM types`),
              names_from = Comparison, values_from = present,
              values_fill = 0)

mat_overview <- as.data.frame(upset_wide[, 5:ncol(upset_wide)])
rownames(mat_overview) <- upset_wide$Identifier

pdf("overview_upset_plot.pdf", width = 14, height = 8)
upset(mat_overview, nsets = ncol(mat_overview), nintersects = 40,
      order.by = "freq",
      main.bar.color = "darkred", sets.bar.color = "darkblue")
dev.off()

# Export summary
overview_summary <- upset_wide %>%
  mutate(
    n_comparisons   = rowSums(select(., starts_with("Variant"))),
    comparisons_list = apply(select(., starts_with("Variant")), 1,
                             function(x) paste(names(x)[x == 1],
                                               collapse = "; "))
  ) %>%
  select(GN, ModSites, `PTM types`, n_comparisons, comparisons_list) %>%
  arrange(desc(n_comparisons))

write.csv(overview_summary, "overview_upset_summary.csv", row.names = FALSE)
cat("  Overview UpSet plot saved.\n")
cat("Step 5 complete.\n\n")

# ============================================================
# 6. Per-Modification-Type UpSet Plots
# ============================================================

cat("Step 6: Per-modification-type UpSet plots...\n")

create_mod_upset <- function(data, mod_type, sig_thr = 0.05,
                             pt_size = 4, row_gap = 3, col_gap = 3) {
  
  mod_data <- data %>%
    filter(`PTM types` == mod_type, adj.P.Val < sig_thr) %>%
    select(Identifier, Comparison, GN, ModSites) %>%
    distinct() %>%
    mutate(present = 1,
           Gene_Modsite = paste(GN, ModSites, sep = "_")) %>%
    pivot_wider(id_cols = c(Identifier, Gene_Modsite),
                names_from = Comparison, values_from = present,
                values_fill = 0)
  
  if (nrow(mod_data) < 2) {
    cat("  Skipping", mod_type, "— insufficient data\n")
    return(invisible(NULL))
  }
  
  mat <- as.matrix(mod_data[, 3:ncol(mod_data)])
  rownames(mat) <- mod_data$Identifier
  m <- make_comb_mat(mat, mode = "distinct")
  set_sz   <- set_size(m)
  comb_sz  <- comb_size(m)
  plot_ord <- order(comb_degree(m), -comb_sz)
  ordered_counts <- comb_sz[plot_ord]
  comp_labels    <- gsub("_", " ", colnames(mat))
  
  mod_color <- ptm_colors[mod_type]
  if (is.na(mod_color)) mod_color <- "steelblue"
  mod_label <- ptm_labels[mod_type]
  if (is.na(mod_label)) mod_label <- mod_type
  
  top_anno <- HeatmapAnnotation(
    "Protein-PTM" = anno_barplot(
      ordered_counts, ylim = c(0, max(ordered_counts) * 1.1),
      border = FALSE, gp = gpar(fill = mod_color),
      height = unit(5, "cm"),
      axis_param = list(side = "left", gp = gpar(fontsize = 10))
    ),
    annotation_name_side = "left", annotation_name_rot = 90
  )
  
  left_anno <- rowAnnotation(
    set_name = anno_text(
      comp_labels, location = 0, just = "left",
      width = max_text_width(comp_labels) + unit(5, "mm"),
      gp = gpar(fontsize = 10))
  )
  
  right_anno <- rowAnnotation(
    "Protein-PTM" = anno_barplot(
      set_sz, border = FALSE, gp = gpar(fill = mod_color),
      width = unit(4, "cm"))
  )
  
  fname <- paste0("upset_", gsub("[^A-Za-z0-9]","_", tolower(mod_type)),
                  ".pdf")
  pdf(fname, width = 14, height = 8)
  ht <- UpSet(
    m,
    set_order   = order(-set_sz),
    comb_order  = plot_ord,
    top_annotation  = top_anno,
    left_annotation = left_anno,
    right_annotation = right_anno,
    show_row_names  = FALSE,
    column_title    = paste0(mod_label, " Modifications\n(",
                             nrow(mod_data),
                             " sites with adj.p < 0.05)"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    pt_size    = unit(pt_size, "mm"),
    comb_col   = "black",
    bg_col     = "grey95", bg_pt_col = "grey80",
    height     = unit(5, "cm"), width = unit(15, "cm"),
    row_gap    = unit(row_gap, "mm"), column_gap = unit(col_gap, "mm")
  )
  draw(ht, padding = unit(c(2, 10, 2, 2), "mm"))
  dev.off()
  
  # Per-pattern summary CSV
  ordered_patterns <- names(comb_sz)[plot_ord]
  
  get_genes_for_pattern <- function(pattern, mat, mod_data) {
    pb  <- as.logical(as.numeric(strsplit(pattern,"")[[1]]))
    rows <- which(rowSums(mat[, pb,  drop = FALSE]) == sum(pb) &
                    rowSums(mat[, !pb, drop = FALSE]) == 0)
    mod_data %>%
      filter(Identifier %in% rownames(mat)[rows]) %>%
      pull(Gene_Modsite) %>% paste(collapse = "; ")
  }
  
  summary_df <- data.frame(
    Column      = seq_along(ordered_patterns),
    Pattern     = ordered_patterns,
    Count       = ordered_counts,
    Comparisons = sapply(ordered_patterns, function(p) {
      pb <- as.logical(as.numeric(strsplit(p,"")[[1]]))
      paste(comp_labels[pb], collapse = " & ")
    }),
    Gene_Modsites = sapply(ordered_patterns, function(p)
      get_genes_for_pattern(p, mat, mod_data))
  )
  
  write.csv(summary_df,
            paste0("upset_", gsub("[^A-Za-z0-9]","_", tolower(mod_type)),
                   "_summary.csv"),
            row.names = FALSE)
  
  cat("  Saved:", fname, "\n")
  invisible(NULL)
}

for (mod in unique(data_long$`PTM types`))
  create_mod_upset(data_long, mod)

cat("Step 6 complete.\n\n")

# ============================================================
# 7. Publication-Quality Scatter Plots
# ============================================================

cat("Step 7: Scatter plots...\n")

# Standardize column names from CSV (dots → spaces issue)
data_scatter <- data_long %>%
  rename_with(~ gsub("\\.", " ", .), everything()) %>%
  rename_with(~ gsub("PTM types", "PTM_types", .), everything()) %>%
  rename_with(~ gsub("Protein ID", "Protein_ID", .), everything()) %>%
  group_by(Identifier, Comparison) %>%
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Context   = case_when(
      grepl("Ratio", Comparison) ~ "Ratio",
      grepl("SOL",   Comparison) ~ "SOL",
      grepl("INSOL", Comparison) ~ "INSOL"),
    Condition = case_when(
      grepl("Group3_vs_Group1", Comparison) ~ "Age",
      grepl("Control",          Comparison) ~ "Genotype")
  ) %>%
  filter(!is.na(Context), !is.na(Condition))

# Core scatter builder
build_scatter <- function(fracA_data, fracB_data, sig_data,
                          title, filename_base,
                          x_label, y_label,
                          label_top_n = 15) {
  
  plot_data <- fracA_data %>%
    inner_join(fracB_data, by = "Identifier", suffix = c("_A","_B")) %>%
    inner_join(sig_data %>% select(Identifier, logFC_sig = logFC,
                                   pval_sig = adj.P.Val),
               by = "Identifier") %>%
    mutate(
      Size        = -log10(pval_sig),
      is_sig      = pval_sig < 0.05,
      GN          = GN_A,
      PTM_types   = PTM_types_A,
      ModSites    = ModSites_A
    )
  
  top_ids <- plot_data %>%
    filter(is_sig) %>%
    slice_min(pval_sig, n = label_top_n, with_ties = FALSE) %>%
    mutate(Label = paste0(GN, " ", ModSites)) %>%
    select(Identifier, Label)
  
  plot_data <- plot_data %>%
    left_join(top_ids, by = "Identifier")
  
  sig_d    <- filter(plot_data, is_sig)
  nonsig_d <- filter(plot_data, !is_sig)
  
  cat("  Matched sites:", nrow(plot_data),
      "| Significant:", nrow(sig_d), "\n")
  
  r2_all <- round(summary(lm(logFC_B ~ logFC_A, plot_data))$r.squared, 3)
  r2_sig <- if (nrow(sig_d) >= 3)
    round(summary(lm(logFC_B ~ logFC_A, sig_d))$r.squared, 3)
  else NA
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_point(data = nonsig_d,
               aes(x = logFC_A, y = logFC_B, size = Size),
               color = "gray70", alpha = 0.1) +
    geom_point(data = sig_d,
               aes(x = logFC_A, y = logFC_B,
                   color = PTM_types, size = Size),
               alpha = 0.6) +
    scale_color_manual(values = ptm_colors, labels = ptm_labels,
                       name = "Modification Type") +
    scale_size_continuous(
      range = c(0.5, 8), name = "Ratio p-value\n(-log10)",
      breaks = c(1.3, 2, 3, 4, 5),
      limits = c(0.5, max(plot_data$Size, na.rm = TRUE)),
      labels = c("1.3 (p<0.05)","2","3","4","5+")) +
    geom_text_repel(
      data = filter(plot_data, !is.na(Label)),
      aes(x = logFC_A, y = logFC_B, label = Label),
      size = 3.5, fontface = "bold", max.overlaps = 30,
      segment.color = "grey50", box.padding = 0.5) +
    annotate("text",
             x = min(plot_data$logFC_A, na.rm = TRUE) + 0.2,
             y = max(plot_data$logFC_B, na.rm = TRUE) - 0.2,
             label = paste0("All: R² = ", r2_all),
             hjust = 0, size = 4.5, fontface = "bold") +
    annotate("text",
             x = min(plot_data$logFC_A, na.rm = TRUE) + 0.2,
             y = max(plot_data$logFC_B, na.rm = TRUE) - 0.6,
             label = paste0("Significant: R² = ",
                            ifelse(is.na(r2_sig), "NA", r2_sig)),
             hjust = 0, color = "darkred", size = 4.5,
             fontface = "bold") +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal(base_size = 11) +
    theme(axis.text      = element_text(size = 12),
          legend.title   = element_text(face = "bold"),
          plot.title     = element_text(face = "bold", size = 14))
  
  ggsave(paste0(filename_base, ".pdf"), p,
         width = 14, height = 10, dpi = 600)
  ggsave(paste0(filename_base, ".png"), p,
         width = 14, height = 10, dpi = 300)
  
  export_data <- plot_data %>%
    select(GN, PTM_types, ModSites,
           Log2FC_FractionA  = logFC_A,
           Log2FC_FractionB  = logFC_B,
           Log2FC_Ratio      = logFC_sig,
           pvalue_Ratio      = pval_sig,
           Size, is_sig, Label)
  write.csv(export_data, paste0(filename_base, "_data.csv"),
            row.names = FALSE)
  
  cat("  Saved:", filename_base, "\n")
  list(plot = p, r2_all = r2_all, r2_sig = r2_sig)
}

get_ctx <- function(condition, context)
  data_scatter %>%
  filter(Condition == condition, Context == context) %>%
  select(Identifier, GN, PTM_types, ModSites, logFC, adj.P.Val)

# Plot 1: FractionA vs FractionB — Age comparison
build_scatter(
  fracA_data = get_ctx("Age","SOL"),
  fracB_data = get_ctx("Age","INSOL"),
  sig_data   = get_ctx("Age","Ratio"),
  title      = "Age Comparison (Group 3 vs Group 1): Fraction A vs Fraction B",
  filename_base = "scatter_FracA_FracB_age",
  x_label    = "Log2 FC (Fraction A)",
  y_label    = "Log2 FC (Fraction B)")

# Plot 2: FractionA vs FractionB — Genotype comparison
build_scatter(
  fracA_data = get_ctx("Genotype","SOL"),
  fracB_data = get_ctx("Genotype","INSOL"),
  sig_data   = get_ctx("Genotype","Ratio"),
  title      = "Genotype Comparison (Variant vs Control): Fraction A vs Fraction B",
  filename_base = "scatter_FracA_FracB_genotype",
  x_label    = "Log2 FC (Fraction A)",
  y_label    = "Log2 FC (Fraction B)")

# Plot 3: FractionA vs Ratio — Age comparison
build_scatter(
  fracA_data = get_ctx("Age","SOL"),
  fracB_data = get_ctx("Age","Ratio"),
  sig_data   = get_ctx("Age","Ratio"),
  title      = "Age Comparison (Group 3 vs Group 1): Fraction A vs Fraction Ratio",
  filename_base = "scatter_FracA_Ratio_age",
  x_label    = "Log2 FC (Fraction A)",
  y_label    = "Log2 FC (Fraction B / Fraction A Ratio)")

cat("Step 7 complete.\n\n")

# ============================================================
# Summary
# ============================================================

summary_stats <- data_long %>%
  group_by(Comparison) %>%
  summarise(Total       = n(),
            Significant = sum(adj.P.Val < 0.05, na.rm = TRUE),
            Sig_pct     = round(Significant / Total * 100, 1),
            .groups = "drop")
print(summary_stats)

cat("\nOutputs:\n")
cat("  clean_data/*.csv           — 6 cleaned datasets\n")
cat("  organized_PTM_data.csv     — combined long-format data\n")
cat("  volcano_*.pdf              — per-comparison volcano plots\n")
cat("  global_PTM_patterns_*.pdf  — barplot and heatmap\n")
cat("  overview_upset_plot.pdf    — all modifications combined\n")
cat("  upset_*.pdf                — per-modification UpSet plots\n")
cat("  upset_*_summary.csv        — site-level UpSet summaries\n")
cat("  scatter_*.pdf/.png         — publication scatter plots\n")
cat("  scatter_*_data.csv         — scatter plot data exports\n")
cat("\nNext: Run 02_GO_enrichment.R\n")