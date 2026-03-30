# ============================================================
# Project: PTM Proteomics — Drosophila Variant × Aging Study
# Script:  02_GO_enrichment.R
# Purpose: GO overrepresentation analysis (BP, MF, CC) for
#          all six comparisons in the PTM dataset. Runs
#          clusterProfiler::enrichGO for significant proteins
#          (adj.p < 0.05) per comparison. Exports enrichment
#          tables (CSV) and dotplots sorted by p-value (PDF).
#          Also identifies and exports proteins that could
#          not be mapped to Entrez IDs for QC purposes.
# Input:   organized_PTM_data.csv
#            (output of Script 01)
# Output:  GO_table_*.csv
#          GO_dotplot_*.pdf
#          unmapped_significant_proteins.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
#          Organism: Drosophila melanogaster (org.Dm.eg.db).
#          Background = all proteins detected in the dataset.
#          p-value-only filtering: adj.p < 0.05, no FC
#          threshold applied in the enrichment step itself.
#          Gene name cleaning (Step 1) handles common
#          Drosophila annotation quirks: "Dmel\" prefixes,
#          isoform suffixes (-RA, -RB), and "Isoform" text.
# ============================================================

library(tidyverse)
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)
library(ggplot2)

# NOTE: Set your working directory to the folder containing input files
# setwd("path/to/your/data")

# ============================================================
# 1. Load Data and Clean Gene Names
# ============================================================

cat("Loading data...\n")

data_long <- read.csv("organized_PTM_data.csv")

# Clean common Drosophila annotation quirks in gene names
data_fixed <- data_long %>%
  mutate(
    GN_original = GN,
    GN = case_when(
      # Remove "Dmel\" prefix
      grepl("^Dmel\\\\", GN) ~ gsub("^Dmel\\\\", "", GN),
      # Remove isoform suffix (e.g. -RA, -RB)
      grepl("-R[A-Z]$", GN)  ~ gsub("-R[A-Z]$", "", GN),
      # Keep only first word if "Isoform" appears
      grepl("Isoform", GN)   ~ word(GN, 1),
      TRUE ~ GN
    )
  )

cat("Gene names cleaned.\n\n")

# ============================================================
# 2. Map Gene Symbols to Entrez IDs
# ============================================================

cat("Mapping gene symbols to Entrez IDs...\n")

all_genes    <- unique(data_fixed$GN)
gene_mapping <- bitr(all_genes,
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Dm.eg.db)

cat("  Total genes:   ", length(all_genes), "\n")
cat("  Mapped:        ", nrow(gene_mapping), "\n")
cat("  Unmapped:      ", length(all_genes) - nrow(gene_mapping), "\n\n")

# ============================================================
# 3. Export Unmapped Significant Proteins (QC)
# ============================================================

cat("Identifying unmapped significant proteins...\n")

# Handle column name format (CSV uses dots, Excel uses spaces)
protein_col <- if ("Protein.ID" %in% names(data_fixed)) "Protein.ID" else "Protein ID"

unmapped_sig <- data_fixed %>%
  filter(adj.P.Val < 0.05,
         !(GN %in% gene_mapping$SYMBOL)) %>%
  select(GN, GN_original,
         Protein_ID = all_of(protein_col)) %>%
  distinct()

write.csv(unmapped_sig, "unmapped_significant_proteins.csv",
          row.names = FALSE)
cat("  Unmapped significant proteins:", nrow(unmapped_sig), "\n")
cat("  Saved: unmapped_significant_proteins.csv\n\n")

# ============================================================
# 4. GO Enrichment per Comparison and Ontology
# ============================================================

cat("Running GO enrichment analysis...\n\n")

background_entrez <- gene_mapping %>%
  filter(SYMBOL %in% unique(data_fixed$GN)) %>%
  pull(ENTREZID)

comparisons   <- unique(data_fixed$Comparison)
ontologies    <- c("BP","MF","CC")
results_count <- 0

for (comp in comparisons) {
  cat("--- Comparison:", comp, "---\n")
  
  sig_genes <- data_fixed %>%
    filter(Comparison == comp, adj.P.Val < 0.05) %>%
    pull(GN) %>% unique()
  
  cat("  Significant genes:", length(sig_genes), "\n")
  
  if (length(sig_genes) < 3) {
    cat("  Skipping: too few genes (<3)\n\n"); next
  }
  
  genes_entrez <- gene_mapping %>%
    filter(SYMBOL %in% sig_genes) %>%
    pull(ENTREZID)
  
  if (length(genes_entrez) < 3) {
    cat("  Skipping: too few mapped genes (<3)\n\n"); next
  }
  
  for (ont in ontologies) {
    
    ego <- tryCatch(
      enrichGO(
        gene          = genes_entrez,
        universe      = background_entrez,
        OrgDb         = org.Dm.eg.db,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      ),
      error = function(e) {
        cat("  Error in", ont, ":", e$message, "\n"); NULL
      })
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      n_terms <- nrow(as.data.frame(ego))
      cat("  ", ont, ":", n_terms, "terms\n")
      
      safe_comp <- gsub("[^A-Za-z0-9]","_", comp)
      
      write.csv(as.data.frame(ego),
                paste0("GO_table_", safe_comp, "_", ont, ".csv"),
                row.names = FALSE)
      
      p <- dotplot(ego, showCategory = 20, font.size = 10,
                   title = paste(gsub("_"," ", comp), "(", ont, ")"))
      ggsave(paste0("GO_dotplot_", safe_comp, "_", ont, ".pdf"),
             p, width = 10, height = 8, dpi = 600)
      
      results_count <- results_count + 1
      
    } else {
      cat("  ", ont, ": no significant terms\n")
    }
  }
  cat("\n")
}

# ============================================================
# 5. Summary
# ============================================================

cat("GO enrichment complete.\n")
cat("  Comparisons tested:    ", length(comparisons), "\n")
cat("  Ontologies per comp:    3 (BP, MF, CC)\n")
cat("  Successful enrichments:", results_count, "\n\n")
cat("Outputs:\n")
cat("  unmapped_significant_proteins.csv\n")
cat("  GO_table_*_[BP|MF|CC].csv\n")
cat("  GO_dotplot_*_[BP|MF|CC].pdf\n")