## ============================
## Pathway Enrichment Script
## Inputs: Excel with columns Ensembl_ID, gene_name
## Databases: GO (BP/MF/CC), KEGG, Reactome, PANTHER
## Outputs: CSVs + dotplots saved alongside the Excel file
## ============================

## 0) Packages ---------------------------------------------------------------
pkg_needed <- c(
  "readxl", "dplyr", "stringr", "tibble", "openxlsx",
  "clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot",
  "ggplot2", "patchwork"   # combined figure
)
# PANTHER is via WebGestaltR (optional)
pkg_optional <- c("WebGestaltR")

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
}
install_if_missing(pkg_needed)
# Bioc packages (clusterProfiler, org.Hs.eg.db, ReactomePA) if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (biopkg in c("clusterProfiler", "org.Hs.eg.db", "ReactomePA")) {
  if (!requireNamespace(biopkg, quietly = TRUE)) BiocManager::install(biopkg, update = FALSE, ask = FALSE)
}
# Optional WebGestaltR for PANTHER
if (!requireNamespace("WebGestaltR", quietly = TRUE)) {
  tryCatch({
    install.packages("WebGestaltR", repos = "https://cloud.r-project.org")
  }, error = function(e) message("Skipping WebGestaltR install (optional): ", e$message))
}

library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(patchwork)

## 1) Paths & gene input -----------------------------------------------------
excel_path <- "/Users/sindypinero/Downloads/TACO/Genes.xlsx"
out_dir <- dirname(excel_path)
prefix  <- file.path(out_dir, "TACO_enrichment")

genes_df <- readxl::read_xlsx(excel_path)

# Sanity check expected columns
stopifnot(all(c("Ensembl_ID", "gene_name") %in% colnames(genes_df)))

# Clean Ensembl IDs (strip versions if any)
genes_df <- genes_df %>%
  mutate(Ensembl_ID = str_replace(Ensembl_ID, "\\.\\d+$", "")) %>%
  distinct(Ensembl_ID, .keep_all = TRUE)

# 2) Map Ensembl -> Entrez ---------------------------------------------------
id_map <- clusterProfiler::bitr(
  genes_df$Ensembl_ID,
  fromType = "ENSEMBL",
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Hs.eg.db
) %>%
  as_tibble()

if (nrow(id_map) == 0) stop("No Ensembl IDs could be mapped to Entrez IDs. Check your input.")
message(sprintf("Mapped %d/%d Ensembl IDs to Entrez.",
                nrow(id_map), nrow(genes_df)))

gene_entrez <- unique(id_map$ENTREZID)

# Optional: background gene universe
# bg_universe <- keys(org.Hs.eg.db, keytype = "ENTREZID")
bg_universe <- NULL

## 3) GO Enrichment (BP/MF/CC) ----------------------------------------------
ego_bp <- enrichGO(
  gene          = gene_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  universe      = bg_universe,
  readable      = TRUE
)

ego_mf <- enrichGO(
  gene          = gene_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
  ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20,
  universe = bg_universe, readable = TRUE
)

ego_cc <- enrichGO(
  gene          = gene_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
  ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20,
  universe = bg_universe, readable = TRUE
)

## 4) KEGG (human = 'hsa') ---------------------------------------------------
ekegg <- enrichKEGG(
  gene          = gene_entrez,
  organism      = "hsa",
  keyType       = "kegg",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.20,
  universe      = if (is.null(bg_universe)) NULL else bg_universe
)

## 5) Reactome ---------------------------------------------------------------
ereact <- ReactomePA::enrichPathway(
  gene          = gene_entrez,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.20,
  universe      = bg_universe,
  readable      = TRUE
)

## 6) (Optional) PANTHER via WebGestaltR ------------------------------------
run_panther <- requireNamespace("WebGestaltR", quietly = TRUE)
epanther <- NULL
if (run_panther) {
  panther_dir <- paste0(prefix, "_PANTHER_WebGestalt")
  dir.create(panther_dir, showWarnings = FALSE, recursive = TRUE)
  message("Running PANTHER (WebGestaltR)...")
  epanther <- tryCatch({
    WebGestaltR::WebGestaltR(
      enrichMethod    = "ORA",
      organism        = "hsapiens",
      interestGene    = gene_entrez,
      interestGeneType= "entrezgene",
      enrichDatabase  = "pathway_PANTHER",
      referenceSet    = if (is.null(bg_universe)) "genome_protein-coding" else bg_universe,
      referenceGeneType = if (is.null(bg_universe)) "ensembl_gene_id" else "entrezgene",
      is.output       = TRUE,
      outputDirectory = panther_dir,
      fdrMethod       = "BH",
      sigMethod       = "fdr",
      fdrThr          = 0.20,
      topThr          = 200
    )
  }, error = function(e) {
    message("WebGestaltR/PANTHER failed (continuing without it): ", e$message)
    NULL
  })
} else {
  message("WebGestaltR not available; skipping PANTHER.")
}

## 7) Save results -----------------------------------------------------------
save_if <- function(ego, name) {
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    out_csv <- paste0(prefix, "_", name, ".csv")
    readr::write_csv(as.data.frame(ego), out_csv)
    message("Saved: ", out_csv)
  }
}
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = "https://cloud.r-project.org")
library(readr)

save_if(ego_bp,  "GO_BP")
save_if(ego_mf,  "GO_MF")
save_if(ego_cc,  "GO_CC")
save_if(ekegg,   "KEGG")
save_if(ereact,  "Reactome")

if (!is.null(epanther) && is.data.frame(epanther) && nrow(epanther) > 0) {
  readr::write_csv(epanther, paste0(prefix, "_PANTHER.csv"))
}

## 8) Quick visualizations (top 20 terms) -----------------------------------
# Label cleaner: remove parentheses content, squeeze spaces, keep first N words,
# and ensure the first character is uppercase. No ellipses.
shorten_labels <- function(x, max_words = 6) {
  vapply(x, function(s) {
    s <- gsub("\\s*\\(.*?\\)", "", s)      # remove parentheticals
    s <- gsub("\\s+", " ", s)              # squeeze spaces
    s <- trimws(s)
    words <- strsplit(s, "\\s+")[[1]]
    if (length(words) > max_words) words <- words[1:max_words]
    s <- paste(words, collapse = " ")
    if (nchar(s) > 0) s <- paste0(toupper(substr(s,1,1)), substr(s,2,nchar(s)))
    s
  }, character(1))
}

dotplot_clean <- function(ego, name, n_show = 20, max_words = 6) {
  enrichplot::dotplot(ego, showCategory = n_show) +
    ggtitle(name) +
    scale_y_discrete(labels = function(lbls) shorten_labels(lbls, max_words = max_words))
}

plot_if <- function(ego, name, n_show = 20, max_words = 6) {
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    p <- dotplot_clean(ego, name, n_show = n_show, max_words = max_words)
    ggsave(filename = paste0(prefix, "_", name, "_dotplot.png"),
           plot = p, width = 8, height = 6, dpi = 300)
    return(p)
  } else {
    return(NULL)
  }
}

# Build each plot (still saving separate PNGs)
p_bp  <- plot_if(ego_bp, "GO_BP", max_words = 6)
p_mf  <- plot_if(ego_mf, "GO_MF", max_words = 6)
p_cc  <- plot_if(ego_cc, "GO_CC", max_words = 6)
p_keg <- plot_if(ekegg,  "KEGG",  max_words = 6)
p_rea <- plot_if(ereact, "Reactome", max_words = 6)

# Combined figure
plist <- Filter(Negate(is.null), list(p_bp, p_mf, p_cc, p_keg, p_rea))
if (length(plist) > 0) {
  layout <- switch(
    length(plist),
    `1` = plist[[1]],
    `2` = plist[[1]] | plist[[2]],
    `3` = (plist[[1]] | plist[[2]]) / plist[[3]],
    `4` = (plist[[1]] | plist[[2]]) / (plist[[3]] | plist[[4]]),
    `5` = ((plist[[1]] | plist[[2]]) / (plist[[3]] | plist[[4]])) / plist[[5]],
    wrap_plots(plist, ncol = 2)
  )
  ggsave(
    filename = paste0(prefix, "_ALL_dotplots.png"),
    plot = layout,
    width = 16, height = 18, dpi = 300
  )
}

## 9) Session info snapshot --------------------------------------------------
sink(paste0(prefix, "_sessionInfo.txt"))
sessionInfo()
sink()

message("Done. Results written to: ", out_dir)