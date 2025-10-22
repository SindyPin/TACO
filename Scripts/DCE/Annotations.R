# Annotate DCE results with COVID/LC relevance & rationale (R version)
# - Edit the paths below if your filenames differ.
#install.packages("data.table")
suppressPackageStartupMessages({
  library(data.table)
})

# ---- Paths ----
input_path      <- "/Users/sindypinero/Downloads/TACO/DCE_COVID_LC_Results.csv"         # your DCE table (CSV or TSV)
relevance_path  <- "/Users/sindypinero/Downloads/TACO/covid_lc_kegg_relevance.csv"       # file generated earlier
output_path     <- "/Users/sindypinero/Downloads/TACO/DCE_COVID_LC_Results_annotated.csv"

# ---- Load (fread auto-detects delimiter) ----
dce <- fread(input_path)
rel <- fread(relevance_path)

# ---- Basic checks ----
if (!"pathway_source" %in% names(dce)) {
  stop("Input file must contain a 'pathway_source' column.")
}
if (!all(c("pathway_source","Relevance","Rationale") %in% names(rel))) {
  stop("Relevance file must have columns: pathway_source, Relevance, Rationale.")
}

# ---- Normalize KEGG IDs (trim/ensure character) ----
dce[, pathway_source := trimws(as.character(pathway_source))]
rel[, pathway_source := trimws(as.character(pathway_source))]

# ---- Keep only needed columns from relevance table and deduplicate ----
rel_small <- unique(rel[, .(pathway_source, Relevance, Rationale)])

# ---- Left join by KEGG ID ----
annot <- merge(dce, rel_small, by = "pathway_source", all.x = TRUE, sort = FALSE)

# ---- Report unmatched KEGG IDs ----
unmatched <- unique(annot[is.na(Relevance) & !is.na(pathway_source), pathway_source])
if (length(unmatched) > 0) {
  message(sprintf("[WARN] %d KEGG IDs had no match in relevance file:", length(unmatched)))
  for (k in unmatched) message("  - ", k)
} else {
  message("All KEGG IDs matched a relevance/rationale entry.")
}

# ---- Save output ----
fwrite(annot, file = output_path)
message(sprintf("✅ Saved annotated table to: %s", normalizePath(output_path)))
