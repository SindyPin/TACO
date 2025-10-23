# KEGG version: Release 113.0, January 1, 2025

##################################################
# Updated Exp1_KEGG_chunking_GADI.R
#
# Modifications:
# - Reduced the number of cores (nCores) from 24 to 12.
# - Reduced the chunk size from 50,000 to 10,000.
# - Added extra garbage collection (gc()) calls and logging.
##################################################

##################################################
# 0. SETUP LOGGING AND MEMORY MONITORING
##################################################
# Create a log file to capture messages and errors
log_file <- file(paste0("error_log_", Sys.getpid(), ".txt"), open = "w")
sink(log_file, type = "message")

# Function to return current memory usage in MB
mem_usage <- function() {
  gc()  # trigger garbage collection
  mem <- gc()
  as.numeric(mem[2,2]) / 1024  # second row, second column, in MB
}

# Function to log memory usage with a timestamp and message
log_memory <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s - Memory usage: %.2f MB\n", 
              timestamp, message, mem_usage()))
}

##################################################
# 1. INSTALL AND LOAD LIBRARIES
##################################################
log_memory("Starting library loading")

# List of required packages
required_packages <- c(
  "Rgraphviz", "dplyr", "igraph", "pheatmap", "ggplot2", 
  "tidyr", "caret", "clusterProfiler", "org.Hs.eg.db", 
  "enrichplot", "tibble", "foreach", "doParallel", "stats", "msigdbr"
)

# Load each package using library() (assumes packages are installed)
for(pkg in required_packages) {
  library(pkg, character.only = TRUE)
}

log_memory("Finished loading libraries")

##################################################
# 2. READ & REFORMAT DATA
##################################################
log_memory("Starting data reading")

# Define a function to safely read and process the input data
safe_read_process <- function(file_path, prefix) {
  tryCatch({
    # Read the tab-delimited file
    df <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    log_memory(paste("Read", prefix, "data"))
    
    # Process the data:
    # - Remove the 'Gene_Type' column
    # - Convert from wide to long format (one row per subject)
    # - Create a 'patientID' using the prefix (e.g., "COVID" or "PASC") and a substring of the original column name
    # - Compute the mean expression per gene per patient and then pivot back to wide format
    df_clean <- df %>%
      dplyr::select(-Gene_Type) %>%
      pivot_longer(
        cols = starts_with("Subj_"), 
        names_to = "original_col_name",
        values_to = "expression"
      ) %>%
      mutate(
        patientID = paste0(prefix, "_", sub("(Subj_[^T]+)T.*", "\\1", original_col_name))
      ) %>%
      group_by(Ensembl_Gene_ID, Gene_Symbol, patientID) %>%
      summarize(
        mean_expression = mean(expression, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_wider(
        id_cols = c(Ensembl_Gene_ID, Gene_Symbol),
        names_from = patientID,
        values_from = mean_expression
      )
    log_memory(paste("Processed", prefix, "data"))
    return(df_clean)
  }, error = function(e) {
    cat(sprintf("Error processing %s data: %s\n", prefix, e$message))
    return(NULL)
  })
}

# Process COVID and PASC data using the safe_read_process function
df_covid_clean <- safe_read_process(
  "/scratch/sq95/sp6154/DCE/Positive_Acute_COVID_control_subset_subjects_samples_pcgenes.txt",
  "COVID"
)

df_pasc_clean <- safe_read_process(
  "/scratch/sq95/sp6154/DCE/PASC_case_subset_subjects_samples_pcgenes.txt",
  "PASC"
)

if (is.null(df_covid_clean) || is.null(df_pasc_clean)) {
  stop("Failed to process input data")
}

##################################################
# 3. MERGE AND FILTER DATA
##################################################
log_memory("Starting data merge")

# Merge the cleaned COVID and PASC data on Ensembl_Gene_ID and Gene_Symbol
df_merged <- full_join(df_covid_clean, df_pasc_clean, 
                       by = c("Ensembl_Gene_ID", "Gene_Symbol"))

# Remove intermediate data to free up memory
rm(df_covid_clean, df_pasc_clean)
gc()

# ---- Local KEGG enrichment using local GMT file ----
log_memory("Performing local KEGG enrichment")

# Set the GMT file path (ensure this file exists at this location)
gmt_file <- "/scratch/sq95/sp6154/DCE/KEGG_gene_sets.gmt"

# Define a function to perform local KEGG enrichment using the GMT file
safe_local_kegg_enrichment <- function(ensembl_genes, prefix, gmt_file) {
  tryCatch({
    # Convert Ensembl IDs to Entrez IDs using clusterProfiler's bitr function
    gene_df <- bitr(ensembl_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    entrez_genes <- gene_df$ENTREZID
    
    # Read the local GMT file (returns a data frame with columns "term" and "gene")
    local_kegg <- read.gmt(gmt_file)
    
    # Perform enrichment analysis using enricher()
    enriched <- enricher(gene = entrez_genes, TERM2GENE = local_kegg, pvalueCutoff = 0.05)
    
    if (!is.null(enriched) && nrow(as.data.frame(enriched)) > 0) {
      # Return the unique set of genes (Entrez IDs) present in enriched KEGG pathways
      return(unique(unlist(strsplit(enriched@result$geneID, "/"))))
    }
    return(character())
    
  }, error = function(e) {
    cat(sprintf("Error in local KEGG enrichment for %s: %s\n", prefix, e$message))
    return(character())
  })
}

# Use the list of Ensembl IDs from the merged data for both groups
genes_in_kegg_covid <- safe_local_kegg_enrichment(unique(df_merged$Ensembl_Gene_ID), "COVID", gmt_file)
genes_in_kegg_pasc  <- safe_local_kegg_enrichment(unique(df_merged$Ensembl_Gene_ID), "PASC", gmt_file)

# Combine the enriched genes from both COVID and PASC groups
union_kegg_genes <- union(genes_in_kegg_covid, genes_in_kegg_pasc)

log_memory("Filtering merged data")

# Filter df_merged to keep only genes that are enriched in KEGG
conversion_df <- bitr(df_merged$Ensembl_Gene_ID,
                      fromType = "ENSEMBL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

df_merged_filtered <- df_merged[df_merged$Ensembl_Gene_ID %in% 
                                  conversion_df$ENSEMBL[conversion_df$ENTREZID %in% union_kegg_genes], ]

rm(df_merged, conversion_df)
gc()

##################################################
# 4. PREPARE DATA FOR ANALYSIS
##################################################
log_memory("Preparing data for analysis")

# Convert the filtered merged data to a data frame and set rownames to Ensembl_Gene_ID
df_merged_filtered <- as.data.frame(df_merged_filtered)
rownames(df_merged_filtered) <- df_merged_filtered$Ensembl_Gene_ID

# Remove the first two columns (Ensembl_Gene_ID and Gene_Symbol)
df_merged2 <- df_merged_filtered[, -c(1, 2)]
# Transpose the data so that rows become samples and columns become genes
expression_data <- t(df_merged2)

# Clean up intermediate objects and free memory
rm(df_merged_filtered, df_merged2)
gc()

##################################################
# 5. SETUP PARALLEL ENVIRONMENT
##################################################
log_memory("Setting up parallel environment")

# Set up a parallel cluster; here nCores is set to 48 (adjust as needed)
nCores <- 48
library(doParallel)
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Increase recursion limit on each worker node
clusterCall(cl, function() options(expressions = 10000))

##################################################
# 6. PREPARE ANALYSIS PARAMETERS
##################################################
# Extract sample names and define condition labels for analysis
sample_names <- rownames(expression_data)
condition <- ifelse(grepl("^COVID", sample_names), "COVID", "PASC")
condition <- factor(condition, levels = c("COVID", "PASC"))

# Extract gene names (columns) from expression_data
genes <- colnames(expression_data)
nGenes <- length(genes)
total_pairs <- nGenes * (nGenes - 1)

# Set chunk size to 10,000 gene pairs per chunk
chunk_size <- 10000
num_chunks <- ceiling(total_pairs / chunk_size)

log_memory(sprintf("Prepared for analysis with %d genes and %d chunks", nGenes, num_chunks))

# Export variables required for parallel processing to worker nodes
clusterExport(cl, varlist = c("nGenes", "genes", "expression_data", "condition"))

##################################################
# 7. PROCESS CHUNKS
##################################################
# Define a function to process a chunk of gene pairs
process_chunk <- function(start_idx, end_idx) {
  # Define a log file for this chunk
  chunk_log_file <- paste0("chunk_log_", start_idx, "_", end_idx, ".txt")
  
  results_chunk <- foreach(idx = start_idx:end_idx, 
                           .combine = 'list', 
                           .packages = c("stats")) %dopar% {
                             # Run garbage collection periodically
                             if (idx %% 1000 == 0) gc()
                             
                             # Convert the linear index to a pair (i, j)
                             i <- ((idx - 1) %/% (nGenes - 1)) + 1
                             j <- ((idx - 1) %% (nGenes - 1)) + 1
                             if (j >= i) j <- j + 1
                             if (j > nGenes) return(NULL)
                             
                             predictor_gene <- genes[i]
                             outcome_gene <- genes[j]
                             
                             # Prepare a data frame for regression
                             df_reg <- data.frame(
                               predictor = expression_data[, predictor_gene],
                               outcome = expression_data[, outcome_gene],
                               condition = condition
                             )
                             
                             # Run linear regression with interaction and extract the coefficient and p-value
                             result <- tryCatch({
                               model <- lm(outcome ~ predictor * condition, data = df_reg)
                               coef_summary <- summary(model)$coefficients
                               
                               if ("predictor:conditionPASC" %in% rownames(coef_summary)) {
                                 beta_int <- coef_summary["predictor:conditionPASC", "Estimate"]
                                 p_val <- coef_summary["predictor:conditionPASC", "Pr(>|t|)"]
                               } else {
                                 beta_int <- NA
                                 p_val <- NA
                               }
                               
                               list(
                                 predictor = predictor_gene,
                                 outcome = outcome_gene,
                                 beta_int = beta_int,
                                 p_value = p_val
                               )
                             }, error = function(e) {
                               cat(sprintf("Error in model fitting: %s\n", e$message), 
                                   file = chunk_log_file, append = TRUE)
                               list(
                                 predictor = predictor_gene,
                                 outcome = outcome_gene,
                                 beta_int = NA,
                                 p_value = NA
                               )
                             })
                             
                             result
                           }
  
  # Remove any NULL elements from the results
  results_chunk <- results_chunk[!sapply(results_chunk, is.null)]
  
  # Ensure each result list contains the required fields
  required_names <- c("predictor", "outcome", "beta_int", "p_value")
  results_chunk <- lapply(results_chunk, function(x) {
    for (nm in required_names) {
      if (is.null(x[[nm]])) x[[nm]] <- NA
    }
    x[required_names]
  })
  
  # Convert the list of results into a data frame
  if (length(results_chunk) > 0) {
    results_chunk <- do.call(rbind, lapply(results_chunk, function(x) {
      as.data.frame(x, stringsAsFactors = FALSE)
    }))
  } else {
    results_chunk <- data.frame(
      predictor = character(),
      outcome = character(),
      beta_int = numeric(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  # Save the results of this chunk to an RDS file
  saveRDS(results_chunk, file = paste0("DCE_results_", start_idx, "_", end_idx, ".rds"))
  rm(results_chunk)
  gc()
}

log_memory("Starting chunk processing")
for (chunk in 1:num_chunks) {
  start_idx <- (chunk - 1) * chunk_size + 1
  end_idx <- min(chunk * chunk_size, total_pairs)
  
  log_memory(sprintf("Processing chunk %d of %d", chunk, num_chunks))
  process_chunk(start_idx, end_idx)
  log_memory(sprintf("Finished processing chunk %d of %d", chunk, num_chunks))
  gc()
}

##################################################
# 8. COMBINE RESULTS
##################################################
log_memory("Combining results")

# Select only the batch files using a regex pattern (e.g., "DCE_results_1_10000.rds")
all_files <- list.files(pattern = "^DCE_results_\\d+_\\d+\\.rds$")
all_results <- do.call(rbind, lapply(all_files, readRDS))

# Convert p_value to numeric and filter out NA values
all_results$p_value <- as.numeric(all_results$p_value)
all_results <- subset(all_results, !is.na(p_value))
all_results$adj_p_value <- p.adjust(all_results$p_value, method = "fdr")

all_results_clean <- na.omit(all_results)

# Save the final results to both RDS and CSV files
saveRDS(all_results, "/scratch/sq95/sp6154/DCE/DCE_results_final_300.rds")
write.csv(all_results, 
          file = "/scratch/sq95/sp6154/DCE/DCE_results_final_clean_300.csv", 
          row.names = FALSE)

log_memory("Saved final results")

##################################################
# 9. OPTIONAL VISUALIZATION
##################################################
log_memory("Creating visualization")

tryCatch({
  # Filter for significant differential causal effects:
  # p-value < 0.05 and absolute beta coefficient > 0.05
  significant_edges <- subset(all_results, 
                              !is.na(adj_p_value) & 
                                adj_p_value < 0.05 & 
                                abs(beta_int) > 0.05)
  
  if (nrow(significant_edges) > 0) {
    # Create an edge list using predictor and outcome columns
    edge_list <- significant_edges[, c("predictor", "outcome")]
    g <- graph_from_data_frame(edge_list, directed = TRUE)
    
    # Save the network plot as a PDF
    pdf("/scratch/sq95/sp6154/DCE/network_visualization.pdf")
    plot(g, vertex.size = 10, 
         vertex.label.cex = 0.7, 
         edge.arrow.size = 0.7,
         main = "Network of Significant Differential Causal Effects")
    dev.off()
  }
}, error = function(e) {
  cat("Error in visualization:", e$message, "\n")
})

##################################################
# 10. CLEANUP
##################################################
log_memory("Starting cleanup")

# Stop the parallel cluster
stopCluster(cl)
sink(type = "message")
close(log_file)

log_memory("Script completed")
