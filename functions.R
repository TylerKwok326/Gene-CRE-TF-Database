library(tidyverse)
library(DESeq2)
library(BiocManager)
library(GEOquery)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(biomaRt)
library(stringr)
library(tidyr)
library(fgsea)

# Getting the metadata from GEO
get_metdata <- function(gse_number, gse_matrix_status = TRUE, list_number) {
  gse_object <- getGEO(gse_number, GSEMatrix = gse_matrix_status)
  gse_object <- gse_object[[list_number]]
  gse_metadata <- pData(gse_object)
  return(gse_metadata)
}

# Getting the count matrix
get_count_matrix <- function(count_file, col_to_row = TRUE, colname = NULL) {
  raw_counts <- read_tsv(count_file)
  if (col_to_row == TRUE) {
    raw_counts <- raw_counts %>%
      column_to_rownames(var = colname)
    raw_count_matrix <- as.matrix(raw_counts)
    return(raw_count_matrix)
  } else {
    raw_count_matrix <- as.matrix(raw_counts)
    return(raw_count_matrix)
  }
}

# Filtering the metadata
filter_metadata <- function(metadata, col_title, col_group, n = 4, filter_type = "head", exclude_row = NULL) {
  filtered_metadata <- metadata %>%
    dplyr::select(title = all_of(col_title), treatment = all_of(col_group)) %>%
    {
      if (filter_type == "head") head(., n) else tail(., n)
    } %>%
    mutate(condition = if_else(treatment == 'resting', 'control', 'treatment'))
  
  if (!is.null(exclude_row)) {
    filtered_metadata <- filtered_metadata[!row.names(filtered_metadata) %in% exclude_row, ]
  }
  
  return(filtered_metadata)
}

# Filtering the count matrix
filter_non_zero <- function(matrix, n=0, start=1, end=NULL) {
  matrix <- matrix[rowSums(matrix > 0) >= n, ]
  matrix <- matrix[,start:end]
  return(matrix)
}
# Run deseq2 and returns results as a dataframe
deseq_results <- function(metadata, count_matrix, design, reference) {
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = metadata,
                                design = as.formula(paste("~", design)))
  dds[[design]] <- relevel(factor(dds[[design]]), ref=reference)
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[order(res$padj), ]
  return(as.data.frame(res))
}

# Plot Volcano Plot
create_volcano_plot <- function(res_dataframe, title = "Volcano Plot", logFC_threshold = 1, padj_threshold=0.05, n_genes_label=5) {
  # Filter out rows with NA padj and padj <= 0
  res_dataframe <- res_dataframe[!is.na(res_dataframe$padj) & res_dataframe$padj > 0, ]
  
  # Add significance column
  res_dataframe <- res_dataframe %>%
    mutate(significance = if_else(abs(log2FoldChange) >= logFC_threshold & padj < padj_threshold, 
                                  "Significant", "Not Significant"))
  
  # Calculate number of significant genes
  n_sig <- sum(res_dataframe$significance == "Significant")
  
  # Find top significant genes for labeling
  top_genes <- res_dataframe %>%
    filter(significance == "Significant") %>%
    arrange(desc(abs(log2FoldChange))) %>%
    head(n_genes_label)
  
  # Create annotation text
  annotation_text <- paste("Significant genes:", n_sig)
  
  # Create volcano plot
  volcano_plot <- ggplot(res_dataframe, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40", alpha = 0.5, linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", alpha = 0.5, linewidth = 0.4) +
    geom_text_repel(
      data = top_genes,
      aes(label = rownames(top_genes)),
      size = 4,
      color = "black",
      box.padding = 0.8,
      point.padding = 0.5,
      min.segment.length = 0.2,
      max.overlaps = 20,
      nudge_x = 0.3,
      nudge_y = 0.3,
      direction = "both",
      force = 1.5,
      max.time = 2,
      max.iter = 10000
    ) +
    annotate("text",
             x = min(res_dataframe$log2FoldChange, na.rm = TRUE) + 0.5,
             y = max(-log10(res_dataframe$padj), na.rm = TRUE) - 0.5,
             label = annotation_text,
             hjust = 0, vjust = 1, size = 5, color = "red") +
    scale_color_manual(values = c("Not Significant" = "lightblue", "Significant" = "red")) +
    labs(title = title,
         x = expression(Log[2]~"Fold Change"),
         y = expression(-Log[10]~"Adjusted p-value")) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "none",
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    coord_cartesian(
      xlim = c(min(res_dataframe$log2FoldChange, na.rm = TRUE) - 0.5, 
               max(res_dataframe$log2FoldChange, na.rm = TRUE) + 0.5),
      ylim = c(0, max(-log10(res_dataframe$padj), na.rm = TRUE) + 0.5)
    )
  
  return(volcano_plot)
}


# Change NCBI IDs to HGNC symbols
convert_ncbi_to_hgnc <- function(res_dataframe) {
  # Extract NCBI gene IDs
  ncbi_ids <- rownames(res_dataframe)
  
  # Connect to Ensembl BioMart
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  # Get gene symbols corresponding to NCBI gene IDs
  gene_conversion <- getBM(
    attributes = c("entrezgene_id", "hgnc_symbol", "external_gene_name", "description"),
    filters = "entrezgene_id",
    values = ncbi_ids,
    mart = mart
  )
  
  # Ensure proper mapping by removing duplicates
  gene_conversion <- gene_conversion[!duplicated(gene_conversion$entrezgene_id), ]
  rownames(gene_conversion) <- gene_conversion$entrezgene_id
  
  # Map NCBI IDs to HGNC symbols
  new_rownames <- gene_conversion[as.character(ncbi_ids), "hgnc_symbol"]
  
  # Fallback to external gene name if HGNC symbol is missing
  new_rownames[is.na(new_rownames) | new_rownames == ""] <- 
    gene_conversion[as.character(ncbi_ids[is.na(new_rownames) | new_rownames == ""]), "external_gene_name"]
  
  # Fallback to NCBI ID if both HGNC and external gene name are missing
  new_rownames[is.na(new_rownames) | new_rownames == ""] <- 
    ncbi_ids[is.na(new_rownames) | new_rownames == ""]
  
  # Identify duplicate gene symbols
  dup_genes <- new_rownames[duplicated(new_rownames)]
  if (length(dup_genes) > 0) {
    message("Duplicate gene symbols found:")
    print(dup_genes)
  }
  
  # Make duplicate gene symbols unique
  new_rownames_unique <- make.unique(new_rownames)
  
  # Assign new unique row names to the dataframe
  rownames(res_dataframe) <- new_rownames_unique
  
  return(res_dataframe)
}

# Plot PCA
perform_pca_plot <- function(count_matrix, metadata, pc_x = 1, pc_y = 2) {
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ title)
  
  # Perform rlog transformation
  rlog_counts <- rlog(dds, blind = TRUE)
  rlog_matrix <- assay(rlog_counts)
  
  # Remove zero-variance genes
  filtered_rlog <- rlog_matrix[apply(rlog_matrix, 1, var) > 0, ]
  
  # Perform PCA
  pca_res <- prcomp(t(filtered_rlog), scale. = TRUE, center = TRUE)
  percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)
  
  # Check if selected PCs exist
  if (pc_x > ncol(pca_res$x) | pc_y > ncol(pca_res$x)) {
    stop("Selected PCs exceed available principal components.")
  }
  
  # Create PCA dataframe
  pca_df <- data.frame(
    PC_x = pca_res$x[, pc_x], 
    PC_y = pca_res$x[, pc_y], 
    condition = metadata$title
  )
  pca_df$condition <- gsub("RNA-seq_rep[0-9]+", "", metadata$title)
  
  # Create PCA plot
  pca_plot <- ggplot(pca_df, aes(x = PC_x, y = PC_y, color = condition)) +
    geom_point(size = 4, alpha = 0.8) +
    labs(
      title = paste0("PCA of RNA-seq Data (PC", pc_x, " vs PC", pc_y, ")"),
      x = paste0("PC", pc_x, " (", percent_var[pc_x], "% variance)"),
      y = paste0("PC", pc_y, " (", percent_var[pc_y], "% variance)"),
      color = "Condition"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(color = "black", size = 12),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 10),
      plot.margin = margin(10, 10, 30, 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")) +
    guides(color = guide_legend(override.aes = list(size = 5), nrow = 2))
  
  return(pca_plot)
}

# Perform GSEA
perform_gsea_plot <- function(res_dataframe, geneset_file, min_size = 15, max_size = 500, top_n = 5) {
  # Load gene sets
  hallmark_geneset <- gmtPathways(geneset_file)
  
  # Rank genes based on log2FoldChange from DESeq2 results
  ranked_genes <- res_dataframe %>%
    filter(!is.na(log2FoldChange)) %>%
    arrange(desc(log2FoldChange)) %>%
    rownames_to_column(var = 'genes') %>%
    dplyr::select(genes, log2FoldChange)
  
  rank <- setNames(ranked_genes$log2FoldChange, ranked_genes$genes)
  
  # Run fgsea
  fgsea_res <- fgsea(
    pathways = hallmark_geneset,
    stats = rank,
    minSize = min_size, 
    maxSize = max_size
  ) %>% arrange(desc(NES))  # Sort results by NES
  
  # Separate positive and negative pathways
  positive_pathways <- fgsea_res %>% filter(NES > 0)
  negative_pathways <- fgsea_res %>% filter(NES < 0)
  
  # Select top pathways
  top_pos <- positive_pathways %>% arrange(desc(NES)) %>% slice_head(n = top_n)
  top_neg <- negative_pathways %>% arrange(NES) %>% slice_head(n = top_n)
  
  # Combine and label pathways
  pathways_tog <- bind_rows(top_pos, top_neg) %>%
    mutate(condition = if_else(NES > 0, "positive", "negative"))
  
  # Generate bar plot
  gsea_plot <- ggplot(pathways_tog, aes(x = reorder(pathway, NES), y = NES, fill = condition)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    coord_flip() +
    scale_fill_manual(
      values = c('positive' = '#FF8282', 'negative' = '#0072B2'),
      labels = c('Negative NES', 'Positive NES'),
      name = "Direction"
    ) +
    labs(
      title = "GSEA Pathway Enrichment",
      x = NULL,
      y = "Normalized Enrichment Score (NES)",
      caption = paste0("Top ", top_n, " significantly enriched pathways per direction")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.text.y = element_text(face = "plain", size = 10, color = "black"),
      axis.text.x = element_text(face = "plain", size = 10, color = "black"),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 9, color = "gray40", hjust = 0),
      plot.margin = unit(c(1,1,1,1), "cm")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40))
  
  return(gsea_plot)
}

# Regulatory GSEA 
regulatory_gsea_plot <- function(res_dataframe, geneset_file, min_size = 15, max_size = 500, top_n = 5) {
  # Load gene sets
  hallmark_geneset <- gmtPathways(geneset_file)
  
  # Rank genes based on log2FoldChange from DESeq2 results
  ranked_genes <- res_dataframe %>%
    filter(!is.na(avg_CRE_log2FC)) %>%
    arrange(desc(avg_CRE_log2FC)) %>%
    dplyr::select(Gene_ID, avg_CRE_log2FC)
  
  rank <- setNames(ranked_genes$avg_CRE_log2FC, ranked_genes$Gene_ID)
  
  # Run fgsea
  fgsea_res <- fgsea(
    pathways = hallmark_geneset,
    stats = rank,
    minSize = min_size, 
    maxSize = max_size
  ) %>% arrange(desc(NES))  # Sort results by NES
  
  if (nrow(fgsea_res) == 0) {
    stop("No pathways met the criteria for GSEA.")
  }
  
  # Separate positive and negative pathways
  positive_pathways <- fgsea_res %>% filter(NES > 0)
  negative_pathways <- fgsea_res %>% filter(NES < 0)
  
  # Select top pathways
  top_pos <- positive_pathways %>% arrange(desc(NES)) %>% slice_head(n = top_n)
  top_neg <- negative_pathways %>% arrange(NES) %>% slice_head(n = top_n)
  
  if (nrow(top_pos) == 0 && nrow(top_neg) == 0) {
    stop("No pathways met the criteria for top NES values.")
  }
  
  # Combine and label pathways
  pathways_tog <- bind_rows(top_pos, top_neg) %>%
    mutate(condition = if_else(NES > 0, "positive", "negative"))
  
  # Generate bar plot
  gsea_plot <- ggplot(pathways_tog, aes(x = reorder(pathway, NES), y = NES, fill = condition)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    coord_flip() +
    scale_fill_manual(
      values = c('positive' = '#FF8282', 'negative' = '#0072B2'),
      labels = c('Positive NES', 'Negative NES'),
      name = "Direction"
    ) +
    labs(
      title = "Regulatory GSEA Pathway Enrichment",
      x = NULL,
      y = "Normalized Enrichment Score (NES)",
      caption = paste0("Top ", top_n, " significantly enriched pathways per direction")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.text.y = element_text(face = "plain", size = 10, color = "black"),
      axis.text.x = element_text(face = "plain", size = 10, color = "black"),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 9, color = "gray40", hjust = 0),
      plot.margin = unit(c(1,1,1,1), "cm")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40))
  
  return(gsea_plot)
}

convert_ensembl_to_entrez <- function(df, gene_column, version) {
  # Connect to Ensembl BioMart (GRCh37)
  mart <- useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     version = version)
  
  # Retrieve mapping from Ensembl Gene ID to Entrez ID
  mapping <- getBM(attributes = c("ensembl_gene_id", "entrezgene"),
                   filters = "ensembl_gene_id",
                   values = df[[gene_column]],
                   mart = mart)
  
  # Merge with the original dataframe to add Entrez IDs
  df <- left_join(df, mapping, by = setNames("ensembl_gene_id", gene_column))
  
  return(df)
}

convert_entrez_to_ensembl <- function(df, gene_column) {
  # Connect to Ensembl BioMart
  mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  
  # Retrieve mapping from Entrez ID to Ensembl Gene ID
  mapping <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"),
                   filters = "entrezgene_id",
                   values = df[[gene_column]],
                   mart = mart)
  
  # Convert gene_column in df to integer to match the mapping's 'entrezgene_id'
  df[[gene_column]] <- as.integer(df[[gene_column]])
  
  # Merge with the original dataframe to add Ensembl Gene IDs
  df <- left_join(df, mapping, by = setNames("entrezgene_id", gene_column))
  
  return(df)
}

convert_ensembl_to_hgnc <- function(df, gene_column, version) {
  # Connect to Ensembl BioMart (GRCh37)
  mart <- useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     version = version)
  
  # Retrieve mapping from Ensembl Gene ID to Entrez ID
  mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = df[[gene_column]],
                   mart = mart)
  
  # Merge with the original dataframe to add Entrez IDs
  df <- left_join(df, mapping, by = setNames("ensembl_gene_id", gene_column))
  
  return(df)
}

get_gene_coordinates <- function(df, gene_column, mart = NULL) {
  # If no mart object is provided, use ensembl human genes by default
  if (is.null(mart)) {
    mart <- tryCatch({
      useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    }, error = function(e) {
      message("Could not connect to Ensembl. Trying alternative mirror...")
      useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://uswest.ensembl.org")
    })
  }
  
  # Define attributes to retrieve
  attributes <- c("entrezgene_id", 
                  "chromosome_name", 
                  "start_position", 
                  "end_position", 
                  "strand")
  
  # Get the data from BioMart
  results <- getBM(attributes = attributes,
                   filters = "entrezgene_id",
                   values = df[[gene_column]],
                   mart = mart)
  
  # Clean up the results
  colnames(results) <- c("entrez_id", "chromosome", "start_position", 
                         "end_position", "strand")
  
  # Convert strand to +/- format
  results$strand <- ifelse(results$strand == 1, "+", "-")
  
  # Filter out non-standard chromosomes (optional)
  results <- results[grep("^[0-9XYMT]+$", results$chromosome), ]
  
  return(results)
}