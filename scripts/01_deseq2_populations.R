#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

repo_root <- normalizePath(".", mustWork = TRUE)
count_path <- file.path(repo_root, "data/input/kallisto.merged.gene_counts.tsv")
metadata_path <- file.path(repo_root, "data/metadata/samples.csv")
heat_genes_path <- file.path(repo_root, "data/metadata/heat_stress_gene_ids.csv")
out_dir <- file.path(repo_root, "results/01_deseq2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(count_path, metadata_path, heat_genes_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required input file(s) are missing: ", paste(missing_files, collapse = ", "))
}

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample", "subpopulation") %in% names(metadata))) {
  stop("samples.csv must contain columns named 'sample' and 'subpopulation'.")
}
if (anyNA(metadata$sample) || any(!nzchar(metadata$sample)) ||
    anyDuplicated(metadata$sample)) {
  stop("Sample IDs must be non-empty and unique.")
}
population_levels <- c("k1", "k5", "k10", "k12", "k14")
if (!setequal(unique(metadata$subpopulation), population_levels)) {
  stop("The sample sheet must contain exactly the five populations: ",
       paste(population_levels, collapse = ", "))
}
if (any(!grepl("^T([1-9]|1[0-9]|20)M[0-9]+$", metadata$sample))) {
  stop("The sample sheet includes an ID outside the T1-T20 analysis scope.")
}
metadata$subpopulation <- factor(metadata$subpopulation, levels = population_levels)
metadata <- metadata[order(metadata$subpopulation, metadata$sample), , drop = FALSE]
rownames(metadata) <- metadata$sample
if (!identical(as.integer(table(metadata$subpopulation)), rep(12L, 5L))) {
  stop("Expected exactly 12 sample libraries in each of k1, k5, k10, k12, and k14.")
}

counts_table <- read.delim(count_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("gene_id", "gene_name") %in% names(counts_table))) {
  stop("Count matrix must begin with 'gene_id' and 'gene_name' columns.")
}
if (anyNA(counts_table$gene_id) || anyDuplicated(counts_table$gene_id)) {
  stop("Count-matrix gene IDs must be present and unique.")
}
matrix_samples <- setdiff(names(counts_table), c("gene_id", "gene_name"))
target_sample_pattern <- "^T([1-9]|1[0-9]|20)M[0-9]+$"
target_matrix_samples <- matrix_samples[grepl(target_sample_pattern, matrix_samples)]
if (!setequal(target_matrix_samples, metadata$sample)) {
  missing_columns <- setdiff(metadata$sample, matrix_samples)
  extra_columns <- setdiff(target_matrix_samples, metadata$sample)
  stop(
    "Count-matrix samples do not exactly match the T1-T20 metadata. Missing: ",
    paste(missing_columns, collapse = ", "),
    "; unlisted T1-T20 samples: ", paste(extra_columns, collapse = ", ")
  )
}

counts <- as.matrix(counts_table[, metadata$sample, drop = FALSE])
storage.mode(counts) <- "numeric"
if (any(!is.finite(counts)) || any(counts < 0)) {
  stop("Count data must be finite and non-negative.")
}
counts <- round(counts)
storage.mode(counts) <- "integer"
rownames(counts) <- counts_table$gene_id

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ subpopulation
)

dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)
lrt <- as.data.frame(results(dds_lrt, alpha = 0.05))
lrt$gene_id <- rownames(lrt)
lrt <- lrt[order(lrt$padj, na.last = TRUE), ]
write.csv(lrt, file.path(out_dir, "population_omnibus_LRT_all_genes.csv"), row.names = FALSE)

significant_lrt <- lrt[!is.na(lrt$padj) & lrt$padj < 0.05, , drop = FALSE]
write.csv(
  significant_lrt,
  file.path(out_dir, "population_omnibus_LRT_significant_genes.csv"),
  row.names = FALSE
)

dds_wald <- DESeq(dds, test = "Wald", quiet = TRUE)
pairwise_dir <- file.path(out_dir, "pairwise")
dir.create(pairwise_dir, recursive = TRUE, showWarnings = FALSE)
pairwise_summary <- list()
contrast_pairs <- combn(population_levels, 2, simplify = FALSE)

for (pair in contrast_pairs) {
  reference <- pair[[1]]
  comparison <- pair[[2]]
  contrast_name <- paste(comparison, "vs", reference, sep = "_")
  result <- as.data.frame(results(
    dds_wald,
    contrast = c("subpopulation", comparison, reference),
    alpha = 0.05
  ))
  result$gene_id <- rownames(result)
  result <- result[order(result$padj, na.last = TRUE), ]
  write.csv(result, file.path(pairwise_dir, paste0(contrast_name, "_all_genes.csv")),
            row.names = FALSE)
  significant <- result[!is.na(result$padj) & result$padj < 0.05, , drop = FALSE]
  write.csv(significant,
            file.path(pairwise_dir, paste0(contrast_name, "_significant_genes.csv")),
            row.names = FALSE)
  pairwise_summary[[contrast_name]] <- data.frame(
    contrast = contrast_name,
    reference_population = reference,
    comparison_population = comparison,
    significant_genes = nrow(significant),
    stringsAsFactors = FALSE
  )
}

pairwise_summary <- do.call(rbind, pairwise_summary)
write.csv(pairwise_summary, file.path(out_dir, "pairwise_contrast_summary.csv"),
          row.names = FALSE)

vsd <- varianceStabilizingTransformation(dds_lrt, blind = FALSE)
sample_order <- rownames(metadata)
vst_matrix <- assay(vsd)[, sample_order, drop = FALSE]
annotation <- data.frame(Subpopulation = metadata$subpopulation)
rownames(annotation) <- rownames(metadata)

pca <- plotPCA(vsd, intgroup = "subpopulation") +
  theme_minimal(base_size = 12) +
  labs(title = "Transcriptomic profiles across five populations",
       x = "PC1", y = "PC2", color = "Population")
ggsave(file.path(out_dir, "pca_five_populations.png"), pca,
       width = 8, height = 6, dpi = 180)
ggsave(file.path(out_dir, "pca_five_populations.pdf"), pca,
       width = 8, height = 6)

plot_heatmap <- function(gene_ids, filename, title, width = 10, height = 8) {
  gene_ids <- intersect(gene_ids, rownames(vst_matrix))
  if (length(gene_ids) < 2L) {
    stop("At least two valid genes are required for heatmap: ", filename)
  }
  matrix <- vst_matrix[gene_ids, , drop = FALSE]
  matrix <- matrix - rowMeans(matrix)
  pheatmap(
    matrix,
    annotation_col = annotation,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = length(gene_ids) <= 50L,
    fontsize_row = 6,
    gaps_col = cumsum(as.integer(table(metadata$subpopulation)))[-length(population_levels)],
    main = title,
    filename = file.path(out_dir, filename),
    width = width,
    height = height
  )
}

top_genes <- head(significant_lrt$gene_id, 100L)
plot_heatmap(top_genes, "top100_population_heatmap.png",
             "Top 100 genes varying across populations", height = 10)

heat_gene_table <- read.csv(heat_genes_path, stringsAsFactors = FALSE)
if (!"gene_id" %in% names(heat_gene_table) || anyDuplicated(heat_gene_table$gene_id)) {
  stop("heat_stress_gene_ids.csv must contain unique gene_id values.")
}
missing_heat_genes <- setdiff(heat_gene_table$gene_id, rownames(vst_matrix))
if (length(missing_heat_genes)) {
  stop("Heat-response IDs absent from count matrix: ",
       paste(missing_heat_genes, collapse = ", "))
}
plot_heatmap(heat_gene_table$gene_id, "heat_response_46_populations.png",
             "Expression of 46 GO:0009408 heat-response genes")

group_means <- vapply(
  population_levels,
  function(population) {
    rowMeans(vst_matrix[, metadata$subpopulation == population, drop = FALSE])
  },
  numeric(nrow(vst_matrix))
)
rownames(group_means) <- rownames(vst_matrix)
candidate_means <- group_means[heat_gene_table$gene_id, , drop = FALSE]
population_mean_sd <- apply(candidate_means, 1, sd)
candidate_ranking <- data.frame(
  gene_id = names(population_mean_sd),
  population_mean_sd = unname(population_mean_sd),
  stringsAsFactors = FALSE
)
candidate_ranking <- candidate_ranking[
  order(candidate_ranking$population_mean_sd, decreasing = TRUE), , drop = FALSE
]
candidate_ranking$rank <- seq_len(nrow(candidate_ranking))
candidate_ranking <- candidate_ranking[, c("rank", "gene_id", "population_mean_sd")]
write.csv(candidate_ranking,
          file.path(out_dir, "heat_response_46_population_variation_ranking.csv"),
          row.names = FALSE)
top12_genes <- head(candidate_ranking$gene_id, 12L)
write.csv(candidate_ranking[seq_along(top12_genes), , drop = FALSE],
          file.path(out_dir, "heat_response_top12_gene_ids.csv"),
          row.names = FALSE)
plot_heatmap(top12_genes, "heat_response_top12_populations.png",
             "Top 12 heat-response genes by population-mean variation")

write.csv(
  data.frame(
    metric = c("samples", "genes_in_matrix", "genes_tested_by_LRT",
               "significant_LRT_genes", "heat_response_genes", "populations"),
    value = c(ncol(counts), nrow(counts), sum(!is.na(lrt$pvalue)),
              nrow(significant_lrt), nrow(heat_gene_table),
              length(population_levels)),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "analysis_summary.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
message(
  "Population analysis complete: ", ncol(counts), " samples; ",
  nrow(significant_lrt), " significant omnibus LRT genes; ",
  length(top12_genes), " ranked heat-response candidates."
)
