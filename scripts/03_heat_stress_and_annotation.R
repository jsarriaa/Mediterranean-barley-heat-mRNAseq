#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(gprofiler2)
})

repo_root <- normalizePath(".", mustWork = TRUE)
heat_genes_path <- file.path(repo_root, "data/metadata/heat_stress_gene_ids.csv")
ranking_path <- file.path(repo_root, "results/01_deseq2/heat_response_top12_gene_ids.csv")
barlex_path <- file.path(repo_root, "references/BARLEX_MorexV3_annotations.tsv")
out_dir <- file.path(repo_root, "results/03_heat_stress")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(heat_genes_path, ranking_path, barlex_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required heat-stress input(s) are missing: ",
       paste(missing_files, collapse = ", "))
}

heat_genes <- read.csv(heat_genes_path, stringsAsFactors = FALSE)
ranking <- read.csv(ranking_path, stringsAsFactors = FALSE)
if (!"gene_id" %in% names(heat_genes) || anyDuplicated(heat_genes$gene_id)) {
  stop("Heat-gene metadata must contain unique gene_id values.")
}
if (!all(c("rank", "gene_id", "population_mean_sd") %in% names(ranking)) ||
    nrow(ranking) != 12L || anyDuplicated(ranking$gene_id)) {
  stop("Expected the ranked top-12 gene table from the population analysis.")
}

barlex <- read.delim(barlex_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("Gene_ID", "BARLEX_Description", "Pfam", "InterPro") %in% names(barlex))) {
  stop("BARLEX file must contain Gene_ID, BARLEX_Description, Pfam, and InterPro.")
}
barlex <- barlex[!duplicated(barlex$Gene_ID), , drop = FALSE]

message("Converting ", nrow(heat_genes), " heat-response IDs with g:Profiler.")
annotations <- gconvert(
  query = heat_genes$gene_id,
  organism = "hvulgare"
)
if (is.null(annotations)) {
  stop("g:Profiler returned no annotation response. Check network access and retry.")
}
if (!all(c("input", "name", "description") %in% names(annotations))) {
  stop("g:Profiler annotation response lacks input/name/description columns.")
}
annotations <- annotations[!duplicated(annotations$input), , drop = FALSE]

annotation_index <- match(heat_genes$gene_id, annotations$input)
barlex_index <- match(heat_genes$gene_id, barlex$Gene_ID)
result <- data.frame(
  gene_id = heat_genes$gene_id,
  gprofiler_name = annotations$name[annotation_index],
  gprofiler_description = annotations$description[annotation_index],
  BARLEX_Description = barlex$BARLEX_Description[barlex_index],
  Pfam = barlex$Pfam[barlex_index],
  InterPro = barlex$InterPro[barlex_index],
  stringsAsFactors = FALSE
)
has_gprofiler_description <- !is.na(result$gprofiler_description) &
  nzchar(trimws(result$gprofiler_description)) &
  result$gprofiler_description != "None"
has_barlex_description <- !is.na(result$BARLEX_Description) &
  nzchar(trimws(result$BARLEX_Description)) &
  result$BARLEX_Description != "None"
result$combined_description <- NA_character_
result$combined_description[has_gprofiler_description] <-
  result$gprofiler_description[has_gprofiler_description]
use_barlex <- !has_gprofiler_description & has_barlex_description
result$combined_description[use_barlex] <- result$BARLEX_Description[use_barlex]
result$annotation_status <- ifelse(
  !is.na(result$combined_description) & nzchar(trimws(result$combined_description)),
  "annotated",
  "no description returned"
)
write.csv(result, file.path(out_dir, "heat_response_46_annotations.csv"),
          row.names = FALSE, na = "")

top12_index <- match(ranking$gene_id, result$gene_id)
if (anyNA(top12_index)) {
  stop("One or more top-12 genes are absent from the 46-gene annotation table.")
}
top12 <- cbind(ranking, result[top12_index, -1, drop = FALSE])
write.csv(top12, file.path(out_dir, "heat_response_top12_annotations.csv"),
          row.names = FALSE, na = "")
write.csv(
  data.frame(
    metric = c("heat_response_genes", "top12_genes",
               "genes_with_description", "genes_without_description"),
    value = c(nrow(result), nrow(top12),
              sum(result$annotation_status == "annotated"),
              sum(result$annotation_status != "annotated"))
  ),
  file.path(out_dir, "annotation_summary.csv"),
  row.names = FALSE
)
write.csv(
  result[result$annotation_status != "annotated", "gene_id", drop = FALSE],
  file.path(out_dir, "genes_without_description.csv"),
  row.names = FALSE
)
writeLines(
  c(
    paste("query_date_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("organism: hvulgare"),
    paste("gprofiler2_version:", as.character(packageVersion("gprofiler2"))),
    paste("barlex_file:", basename(barlex_path))
  ),
  file.path(out_dir, "annotation_query_metadata.txt")
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
message(
  "Annotation complete: retained all ", nrow(result), " heat genes; ",
  sum(result$annotation_status == "annotated"), " have a description."
)
