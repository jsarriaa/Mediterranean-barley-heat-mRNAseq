#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(gprofiler2)
  library(ggplot2)
})

repo_root <- normalizePath(".", mustWork = TRUE)
de_path <- file.path(repo_root, "results/01_deseq2/population_omnibus_LRT_significant_genes.csv")
out_dir <- file.path(repo_root, "results/02_go_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(de_path)) {
  stop("Run scripts/01_deseq2_populations.R before GO enrichment.")
}

de <- read.csv(de_path, stringsAsFactors = FALSE)
if (!all(c("gene_id", "padj") %in% names(de)) || !nrow(de)) {
  stop("Omnibus DE file is empty or lacks gene_id/padj columns.")
}
query_genes <- unique(de$gene_id[!is.na(de$gene_id) & nzchar(de$gene_id)])
if (!length(query_genes)) {
  stop("No valid gene IDs are available for enrichment.")
}

message("Querying g:Profiler for ", length(query_genes), " barley genes.")
go <- gost(
  query = query_genes,
  organism = "hvulgare",
  significant = TRUE,
  evcodes = TRUE
)
if (is.null(go)) {
  stop("g:Profiler returned no response. Check network access and try again.")
}

if (is.null(go$result) || !nrow(go$result)) {
  write.csv(data.frame(), file.path(out_dir, "GO_Enrichment_Table.csv"),
            row.names = FALSE)
  message("No significant GO terms were returned for this query.")
} else {
  go_table <- go$result
  for (column in names(go_table)) {
    if (is.list(go_table[[column]])) {
      go_table[[column]] <- vapply(
        go_table[[column]],
        function(value) paste(as.character(value), collapse = ","),
        character(1)
      )
    }
  }
  write.csv(go_table, file.path(out_dir, "GO_Enrichment_Table.csv"),
            row.names = FALSE)
  plot <- gostplot(go, capped = FALSE, interactive = FALSE)
  ggsave(file.path(out_dir, "GO_Enrichment_Plot.pdf"), plot,
         width = 12, height = 10)
  ggsave(file.path(out_dir, "GO_Enrichment_Plot.png"), plot,
         width = 12, height = 10, dpi = 180)
}

heat_term <- if (!is.null(go$result) && nrow(go$result)) {
  go_table[go_table$term_id == "GO:0009408", , drop = FALSE]
} else {
  data.frame()
}
curated_genes <- read.csv(
  file.path(repo_root, "data/metadata/heat_stress_gene_ids.csv"),
  stringsAsFactors = FALSE
)$gene_id
if (nrow(heat_term)) {
  heat_gene_ids <- unique(trimws(unlist(strsplit(
    as.character(heat_term$intersection[[1]]), ",", fixed = TRUE
  ))))
  heat_gene_ids <- heat_gene_ids[nzchar(heat_gene_ids)]
  write.csv(data.frame(gene_id = heat_gene_ids),
            file.path(out_dir, "GO_0009408_heat_response_genes.csv"),
            row.names = FALSE)
  write.csv(
    data.frame(
      term_id = heat_term$term_id[[1]],
      term_name = heat_term$term_name[[1]],
      p_value = heat_term$p_value[[1]],
      intersection_size = length(heat_gene_ids),
      curated_overlap = length(intersect(heat_gene_ids, curated_genes)),
      same_gene_set = setequal(heat_gene_ids, curated_genes),
      query_date_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    file.path(out_dir, "heat_response_term_summary.csv"),
    row.names = FALSE
  )
  if ("GO:0009408" %in% heat_term$term_id) {
    original_directory <- getwd()
    tryCatch(
      {
        setwd(out_dir)
        publish_gostplot(
          plot,
          highlight_terms = "GO:0009408",
          width = 10,
          height = 8,
          filename = "GO_0009408_highlighted.pdf"
        )
      },
      finally = setwd(original_directory)
    )
  }
  if (!setequal(heat_gene_ids, curated_genes)) {
    warning(
      "The current g:Profiler GO:0009408 intersection differs from the ",
      "versioned 46-gene input; see heat_response_term_summary.csv."
    )
  }
  message("GO:0009408 intersection: ", length(heat_gene_ids), " genes.")
} else {
  message("GO:0009408 was not significant in this live query; the archived 46-gene set remains unchanged.")
}

writeLines(
  c(
    paste("query_date_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("organism: hvulgare"),
    paste("query_gene_count:", length(query_genes)),
    paste("significant_terms:", if (is.null(go$result)) 0L else nrow(go$result)),
    paste("gprofiler2_version:", as.character(packageVersion("gprofiler2")))
  ),
  file.path(out_dir, "query_metadata.txt")
)
