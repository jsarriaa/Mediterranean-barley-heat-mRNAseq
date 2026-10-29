# Part B: GO enrichment of groups of population DEGs.
#   1. Gene-to-GO table for HC genes (Morex V3 + Ensembl 63).
#   2. DEG groups:
#      - expression-pattern clusters of the omnibus DEGs (hierarchical
#        clustering, number of clusters by silhouette width);
#      - per pairwise contrast, DEGs with |log2FC| >= 1, up and down separately.
#   3. Each group tested for GO enrichment: hypergeometric test, background =
#      tested HC genes with GO in the namespace, BH within each namespace.
# Run from the repository root:  Rscript scripts/B_go_enrichment.R

library(ggplot2)

ann <- "references/annotation"
out <- "results/B_go"
dir.create(file.path(out, "clusters"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "contrasts"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Gene-to-GO table
# ---------------------------------------------------------------------------

# GO terms from go-basic.obo: one block per term.
lines <- readLines(file.path(ann, "go-basic.obo"))
block <- cumsum(lines == "[Term]" | lines == "[Typedef]")
in_term <- block %in% unique(block[lines == "[Term]"])
field <- function(prefix) {
  i <- which(in_term & startsWith(lines, prefix))
  data.frame(block = block[i], value = substring(lines[i], nchar(prefix) + 1))
}
ids <- field("id: ")
terms <- data.frame(block = ids$block, go_id = ids$value)
f <- field("name: ");        terms$name <- f$value[match(terms$block, f$block)]
f <- field("namespace: ");   terms$namespace <- f$value[match(terms$block, f$block)]
f <- field("is_obsolete: "); terms$obsolete <- terms$block %in% f$block[f$value == "true"]
f <- field("replaced_by: "); terms$replaced_by <- f$value[match(terms$block, f$block)]
alt <- field("alt_id: ")
alt$go_id <- terms$go_id[match(alt$block, terms$block)]
parents <- rbind(field("is_a: "), field("relationship: part_of "))
edges <- data.frame(child = terms$go_id[match(parents$block, terms$block)],
                    parent = substring(parents$value, 1, 10))

# GO assignments of HC genes from both sources.
confidence <- read.delim("references/MorexV3_gene_confidence.tsv")
hc <- confidence$gene_id[confidence$confidence == "high"]
morex <- read.delim(file.path(ann, "morexv3_go.tsv"), header = FALSE, col.names = c("gene_id", "go_id"))
morex$source <- "MorexV3"
morex$evidence <- "not given"
ensembl <- read.delim(file.path(ann, "ensembl63_go.tsv"), header = FALSE,
                      col.names = c("gene_id", "go_id", "evidence", "namespace"))
ensembl$source <- "Ensembl63"
gene2go <- rbind(morex[, c("gene_id", "go_id", "source", "evidence")],
                 ensembl[, c("gene_id", "go_id", "source", "evidence")])
gene2go <- gene2go[gene2go$gene_id %in% hc, ]

# Old GO IDs: alternative IDs -> main ID; obsolete terms -> replacement or dropped.
is_alt <- gene2go$go_id %in% alt$value
gene2go$go_id[is_alt] <- alt$go_id[match(gene2go$go_id[is_alt], alt$value)]
obsolete <- gene2go$go_id %in% terms$go_id[terms$obsolete]
replacement <- terms$replaced_by[match(gene2go$go_id, terms$go_id)]
n_replaced <- sum(obsolete & !is.na(replacement))
n_obsolete_dropped <- sum(obsolete & is.na(replacement))
gene2go$go_id[obsolete & !is.na(replacement)] <- replacement[obsolete & !is.na(replacement)]
gene2go <- gene2go[gene2go$go_id %in% terms$go_id[!terms$obsolete], ]

# Union: one row per gene and GO term, recording which source(s) gave it.
key <- paste(gene2go$gene_id, gene2go$go_id)
union <- data.frame(gene_id = gene2go$gene_id, go_id = gene2go$go_id)[!duplicated(key), ]
ukey <- paste(union$gene_id, union$go_id)
union$source <- tapply(gene2go$source, key, function(s) paste(sort(unique(s)), collapse = "+"))[ukey]
union$evidence <- tapply(gene2go$evidence, key, function(s) paste(sort(unique(s)), collapse = "+"))[ukey]
union$go_name <- terms$name[match(union$go_id, terms$go_id)]
union$namespace <- terms$namespace[match(union$go_id, terms$go_id)]
write.table(union, file.path(ann, "gene2go_HC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

coverage <- data.frame(
  metric = c("HC genes", "HC genes with any GO", "HC genes without GO",
             "obsolete GO assignments replaced", "obsolete GO assignments dropped (no replacement)",
             "gene-GO pairs (union)"),
  value = c(length(hc), length(unique(union$gene_id)), length(setdiff(hc, union$gene_id)),
            n_replaced, n_obsolete_dropped, nrow(union)))
write.table(coverage, file.path(out, "go_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Pass each annotation up to all parent terms: a gene annotated to a term
# also counts for every ancestor of that term.
ancestors <- data.frame(term = terms$go_id, ancestor = terms$go_id)
repeat {
  step <- merge(ancestors, edges, by.x = "ancestor", by.y = "child")
  step <- data.frame(term = step$term, ancestor = step$parent)
  before <- nrow(ancestors)
  ancestors <- unique(rbind(ancestors, step))
  if (nrow(ancestors) == before) break
}
propagated <- merge(union[, c("gene_id", "go_id")], ancestors, by.x = "go_id", by.y = "term")
propagated <- unique(data.frame(gene_id = propagated$gene_id, go_id = propagated$ancestor))
propagated$namespace <- terms$namespace[match(propagated$go_id, terms$go_id)]

# Background = HC genes tested in Part A; omnibus DEGs = adjusted p < 0.05.
lrt <- read.delim("results/A_deg/omnibus_LRT_all_genes.tsv")
universe <- lrt$gene_id[!is.na(lrt$padj)]
degs <- lrt$gene_id[!is.na(lrt$padj) & lrt$padj < 0.05]
propagated <- propagated[propagated$gene_id %in% universe, ]

# ---------------------------------------------------------------------------
# Enrichment test and plot for one gene list
# ---------------------------------------------------------------------------

# Per namespace: N = background genes with GO in the namespace, n = list genes
# among them, K = background genes in the term, k = list genes in the term.
# p = probability of k or more list genes in the term by chance.
enrich <- function(query) {
  results <- data.frame()
  for (ns in c("molecular_function", "biological_process", "cellular_component")) {
    p_ns <- propagated[propagated$namespace == ns, ]
    N <- length(unique(p_ns$gene_id))
    n <- length(intersect(unique(p_ns$gene_id), query))
    K <- tapply(p_ns$gene_id, p_ns$go_id, length)
    k <- tapply(p_ns$gene_id %in% query, p_ns$go_id, sum)
    genes <- tapply(p_ns$gene_id, p_ns$go_id, function(g) paste(sort(g[g %in% query]), collapse = ","))
    res <- data.frame(namespace = ns, go_id = names(K), go_name = terms$name[match(names(K), terms$go_id)],
                      term_size = as.integer(K), genes_in_term = as.integer(k[names(K)]),
                      list_annotated = n, universe_annotated = N)
    res$p_value <- phyper(res$genes_in_term - 1, res$term_size, N - res$term_size, n, lower.tail = FALSE)
    res$p_adj <- p.adjust(res$p_value, method = "BH")
    res$genes <- genes[names(K)]
    results <- rbind(results, res)
  }
  results[order(results$namespace, results$p_adj), ]
}

# g:Profiler-style plot: one dot per significant term; x = term position by
# GO ID among all tested terms of its namespace; y = -log10(adjusted p),
# capped at 16; dot size = number of background genes in the term.
plot_enrich <- function(results, title, file) {
  d <- results
  d$namespace <- factor(d$namespace, levels = c("molecular_function", "biological_process", "cellular_component"),
                        labels = c("GO:MF", "GO:BP", "GO:CC"))
  d <- d[order(d$namespace, d$go_id), ]
  d$position <- ave(seq_len(nrow(d)), d$namespace, FUN = seq_along)
  d$log_p <- pmin(-log10(d$p_adj), 16)
  ends <- aggregate(position ~ namespace, d, max)
  p <- ggplot(d[d$p_adj < 0.05, ], aes(x = position, y = log_p, colour = namespace, size = term_size)) +
    geom_blank(data = rbind(transform(ends, position = 1), ends), aes(x = position), inherit.aes = FALSE) +
    geom_point(alpha = 0.7) +
    facet_grid(. ~ namespace, scales = "free_x", space = "free_x", switch = "x", drop = FALSE) +
    scale_colour_manual(values = c("GO:MF" = "#dc3912", "GO:BP" = "#ff9900", "GO:CC" = "#109618"), guide = "none") +
    scale_size_area(max_size = 6, name = "Term size") +
    scale_y_continuous(limits = c(0, 16.5), breaks = seq(0, 16, 2), labels = c(seq(0, 14, 2), ">16")) +
    labs(x = NULL, y = expression(-log[10](p[adj])), title = title) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          strip.placement = "outside", strip.background = element_blank())
  ggsave(file, p, width = 10, height = 5, dpi = 200)
}

# Tests one list, writes its significant terms and plot, returns a summary row.
run_list <- function(name, genes, folder) {
  res <- enrich(genes)
  sig <- res[res$p_adj < 0.05, ]
  write.table(sig, file.path(out, folder, paste0(name, "_GO_significant.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  plot_enrich(res, paste0(name, " (", length(genes), " genes)"),
              file.path(out, folder, paste0(name, "_GO.png")))
  data.frame(list = name, genes = length(genes),
             significant_MF = sum(sig$namespace == "molecular_function"),
             significant_BP = sum(sig$namespace == "biological_process"),
             significant_CC = sum(sig$namespace == "cellular_component"))
}

# ---------------------------------------------------------------------------
# 2a. DEG groups by expression pattern (clusters)
# ---------------------------------------------------------------------------

# Mean expression per population (4 genotypes each), then each gene is centred
# and scaled across the 5 populations, so clusters group shapes, not levels.
vst <- read.delim("results/A_deg/vst_genotypes.tsv", check.names = FALSE)
rownames(vst) <- vst$gene_id
library_table <- read.delim("results/A_deg/library_table.tsv")
genotype_population <- unique(library_table[, c("genotype", "population")])
populations <- c("k1", "k5", "k10", "k12", "k14")
means <- sapply(populations, function(p)
  rowMeans(vst[degs, genotype_population$genotype[genotype_population$population == p]]))
z <- t(scale(t(means)))

# Hierarchical clustering: Euclidean distance, Ward linkage.
distances <- dist(z)
tree <- hclust(distances, method = "ward.D2")
D <- as.matrix(distances)

# Mean silhouette width: for each gene, a = mean distance to its own cluster,
# b = mean distance to the nearest other cluster, s = (b - a) / max(a, b).
mean_silhouette <- function(cluster) {
  k <- max(cluster)
  size <- tabulate(cluster, k)
  sums <- sapply(1:k, function(j) rowSums(D[, cluster == j, drop = FALSE]))
  own <- cbind(seq_along(cluster), cluster)
  a <- sums[own] / pmax(size[cluster] - 1, 1)
  others <- sweep(sums, 2, size, "/")
  others[own] <- Inf
  b <- apply(others, 1, min)
  s <- (b - a) / pmax(a, b)
  s[size[cluster] == 1] <- 0
  mean(s)
}
k_values <- 2:10
silhouette <- data.frame(k = k_values,
                         mean_silhouette = sapply(k_values, function(k) mean_silhouette(cutree(tree, k))))
k_best <- silhouette$k[which.max(silhouette$mean_silhouette)]
write.table(silhouette, file.path(out, "clusters", "silhouette.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p <- ggplot(silhouette, aes(x = k, y = mean_silhouette)) +
  geom_line() + geom_point(size = 2) +
  geom_point(data = silhouette[silhouette$k == k_best, ], colour = "red", size = 3) +
  scale_x_continuous(breaks = k_values) +
  labs(x = "Number of clusters (k)", y = "Mean silhouette width", title = "Choice of cluster number (red = chosen)") +
  theme_bw(base_size = 10)
ggsave(file.path(out, "clusters", "silhouette.png"), p, width = 6, height = 4, dpi = 200)

cluster <- cutree(tree, k_best)
assignment <- data.frame(gene_id = degs, cluster = cluster, round(z, 4), check.names = FALSE)
write.table(assignment, file.path(out, "clusters", "cluster_assignment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Pattern plot: one thin line per gene, thick red line = cluster mean.
long <- data.frame(gene_id = rep(degs, times = 5),
                   population = factor(rep(populations, each = length(degs)), levels = populations),
                   z = as.vector(z),
                   cluster = rep(cluster, times = 5))
long$panel <- paste0("Cluster ", long$cluster, " (", table(cluster)[as.character(long$cluster)], " genes)")
cluster_means <- aggregate(z ~ panel + population, long, mean)
p <- ggplot(long, aes(x = population, y = z, group = gene_id)) +
  geom_line(alpha = 0.03, colour = "grey40") +
  geom_line(data = cluster_means, aes(group = panel), colour = "red", linewidth = 1) +
  facet_wrap(~ panel) +
  labs(x = "Population", y = "Scaled expression (z-score of population means)") +
  theme_bw(base_size = 10)
ggsave(file.path(out, "clusters", "cluster_patterns.png"), p, width = 10, height = 7, dpi = 200)

summary_rows <- data.frame()
for (j in sort(unique(cluster))) {
  summary_rows <- rbind(summary_rows, run_list(paste0("cluster_", j), degs[cluster == j], "clusters"))
}

# ---------------------------------------------------------------------------
# 2b. DEG groups by pairwise contrast: |log2FC| >= 1, up and down separately
# ---------------------------------------------------------------------------

for (file in list.files("results/A_deg", pattern = "^pairwise_k.*_vs_k.*\\.tsv$", full.names = TRUE)) {
  contrast <- sub("^pairwise_(.*)\\.tsv$", "\\1", basename(file))
  first <- sub("_vs_.*", "", contrast)
  second <- sub(".*_vs_", "", contrast)
  r <- read.delim(file)
  ok <- !is.na(r$padj) & r$padj < 0.05
  up <- r$gene_id[ok & r$log2FoldChange >= 1]      # higher in the first population
  down <- r$gene_id[ok & r$log2FoldChange <= -1]   # higher in the second population
  if (length(up) > 0) summary_rows <- rbind(summary_rows, run_list(paste0(contrast, "_up_in_", first), up, "contrasts"))
  if (length(down) > 0) summary_rows <- rbind(summary_rows, run_list(paste0(contrast, "_up_in_", second), down, "contrasts"))
}

# ---------------------------------------------------------------------------
# 3. Summary
# ---------------------------------------------------------------------------

write.table(summary_rows, file.path(out, "summary_significant_terms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(silhouette)
cat("Chosen number of clusters:", k_best, "\n")
print(summary_rows)
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
