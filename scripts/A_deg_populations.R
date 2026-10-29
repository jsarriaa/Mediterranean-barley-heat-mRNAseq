# Part A: genes differentially expressed between the 5 populations (HC genes).
# Run from the repository root:  Rscript scripts/A_deg_populations.R

library(DESeq2)
library(ggplot2)
library(pheatmap)

out <- "results/A_deg"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# 1. Count matrix: HC genes x 60 libraries.
# kallisto counts are estimates with decimals; DESeq2 needs whole numbers,
# so they are rounded (the standard step for kallisto counts in DESeq2).
tab <- read.delim("data/input/kallisto.merged.gene_counts.HC.tsv", check.names = FALSE)
counts <- as.matrix(tab[, -(1:2)])
rownames(counts) <- tab$gene_id
counts <- round(counts)

# 2. Library table: genotype from the library ID (T1M2 -> T1), population
# from samples_metadata.tsv (row n = genotype Tn, column K17_pop).
meta <- read.delim("data/input/samples_metadata.tsv", check.names = FALSE)
samples <- data.frame(library = colnames(counts))
samples$genotype <- sub("M[0-9]+$", "", samples$library)
genotype_number <- as.integer(sub("^T", "", samples$genotype))
samples$population <- paste0("k", meta$K17_pop[match(genotype_number, meta[[1]])])
stopifnot(!anyNA(match(genotype_number, meta[[1]])))
samples$population <- factor(samples$population, levels = c("k1", "k5", "k10", "k12", "k14"))
rownames(samples) <- samples$library
write.table(samples, file.path(out, "library_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(samples$population))   # expect 12 libraries per population

# 3. Sum the 3 plants of each genotype: 20 genotypes, 4 per population.
dds_lib <- DESeqDataSetFromMatrix(counts, samples, design = ~ population)
dds <- collapseReplicates(dds_lib, groupby = dds_lib$genotype)
stopifnot(sum(counts(dds)) == sum(counts(dds_lib)))
print(table(dds$population))       # expect 4 genotypes per population

# 4. Omnibus test: does expression differ anywhere among the 5 populations?
# Likelihood-ratio test of ~ population against ~ 1. DEG: adjusted p < 0.05.
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res <- results(dds_lrt, alpha = 0.05)
res <- data.frame(gene_id = rownames(res), as.data.frame(res))
res <- res[order(res$padj), ]
write.table(res, file.path(out, "omnibus_LRT_all_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
deg <- res[!is.na(res$padj) & res$padj < 0.05, ]
write.table(deg, file.path(out, "omnibus_LRT_DEGs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Note: in an LRT table, log2FoldChange is only the last level vs the first
# (k14 vs k1); the test itself covers all 5 populations.

# Genotype-level expression for pattern analyses (Part B): variance-stabilised
# counts computed with the design (blind = FALSE), as DESeq2 recommends for
# downstream use.
vst_design <- assay(vst(dds_lrt, blind = FALSE))
write.table(data.frame(gene_id = rownames(vst_design), vst_design, check.names = FALSE),
            file.path(out, "vst_genotypes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 5. Pairwise tests between populations (Wald test), reported for reference.
dds_wald <- DESeq(dds)
populations <- levels(dds$population)
summary_rows <- data.frame()
pairwise_degs <- data.frame()   # one row per DEG per contrast, kept for the UpSet plot
for (i in 1:4) {
  for (j in (i + 1):5) {
    a <- populations[i]
    b <- populations[j]
    r <- results(dds_wald, contrast = c("population", b, a), alpha = 0.05)
    r <- data.frame(gene_id = rownames(r), as.data.frame(r))
    r <- r[order(r$padj), ]
    write.table(r, file.path(out, paste0("pairwise_", b, "_vs_", a, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    sig <- !is.na(r$padj) & r$padj < 0.05
    pairwise_degs <- rbind(pairwise_degs, data.frame(
      gene_id = r$gene_id[sig],
      contrast = paste0(b, "_vs_", a),
      direction = ifelse(r$log2FoldChange[sig] > 0, paste("up in", b), paste("up in", a))))
    summary_rows <- rbind(summary_rows, data.frame(
      contrast = paste0(b, "_vs_", a),
      DEGs = sum(sig),
      up_in_first = sum(sig & r$log2FoldChange > 0),
      down_in_first = sum(sig & r$log2FoldChange < 0)))
  }
}
write.table(summary_rows, file.path(out, "pairwise_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 6. Summary of the omnibus test.
summary_omnibus <- data.frame(
  metric = c("HC genes in matrix", "genes with an adjusted p-value", "omnibus DEGs (padj < 0.05)"),
  value = c(nrow(res), sum(!is.na(res$padj)), nrow(deg)))
write.table(summary_omnibus, file.path(out, "omnibus_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary_omnibus)
print(summary_rows)

# 7. Quality control (descriptive only).
# Variance-stabilised expression computed blind to the design, as DESeq2
# recommends for quality control.
vsd_genotype <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd_genotype, intgroup = "population", returnData = TRUE)
percent <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(x = PC1, y = PC2, colour = population, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.9, size = 3, show.legend = FALSE) +   # genotype ID next to each dot
  labs(x = paste0("PC1: ", percent[1], "% variance"), y = paste0("PC2: ", percent[2], "% variance"),
       title = "PCA, 20 genotypes (HC genes)") +
  theme_bw()
ggsave(file.path(out, "pca_genotypes.png"), p, width = 7, height = 5, dpi = 200)

vsd_library <- vst(dds_lib, blind = TRUE)
p <- plotPCA(vsd_library, intgroup = "population") +
  ggtitle("PCA, 60 libraries (HC genes)") + theme_bw()
ggsave(file.path(out, "pca_libraries.png"), p, width = 7, height = 5, dpi = 200)

# Distances between all 60 libraries. No clustering: libraries are ordered by
# population, then genotype, then plant, so the 3 plants of each genotype sit
# together; white gaps separate genotypes.
order_libs <- order(samples$population, as.integer(sub("^T", "", samples$genotype)), samples$library)
libs <- samples$library[order_libs]
distances <- as.matrix(dist(t(assay(vsd_library))))[libs, libs]
annotation <- data.frame(population = samples$population[order_libs], row.names = libs)
genotype_breaks <- which(head(samples$genotype[order_libs], -1) != tail(samples$genotype[order_libs], -1))
pheatmap(distances, annotation_row = annotation, annotation_col = annotation,
         cluster_rows = FALSE, cluster_cols = FALSE,
         gaps_row = genotype_breaks, gaps_col = genotype_breaks,
         fontsize = 6, filename = file.path(out, "library_distances.png"),
         width = 10, height = 9)

# Distances between the 20 genotypes (3 plants summed, as in the tests).
# Genotypes are clustered on these distances (complete linkage); the trees
# show how the genotypes group.
genotype_order <- paste0("T", 1:20)
genotype_population <- dds$population[match(genotype_order, colnames(dds))]
distances_genotype <- as.matrix(dist(t(assay(vsd_genotype))))[genotype_order, genotype_order]
annotation_genotype <- data.frame(population = genotype_population, row.names = genotype_order)
pheatmap(distances_genotype, annotation_row = annotation_genotype, annotation_col = annotation_genotype,
         clustering_distance_rows = as.dist(distances_genotype),
         clustering_distance_cols = as.dist(distances_genotype),
         clustering_method = "complete",
         fontsize = 9, filename = file.path(out, "genotype_distances.png"),
         width = 7.5, height = 6.5)

# 8. UpSet plot of the pairwise DEG sets (each set = one contrast).
# Top: the 30 largest intersections. Middle: which contrasts form each
# intersection. Right: total DEGs per contrast, coloured by direction.
library(gridExtra)
write.table(pairwise_degs, file.path(out, "pairwise_DEGs_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
contrasts <- unique(pairwise_degs$contrast)
membership <- table(pairwise_degs$gene_id, factor(pairwise_degs$contrast, levels = contrasts)) > 0
pattern <- apply(membership, 1, function(x) paste(contrasts[x], collapse = " & "))
inter <- as.data.frame(sort(table(pattern), decreasing = TRUE)[1:30])
names(inter) <- c("pattern", "genes")
inter$id <- seq_len(nrow(inter))
write.table(inter[, c("id", "pattern", "genes")], file.path(out, "upset_top30_intersections.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

dots <- expand.grid(id = inter$id, contrast = contrasts)
dots$member <- mapply(function(id, set) set %in% strsplit(as.character(inter$pattern[id]), " & ")[[1]],
                      dots$id, as.character(dots$contrast))
dots$contrast <- factor(dots$contrast, levels = rev(contrasts))

sizes <- as.data.frame(table(contrast = pairwise_degs$contrast,
                             direction = ifelse(sub("up in ", "", pairwise_degs$direction) ==
                                                sub("_vs_.*", "", pairwise_degs$contrast),
                                                "up in first population", "up in second population")))
sizes$contrast <- factor(sizes$contrast, levels = rev(contrasts))

x_scale <- scale_x_continuous(limits = c(0.5, 30.5), expand = c(0, 0))
y_scale <- scale_y_discrete(limits = rev(contrasts))   # same row order in matrix and set bars
p_top <- ggplot(inter, aes(x = id, y = genes)) +
  geom_col(fill = "grey30") +
  geom_text(aes(label = genes), vjust = -0.3, size = 2.2) +
  x_scale + labs(x = NULL, y = "Genes in intersection") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
p_mid <- ggplot(dots, aes(x = id, y = contrast)) +
  geom_line(data = dots[dots$member, ], aes(group = id), colour = "grey30") +
  geom_point(aes(colour = member), size = 2.2, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "grey20", `FALSE` = "grey85")) +
  x_scale + y_scale + labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
p_right <- ggplot(sizes, aes(x = Freq, y = contrast, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("up in first population" = "#D55E00", "up in second population" = "#0072B2")) +
  y_scale + labs(x = "DEGs per contrast", y = NULL, fill = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.text = element_text(size = 6), legend.key.size = unit(0.3, "cm"))

# Line up the panels: same widths for top and middle, same heights for middle and right.
g_top <- ggplotGrob(p_top); g_mid <- ggplotGrob(p_mid); g_right <- ggplotGrob(p_right)
widths <- grid::unit.pmax(g_top$widths, g_mid$widths)
g_top$widths <- widths; g_mid$widths <- widths
heights <- grid::unit.pmax(g_mid$heights, g_right$heights)
g_mid$heights <- heights; g_right$heights <- heights
upset <- arrangeGrob(g_top, grid::nullGrob(), g_mid, g_right, ncol = 2,
                     widths = c(3, 1), heights = c(1.2, 1))
ggsave(file.path(out, "pairwise_upset.png"), upset, width = 12, height = 7, dpi = 200)

write.table(data.frame(genotype = colnames(dds), size_factor = sizeFactors(dds_lrt)),
            file.path(out, "size_factors_genotypes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
