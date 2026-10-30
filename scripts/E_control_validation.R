# Part E: control validation. Genotypes 21 (heat-resistant) and 22
# (heat-susceptible), grown with (T) and without (C) heat; HC genes.
# Run from the repository root:  Rscript scripts/E_control_validation.R
# Needs Part C (results/C_cherrypicking/heat_response_genes.tsv).

library(DESeq2)
library(ggplot2)
library(pheatmap)

out <- "results/E_control"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# 1. Count matrix: HC genes x 15 libraries, rounded for DESeq2 (as in Part A).
tab <- read.delim("data/input/kallisto.E.gene_counts.HC.tsv", check.names = FALSE)
counts <- as.matrix(tab[, -(1:2)])
rownames(counts) <- tab$gene_id
counts <- round(counts)

# 2. Library table from the library ID: C21M2 -> treatment control,
# genotype 21, plant 2.
samples <- data.frame(library = colnames(counts))
samples$treatment <- factor(ifelse(startsWith(samples$library, "C"), "control", "heat"),
                            levels = c("control", "heat"))
samples$genotype <- factor(sub("^[CT]([0-9]+)M[0-9]+$", "\\1", samples$library), levels = c("21", "22"))
samples$plant <- sub("^.*M", "M", samples$library)
rownames(samples) <- samples$library
write.table(samples, file.path(out, "library_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(samples$genotype, samples$treatment))   # expect 3/4 (21) and 4/4 (22)

# 3. Quality control (descriptive only): variance-stabilised expression,
# blind to the design, all 15 libraries.
dds_all <- DESeqDataSetFromMatrix(counts, samples, design = ~ 1)
vsd <- vst(dds_all, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("genotype", "treatment"), returnData = TRUE)
percent <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(x = PC1, y = PC2, colour = treatment, shape = genotype, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.9, size = 3, show.legend = FALSE) +   # library ID next to each dot
  scale_colour_manual(values = c(control = "#0072B2", heat = "#D55E00")) +
  labs(x = paste0("PC1: ", percent[1], "% variance"), y = paste0("PC2: ", percent[2], "% variance"),
       title = "PCA, 15 libraries (HC genes)") +
  theme_bw()
ggsave(file.path(out, "pca_libraries.png"), p, width = 7, height = 5, dpi = 200)

# Distances between libraries, clustered (complete linkage on the distances).
distances <- as.matrix(dist(t(assay(vsd))))
annotation <- samples[, c("genotype", "treatment")]
annotation_colours <- list(treatment = c(control = "#0072B2", heat = "#D55E00"),
                           genotype = c(`21` = "grey30", `22` = "grey75"))
pheatmap(distances, annotation_row = annotation, annotation_col = annotation,
         annotation_colors = annotation_colours,
         clustering_distance_rows = as.dist(distances),
         clustering_distance_cols = as.dist(distances),
         clustering_method = "complete",
         fontsize = 9, filename = file.path(out, "library_distances.png"),
         width = 7.5, height = 6.5)

# 4. Heat vs control within each genotype: a separate DESeq2 model per
# genotype, design ~ treatment (control = reference), Wald test; each plant
# is one replicate. DEG: adjusted p < 0.05. log2FoldChange > 0 = higher in heat.
results_list <- list()
summary_rows <- data.frame()
for (g in c("21", "22")) {
  keep <- samples$genotype == g
  dds <- DESeqDataSetFromMatrix(counts[, keep], samples[keep, ], design = ~ treatment)
  dds <- DESeq(dds)
  r <- results(dds, contrast = c("treatment", "heat", "control"), alpha = 0.05)
  r <- data.frame(gene_id = rownames(r), as.data.frame(r))
  r <- r[order(r$padj), ]
  write.table(r, file.path(out, paste0("heat_vs_control_G", g, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  results_list[[g]] <- r
  sig <- !is.na(r$padj) & r$padj < 0.05
  summary_rows <- rbind(summary_rows, data.frame(
    genotype = g, control_plants = sum(keep & samples$treatment == "control"),
    heat_plants = sum(keep & samples$treatment == "heat"),
    genes_tested = sum(!is.na(r$padj)), DEGs = sum(sig),
    up_in_heat = sum(sig & r$log2FoldChange > 0), down_in_heat = sum(sig & r$log2FoldChange < 0)))
}
write.table(summary_rows, file.path(out, "heat_vs_control_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary_rows)

# Volcano plots: one panel per genotype. Each dot is a gene; red = DEG
# higher in heat, blue = DEG lower in heat, grey = not DEG.
volcano <- do.call(rbind, lapply(names(results_list), function(g) {
  r <- results_list[[g]]
  r <- r[!is.na(r$padj), ]
  data.frame(genotype = paste("Genotype", g), log2FoldChange = r$log2FoldChange, padj = r$padj,
             class = ifelse(r$padj >= 0.05, "not DEG", ifelse(r$log2FoldChange > 0, "up in heat", "down in heat")))
}))
p <- ggplot(volcano, aes(x = log2FoldChange, y = -log10(padj), colour = class)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(values = c("up in heat" = "#D55E00", "down in heat" = "#0072B2", "not DEG" = "grey70")) +
  facet_wrap(~ genotype) +
  labs(x = "log2 fold change (heat / control)", y = "-log10 adjusted p", colour = NULL,
       title = "Heat vs control within each genotype (HC genes)") +
  theme_bw()
ggsave(file.path(out, "volcano.png"), p, width = 10, height = 5, dpi = 200)

# 5. The 76 heat-response genes of Part C: log2 fold change and adjusted p in
# each genotype, and a heatmap of their expression in the 15 libraries.
heat_genes <- read.delim("results/C_cherrypicking/heat_response_genes.tsv", quote = "")
heat_table <- heat_genes[, c("gene_id", "description")]
for (g in c("21", "22")) {
  r <- results_list[[g]]
  heat_table[[paste0("log2FC_G", g)]] <- r$log2FoldChange[match(heat_table$gene_id, r$gene_id)]
  heat_table[[paste0("padj_G", g)]] <- r$padj[match(heat_table$gene_id, r$gene_id)]
}
write.table(heat_table, file.path(out, "heat_response_genes_E.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap: variance-stabilised expression (blind, from step 3), scaled per
# gene across the 15 libraries (z-score). Columns: genotype 21 control, 21
# heat, 22 control, 22 heat. Rows clustered; row bars = DEG in each genotype.
lib_order <- samples$library[order(samples$genotype, samples$treatment, samples$plant)]
z <- t(scale(t(assay(vsd)[heat_genes$gene_id, lib_order])))
z[is.na(z)] <- 0   # genes with no variation (not expected; kept as 0)
deg_label <- function(g) {
  p <- heat_table[[paste0("padj_G", g)]]; fc <- heat_table[[paste0("log2FC_G", g)]]
  factor(ifelse(is.na(p) | p >= 0.05, "not DEG", ifelse(fc > 0, "up in heat", "down in heat")),
         levels = c("up in heat", "down in heat", "not DEG"))
}
annotation_row <- data.frame(`DEG genotype 21` = deg_label("21"), `DEG genotype 22` = deg_label("22"),
                             row.names = paste0(heat_genes$gene_id, "  ", heat_genes$description),
                             check.names = FALSE)
rownames(z) <- rownames(annotation_row)
deg_colours <- c("up in heat" = "#D55E00", "down in heat" = "#0072B2", "not DEG" = "grey90")
gaps <- which(head(paste(samples[lib_order, "genotype"], samples[lib_order, "treatment"]), -1) !=
              tail(paste(samples[lib_order, "genotype"], samples[lib_order, "treatment"]), -1))
pheatmap(z, cluster_cols = FALSE, cluster_rows = TRUE, gaps_col = gaps,
         annotation_col = annotation, annotation_row = annotation_row,
         annotation_colors = c(annotation_colours,
                               list(`DEG genotype 21` = deg_colours, `DEG genotype 22` = deg_colours)),
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-3, 3, length.out = 101),
         main = paste0("Response to heat (GO:0009408): ", nrow(z), " Part C genes"),
         fontsize_row = 6, fontsize_col = 8,
         filename = file.path(out, "heat_response_heatmap.png"), width = 11, height = 12)

# 6. The 5 HSFs studied by Pratx et al. (2025): normalised counts per plant
# (DESeq2 size factors over all 15 libraries).
hsf <- data.frame(gene_id = c("HORVU.MOREX.r3.7HG0717290", "HORVU.MOREX.r3.2HG0134940",
                              "HORVU.MOREX.r3.5HG0519000", "HORVU.MOREX.r3.5HG0509970",
                              "HORVU.MOREX.r3.4HG0402170"),
                  name = c("HvHSFA2", "HvHSFA3", "HSFA6/7 (5HG0519000)",
                           "HSFA6/7 (5HG0509970)", "putative HSFA6 (4HG0402170)"))
dds_all <- estimateSizeFactors(dds_all)
norm <- counts(dds_all, normalized = TRUE)
hsf_long <- data.frame(gene = factor(rep(hsf$name, each = ncol(norm)), levels = hsf$name),
                       library = rep(colnames(norm), times = nrow(hsf)),
                       count = as.vector(t(norm[hsf$gene_id, ])))
hsf_long$genotype <- paste("Genotype", samples[hsf_long$library, "genotype"])
hsf_long$treatment <- samples[hsf_long$library, "treatment"]
write.table(hsf_long, file.path(out, "hsf_normalised_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p <- ggplot(hsf_long, aes(x = genotype, y = count + 1, colour = treatment)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1), size = 2) +
  stat_summary(aes(group = treatment), fun = mean, geom = "crossbar", width = 0.4,
               position = position_dodge(width = 0.6), linewidth = 0.3) +   # bar = mean
  scale_colour_manual(values = c(control = "#0072B2", heat = "#D55E00")) +
  scale_y_log10() +
  facet_wrap(~ gene, nrow = 1) +
  labs(x = NULL, y = "Normalised count + 1 (log scale)", colour = NULL,
       title = "HSFs from Pratx et al. (2025)") +
  theme_bw(base_size = 9)
ggsave(file.path(out, "hsf_expression.png"), p, width = 12, height = 4, dpi = 200)

write.table(data.frame(library = colnames(dds_all), size_factor = sizeFactors(dds_all)),
            file.path(out, "size_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
