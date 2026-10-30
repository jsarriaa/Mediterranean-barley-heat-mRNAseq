# Part C: gene cherry-picking. Heatmap of the omnibus DEGs annotated to
# GO:0009408 "response to heat" (directly or through a more specific term).
# Run from the repository root:  Rscript scripts/C_cherrypicking.R

library(pheatmap)

out <- "results/C_cherrypicking"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
go_term <- "GO:0009408"

# 1. Genes: omnibus DEGs annotated to the term or to any term below it.
lines <- readLines("references/annotation/go-basic.obo")
block <- cumsum(lines == "[Term]" | lines == "[Typedef]")
is_id <- startsWith(lines, "id: ")
ids <- data.frame(block = block[is_id], go_id = sub("^id: ", "", lines[is_id]))
i <- which(startsWith(lines, "is_a: ") | startsWith(lines, "relationship: part_of "))
edges <- data.frame(child = ids$go_id[match(block[i], ids$block)],
                    parent = substring(sub("^(is_a|relationship: part_of): ", "", lines[i]), 1, 10))
below <- go_term
repeat {
  new <- unique(c(below, edges$child[edges$parent %in% below]))
  if (length(new) == length(below)) break
  below <- new
}
gene2go <- read.delim("references/annotation/gene2go_HC.tsv", quote = "")
lrt <- read.delim("results/A_deg/omnibus_LRT_all_genes.tsv")
degs <- lrt$gene_id[!is.na(lrt$padj) & lrt$padj < 0.05]
genes <- intersect(unique(gene2go$gene_id[gene2go$go_id %in% below]), degs)

# Gene descriptions from the Morex V3 annotation.
gtf <- readLines("references/MorexV3.gtf")
gtf <- gtf[grepl("\ttranscript\t", gtf)]
gtf_gene <- sub('.*gene_id "([^"]*)".*', "\\1", gtf)
gtf_desc <- ifelse(grepl('description "', gtf), sub('.*description "([^"]*)".*', "\\1", gtf), "")
description <- gtf_desc[match(genes, gtf_gene)]

# 2. Expression: one value per genotype. The 3 plants' counts are summed and
# variance-stabilised with the design (results/A_deg/vst_genotypes.tsv, from
# Part A), then scaled per gene (z-score across the 20 genotypes).
vst <- read.delim("results/A_deg/vst_genotypes.tsv", check.names = FALSE)
rownames(vst) <- vst$gene_id
genotypes <- paste0("T", 1:20)
z <- t(scale(t(as.matrix(vst[genes, genotypes]))))

# 3. Columns: genotypes T1-T20 in order; gaps between populations.
library_table <- read.delim("results/A_deg/library_table.tsv")
population <- library_table$population[match(genotypes, library_table$genotype)]
population <- factor(population, levels = c("k1", "k5", "k10", "k12", "k14"))
gaps <- which(head(as.character(population), -1) != tail(as.character(population), -1))
annotation_col <- data.frame(population = population, row.names = genotypes)
rownames(z) <- paste0(genes, "  ", description)

# 4. Heatmap: rows clustered (Euclidean distance, complete linkage, pheatmap default).
pheatmap(z, cluster_cols = FALSE, cluster_rows = TRUE, gaps_col = gaps,
         annotation_col = annotation_col,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-3, 3, length.out = 101),
         main = paste0("Response to heat (", go_term, "): ", length(genes), " DEGs"),
         fontsize_row = 6, fontsize_col = 8, angle_col = 0,
         filename = file.path(out, "heat_response_heatmap.png"), width = 10, height = 12)

write.table(data.frame(gene_id = genes, description = description,
                       omnibus_padj = lrt$padj[match(genes, lrt$gene_id)]),
            file.path(out, "heat_response_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 5. Pattern plot: z-scores of the population means and Part B cluster of
# each gene (results/B_go/clusters/cluster_assignment.tsv). One thin line per
# gene, coloured by cluster; thick line = cluster mean.
library(ggplot2)
populations <- c("k1", "k5", "k10", "k12", "k14")
assignment <- read.delim("results/B_go/clusters/cluster_assignment.tsv", check.names = FALSE)
assignment <- assignment[assignment$gene_id %in% genes, ]
long <- data.frame(gene_id = rep(assignment$gene_id, times = 5),
                   population = factor(rep(populations, each = nrow(assignment)), levels = populations),
                   z = unlist(assignment[, populations]),
                   cluster = rep(assignment$cluster, times = 5))
long$cluster <- paste0("Cluster ", long$cluster, " (", table(assignment$cluster)[as.character(long$cluster)], " genes)")
cluster_means <- aggregate(z ~ cluster + population, long, mean)
p <- ggplot(long, aes(x = population, y = z, group = gene_id, colour = cluster)) +
  geom_line(alpha = 0.35) +
  geom_line(data = cluster_means, aes(group = cluster), linewidth = 1.5) +
  scale_colour_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  facet_wrap(~ cluster) +   # one panel per Part B cluster
  labs(x = "Population", y = "Scaled expression (z-score of population means)",
       title = paste0("Response to heat (", go_term, "): ", nrow(assignment), " omnibus DEGs by Part B cluster")) +
  theme_bw(base_size = 10)
ggsave(file.path(out, "heat_response_patterns.png"), p, width = 10, height = 5, dpi = 200)

# 6. UpSet plot of the pairwise DEG sets (Part A, padj < 0.05), restricted to
# these genes. All non-empty intersections are shown. Layout as in Part A.
library(gridExtra)
pairwise <- read.delim("results/A_deg/pairwise_DEGs_long.tsv")
contrasts <- unique(pairwise$contrast)          # Part A order, all 10 contrasts
pairwise <- pairwise[pairwise$gene_id %in% genes, ]
membership <- table(pairwise$gene_id, factor(pairwise$contrast, levels = contrasts)) > 0
pattern <- apply(membership, 1, function(x) paste(contrasts[x], collapse = " & "))
inter <- as.data.frame(sort(table(pattern), decreasing = TRUE), stringsAsFactors = FALSE)
names(inter) <- c("pattern", "genes")
inter$id <- seq_len(nrow(inter))
write.table(inter[, c("id", "pattern", "genes")], file.path(out, "heat_response_upset_intersections.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

dots <- expand.grid(id = inter$id, contrast = contrasts)
dots$member <- mapply(function(id, set) set %in% strsplit(inter$pattern[id], " & ")[[1]],
                      dots$id, as.character(dots$contrast))
sizes <- as.data.frame(table(contrast = factor(pairwise$contrast, levels = contrasts),
                             direction = factor(ifelse(sub("up in ", "", pairwise$direction) ==
                                                       sub("_vs_.*", "", pairwise$contrast),
                                                       "up in first population", "up in second population"),
                                                levels = c("up in first population", "up in second population"))))

x_scale <- scale_x_continuous(limits = c(0.5, nrow(inter) + 0.5), expand = c(0, 0))
y_scale <- scale_y_discrete(limits = rev(contrasts))
p_top <- ggplot(inter, aes(x = id, y = genes)) +
  geom_col(fill = "grey30") +
  geom_text(aes(label = genes), vjust = -0.3, size = 2.5) +
  x_scale + labs(x = NULL, y = "Genes in intersection",
                 title = paste0("Response to heat (", go_term, "): pairwise DEGs among ", length(genes), " genes")) +
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
  scale_fill_manual(values = c("up in first population" = "#D55E00", "up in second population" = "#0072B2"),
                    drop = FALSE) +
  y_scale + labs(x = "DEGs per contrast", y = NULL, fill = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.text = element_text(size = 6), legend.key.size = unit(0.3, "cm"))

g_top <- ggplotGrob(p_top); g_mid <- ggplotGrob(p_mid); g_right <- ggplotGrob(p_right)
widths <- grid::unit.pmax(g_top$widths, g_mid$widths)
g_top$widths <- widths; g_mid$widths <- widths
heights <- grid::unit.pmax(g_mid$heights, g_right$heights)
g_mid$heights <- heights; g_right$heights <- heights
upset <- arrangeGrob(g_top, grid::nullGrob(), g_mid, g_right, ncol = 2,
                     widths = c(3, 1), heights = c(1.2, 1))
ggsave(file.path(out, "heat_response_upset.png"), upset, width = 12, height = 7, dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
