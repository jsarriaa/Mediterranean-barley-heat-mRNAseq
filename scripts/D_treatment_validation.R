# Part D: treatment validation. TPM of the 5 HSFs studied by Pratx et al.
# (2025): their control and heat shock libraries next to our 60 libraries.
# Run from the repository root:  Rscript scripts/D_treatment_validation.R
# Needs the TPM matrices from scripts/gene_counts.R (README Step 4 and Step D3).

library(ggplot2)

out <- "results/D_validation"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# 1. The 5 HSFs (Part C table).
hsf <- data.frame(gene_id = c("HORVU.MOREX.r3.7HG0717290", "HORVU.MOREX.r3.2HG0134940",
                              "HORVU.MOREX.r3.5HG0519000", "HORVU.MOREX.r3.5HG0509970",
                              "HORVU.MOREX.r3.4HG0402170"),
                  name = c("HvHSFA2", "HvHSFA3", "HSFA6/7 (5HG0519000)",
                           "HSFA6/7 (5HG0509970)", "putative HSFA6 (4HG0402170)"))

# 2. Gene-level TPM: reference and our libraries.
ref <- read.delim("data/processed/kallisto.pratx2025.gene_tpm.tsv", check.names = FALSE)
ours <- read.delim("data/processed/kallisto.merged.gene_tpm.tsv", check.names = FALSE)
stopifnot(identical(ref$gene_id, ours$gene_id))
rownames(ref) <- ref$gene_id
rownames(ours) <- ours$gene_id

# 3. Library tables. Reference: control or heat shock. Ours: genotype and
# population (Part A), with genotypes ordered T1-T20 inside each population.
ref_samples <- read.delim("data/input/samplesheet_D.tsv")
our_samples <- read.delim("results/A_deg/library_table.tsv")
genotypes <- paste0("T", 1:20)
populations <- c("k1", "k5", "k10", "k12", "k14")
pop_of <- our_samples$population[match(genotypes, our_samples$genotype)]
genotype_order <- genotypes[order(match(pop_of, populations))]

libraries <- rbind(
  data.frame(library = ref_samples$library_id, sample = ref_samples$sample,
             x = ifelse(ref_samples$condition == "control", "Control", "Heat shock"),
             panel = "Pratx 2025", group = ref_samples$condition),
  data.frame(library = our_samples$library, sample = our_samples$library,
             x = our_samples$genotype, panel = our_samples$population,
             group = our_samples$population))
libraries$x <- factor(libraries$x, levels = c("Control", "Heat shock", genotype_order))
libraries$panel <- factor(libraries$panel, levels = c("Pratx 2025", populations))

# 4. One row per gene and library.
tpm <- cbind(ref[hsf$gene_id, ref_samples$library_id], ours[hsf$gene_id, our_samples$library])
long <- data.frame(gene_id = rep(hsf$gene_id, times = ncol(tpm)),
                   gene = factor(rep(hsf$name, times = ncol(tpm)), levels = hsf$name),
                   library = rep(colnames(tpm), each = nrow(tpm)),
                   TPM = as.vector(as.matrix(tpm)))
long <- cbind(long, libraries[match(long$library, libraries$library), c("sample", "x", "panel", "group")])
write.table(long[order(long$gene, long$panel, long$x, long$library),
                 c("gene_id", "gene", "panel", "x", "sample", "library", "TPM")],
            file.path(out, "hsf_tpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Mean TPM per gene: reference control, reference heat shock, each population.
long$summary_group <- ifelse(long$panel == "Pratx 2025", as.character(long$x), as.character(long$panel))
means <- aggregate(TPM ~ gene + summary_group, long, mean)
means <- reshape(means, idvar = "gene", timevar = "summary_group", direction = "wide")
names(means) <- sub("^TPM\\.", "", names(means))
means <- means[, c("gene", "Control", "Heat shock", populations)]
write.table(means, file.path(out, "hsf_tpm_means.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(means)

# 5. Plot: one row per HSF, one dot per library. Bar = mean of each x group,
# computed on the log scale of the axis (geometric mean of TPM + 1).
colours <- c(control = "#0072B2", `heat shock` = "#D55E00", k1 = "#E69F00", k5 = "#56B4E9",
             k10 = "#009E73", k12 = "#CC79A7", k14 = "grey40")
p <- ggplot(long, aes(x = x, y = TPM + 1, colour = group)) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 1.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, linewidth = 0.25, show.legend = FALSE) +
  scale_colour_manual(values = colours, breaks = names(colours)) +
  scale_y_log10() +
  facet_grid(gene ~ panel, scales = "free", space = "free_x",
             labeller = labeller(gene = label_wrap_gen(14))) +
  labs(x = NULL, y = "TPM + 1 (log scale)", colour = NULL,
       title = "HSFs from Pratx et al. (2025): reference libraries and our 20 genotypes") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0))
ggsave(file.path(out, "hsf_tpm.png"), p, width = 12, height = 9, dpi = 200)

writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
