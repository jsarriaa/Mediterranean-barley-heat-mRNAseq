# Step 4: sum kallisto transcript counts into one gene-level count matrix.
# Run from the repository root:  Rscript scripts/gene_counts.R
# Optional arguments: samplesheet and output name (Part E):
#   Rscript scripts/gene_counts.R data/input/samplesheet_E.tsv kallisto.E
# A third argument sets the kallisto folder (Part D):
#   Rscript scripts/gene_counts.R data/input/samplesheet_D.tsv kallisto.pratx2025 results/D_validation/upstream/kallisto
# Environment: barley-tximport

library(tximport)

args <- commandArgs(trailingOnly = TRUE)
samplesheet <- if (length(args) >= 1) args[1] else "data/input/samplesheet.tsv"
name <- if (length(args) >= 2) args[2] else "kallisto.merged"
kallisto_dir <- if (length(args) >= 3) args[3] else "results/upstream/kallisto"

# Libraries to include, from the samplesheet. Sorted in C-locale order
# (T10M1, T10M3, ..., T1M1, ...) so the column order does not depend on the
# computer's language settings.
samples <- read.delim(samplesheet)
ids <- sort(samples$library_id, method = "radix")

# One kallisto abundance file per library
files <- file.path(kallisto_dir, ids, "abundance.h5")
names(files) <- ids
stopifnot(all(file.exists(files)))

# Transcript -> gene table (column 1: transcript ID, column 2: gene ID)
tx2gene <- read.delim("references/MorexV3_tx2gene.tsv", header = FALSE)

# Sum the estimated counts of all transcripts of each gene.
# Counts are not scaled by transcript length or library size.
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)

# Write the matrix: gene_id, gene_name, one column per library
counts <- data.frame(gene_id = rownames(txi$counts),
                     gene_name = rownames(txi$counts),
                     txi$counts,
                     check.names = FALSE)
write.table(counts, file.path("data/processed", paste0(name, ".gene_counts.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Gene-level TPM (transcripts per million): sum of the kallisto TPM of all
# transcripts of each gene. Used for absolute expression levels (Part D).
tpm <- data.frame(gene_id = rownames(txi$abundance),
                  gene_name = rownames(txi$abundance),
                  txi$abundance,
                  check.names = FALSE)
write.table(tpm, file.path("data/processed", paste0(name, ".gene_tpm.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Genes:", nrow(counts), " Libraries:", length(ids), "\n")
