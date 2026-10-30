# Part A — Genes differentially expressed between populations

## Index of figures

| Figure | Content |
|---|---|
| [`genotype_distances.png`](#genotype_distancespng) | Same distance with the 3 plants of each genotype summed. |
| [`library_distances.png`](#library_distancespng) | Euclidean distance between the 60 libraries (variance-stabilised, all HC genes). |
| [`pairwise_upset.png`](#pairwise_upsetpng) | Overlap of the 10 pairwise DEG sets. |
| [`pca_genotypes.png`](#pca_genotypespng) | PCA of all HC genes after variance-stabilising transformation (blind to the design). |
| [`pca_libraries.png`](#pca_librariespng) | Same PCA with one dot per plant (60 libraries). |

## Index of tables

| Table | Rows | Columns |
|---|---:|---:|
| [`library_table.tsv`](#library_tabletsv) | 60 | 3 |
| [`omnibus_LRT_DEGs.tsv`](#omnibus_lrt_degstsv) | 10,497 | 7 |
| [`omnibus_LRT_all_genes.tsv`](#omnibus_lrt_all_genestsv) | 35,827 | 7 |
| [`omnibus_summary.tsv`](#omnibus_summarytsv) | 3 | 2 |
| [`pairwise_DEGs_long.tsv`](#pairwise_degs_longtsv) | 34,877 | 3 |
| [`pairwise_k10_vs_k1.tsv`](#pairwise_k10_vs_k1tsv) | 35,827 | 7 |
| [`pairwise_k10_vs_k5.tsv`](#pairwise_k10_vs_k5tsv) | 35,827 | 7 |
| [`pairwise_k12_vs_k1.tsv`](#pairwise_k12_vs_k1tsv) | 35,827 | 7 |
| [`pairwise_k12_vs_k10.tsv`](#pairwise_k12_vs_k10tsv) | 35,827 | 7 |
| [`pairwise_k12_vs_k5.tsv`](#pairwise_k12_vs_k5tsv) | 35,827 | 7 |
| [`pairwise_k14_vs_k1.tsv`](#pairwise_k14_vs_k1tsv) | 35,827 | 7 |
| [`pairwise_k14_vs_k10.tsv`](#pairwise_k14_vs_k10tsv) | 35,827 | 7 |
| [`pairwise_k14_vs_k12.tsv`](#pairwise_k14_vs_k12tsv) | 35,827 | 7 |
| [`pairwise_k14_vs_k5.tsv`](#pairwise_k14_vs_k5tsv) | 35,827 | 7 |
| [`pairwise_k5_vs_k1.tsv`](#pairwise_k5_vs_k1tsv) | 35,827 | 7 |
| [`pairwise_summary.tsv`](#pairwise_summarytsv) | 10 | 4 |
| [`size_factors_genotypes.tsv`](#size_factors_genotypestsv) | 20 | 2 |
| [`upset_top30_intersections.tsv`](#upset_top30_intersectionstsv) | 30 | 3 |
| [`vst_genotypes.tsv`](#vst_genotypestsv) | 35,827 | 21 |

## Figures

### genotype_distances.png

![genotype_distances.png](genotype_distances.png)

Same distance with the 3 plants of each genotype summed. Genotypes are clustered (complete linkage); the trees show how they group and the colour strips show the population.

### library_distances.png

![library_distances.png](library_distances.png)

Euclidean distance between the 60 libraries (variance-stabilised, all HC genes). Blue: similar; red: different. Ordered by population, genotype and plant, without clustering. White gaps separate genotypes; the side strip shows the population.

### pairwise_upset.png

![pairwise_upset.png](pairwise_upset.png)

Overlap of the 10 pairwise DEG sets. Each column of the dot matrix is a group of genes that are DEGs in exactly the marked contrasts. Top bars: genes per group (30 largest). Right bars: DEGs per contrast; orange is up in the first population, blue in the second.

### pca_genotypes.png

![pca_genotypes.png](pca_genotypes.png)

PCA of all HC genes after variance-stabilising transformation (blind to the design). Each dot is one genotype (3 plants summed), labelled T1–T20 and coloured by population. Axis labels give the variance explained.

### pca_libraries.png

![pca_libraries.png](pca_libraries.png)

Same PCA with one dot per plant (60 libraries).

## Tables

First 5 rows of each table. Values are rounded and long text is cut for display; at most 8 columns are shown. The files hold the full values.

### library_table.tsv

Library, genotype and population.

[`library_table.tsv`](library_table.tsv): 60 rows, 3 columns.

| library | genotype | population |
| --- | --- | --- |
| T10M1 | T10 | k10 |
| T10M3 | T10 | k10 |
| T10M4 | T10 | k10 |
| T11M1 | T11 | k10 |
| T11M3 | T11 | k10 |

### omnibus_LRT_DEGs.tsv

Genes with adjusted p < 0.05 in the omnibus test.

[`omnibus_LRT_DEGs.tsv`](omnibus_LRT_DEGs.tsv): 10,497 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.1HG0021100 | 624.8 | 11.1 | 1.042 | 823.1 | 7.494e-177 | 1.672e-172 |
| HORVU.MOREX.r3.7HG0717880 | 132 | 2.562 | 1.205 | 385.1 | 4.653e-82 | 5.189e-78 |
| HORVU.MOREX.r3.1HG0021110 | 55.49 | 8.633 | 1.086 | 282.6 | 6.172e-60 | 4.589e-56 |
| HORVU.MOREX.r3.3HG0247960 | 31.22 | -8.528 | 1.075 | 258.8 | 8.344e-55 | 4.652e-51 |
| HORVU.MOREX.r3.2HG0138060 | 80.41 | 10.87 | 1.205 | 244.3 | 1.108e-51 | 4.945e-48 |

### omnibus_LRT_all_genes.tsv

Omnibus likelihood-ratio test (`~ population` vs `~ 1`) for every HC gene. Its `log2FoldChange` column is k14 vs k1 only, not an omnibus effect.

[`omnibus_LRT_all_genes.tsv`](omnibus_LRT_all_genes.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.1HG0021100 | 624.8 | 11.1 | 1.042 | 823.1 | 7.494e-177 | 1.672e-172 |
| HORVU.MOREX.r3.7HG0717880 | 132 | 2.562 | 1.205 | 385.1 | 4.653e-82 | 5.189e-78 |
| HORVU.MOREX.r3.1HG0021110 | 55.49 | 8.633 | 1.086 | 282.6 | 6.172e-60 | 4.589e-56 |
| HORVU.MOREX.r3.3HG0247960 | 31.22 | -8.528 | 1.075 | 258.8 | 8.344e-55 | 4.652e-51 |
| HORVU.MOREX.r3.2HG0138060 | 80.41 | 10.87 | 1.205 | 244.3 | 1.108e-51 | 4.945e-48 |

### omnibus_summary.tsv

Gene counts of the omnibus test.

[`omnibus_summary.tsv`](omnibus_summary.tsv): 3 rows, 2 columns.

| metric | value |
| --- | --- |
| HC genes in matrix | 35,827 |
| genes with an adjusted p-value | 22,304 |
| omnibus DEGs (padj < 0.05) | 10,497 |

### pairwise_DEGs_long.tsv

Every pairwise DEG with its contrast and direction.

[`pairwise_DEGs_long.tsv`](pairwise_DEGs_long.tsv): 34,877 rows, 3 columns.

| gene_id | contrast | direction |
| --- | --- | --- |
| HORVU.MOREX.r3.4HG0416320 | k5_vs_k1 | up in k1 |
| HORVU.MOREX.r3.1HG0021100 | k5_vs_k1 | up in k5 |
| HORVU.MOREX.r3.3HG0254610 | k5_vs_k1 | up in k5 |
| HORVU.MOREX.r3.5HG0508220 | k5_vs_k1 | up in k1 |
| HORVU.MOREX.r3.7HG0717880 | k5_vs_k1 | up in k5 |

### pairwise_k10_vs_k1.tsv

Wald test, k10 vs k1. DESeq2 results for every gene: mean normalised count (`baseMean`), log2 fold change, standard error, test statistic, p and adjusted p. `log2FoldChange` > 0 means higher in k10.

[`pairwise_k10_vs_k1.tsv`](pairwise_k10_vs_k1.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.1HG0092550 | 86.55 | -24.03 | 2.106 | -11.41 | 3.825e-30 | 5.646e-26 |
| HORVU.MOREX.r3.2HG0106750 | 472.6 | -7.601 | 0.6676 | -11.39 | 4.946e-30 | 5.646e-26 |
| HORVU.MOREX.r3.3HG0237440 | 267.4 | -6.338 | 0.56 | -11.32 | 1.068e-29 | 8.127e-26 |
| HORVU.MOREX.r3.4HG0415680 | 178.4 | -24.81 | 2.585 | -9.596 | 8.304e-22 | 4.74e-18 |
| HORVU.MOREX.r3.2HG0130950 | 2200 | -1.134 | 0.1208 | -9.388 | 6.14e-21 | 2.804e-17 |

### pairwise_k10_vs_k5.tsv

Wald test, k10 vs k5. DESeq2 results for every gene: mean normalised count (`baseMean`), log2 fold change, standard error, test statistic, p and adjusted p. `log2FoldChange` > 0 means higher in k10.

[`pairwise_k10_vs_k5.tsv`](pairwise_k10_vs_k5.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.2HG0102640 | 27.27 | 33.15 | 2.342 | 14.15 | 1.784e-45 | 3.979e-41 |
| HORVU.MOREX.r3.4HG0404600 | 13.46 | 33.31 | 2.518 | 13.23 | 5.88e-40 | 6.557e-36 |
| HORVU.MOREX.r3.2HG0132100 | 213.4 | 32.93 | 2.536 | 12.99 | 1.463e-38 | 1.088e-34 |
| HORVU.MOREX.r3.5HG0513750 | 233.6 | -39.38 | 3.16 | -12.46 | 1.243e-35 | 6.93e-32 |
| HORVU.MOREX.r3.1HG0021100 | 624.8 | -10.36 | 0.8802 | -11.77 | 5.829e-32 | 2.6e-28 |

### pairwise_k12_vs_k1.tsv

Wald test, k12 vs k1. DESeq2 results for every gene: mean normalised count (`baseMean`), log2 fold change, standard error, test statistic, p and adjusted p. `log2FoldChange` > 0 means higher in k12.

[`pairwise_k12_vs_k1.tsv`](pairwise_k12_vs_k1.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.1HG0021100 | 624.8 | 11.06 | 1.042 | 10.62 | 2.38e-26 | 2.593e-22 |
| HORVU.MOREX.r3.2HG0216230 | 547.1 | -6.677 | 0.6281 | -10.63 | 2.134e-26 | 2.593e-22 |
| HORVU.MOREX.r3.2HG0132100 | 213.4 | 23.58 | 2.535 | 9.299 | 1.424e-20 | 1.034e-16 |
| HORVU.MOREX.r3.5HG0482120 | 298.1 | 25.72 | 2.863 | 8.986 | 2.569e-19 | 1.399e-15 |
| HORVU.MOREX.r3.1HG0021110 | 55.49 | 9.113 | 1.085 | 8.399 | 4.504e-17 | 1.962e-13 |

### pairwise_k12_vs_k10.tsv

Wald test, k12 vs k10. DESeq2 results for every gene: mean normalised count (`baseMean`), log2 fold change, standard error, test statistic, p and adjusted p. `log2FoldChange` > 0 means higher in k12.

[`pairwise_k12_vs_k10.tsv`](pairwise_k12_vs_k10.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.3HG0237440 | 267.4 | 6.972 | 0.5595 | 12.46 | 1.219e-35 | 2.719e-31 |
| HORVU.MOREX.r3.5HG0513750 | 233.6 | 38.36 | 3.161 | 12.14 | 6.762e-34 | 7.541e-30 |
| HORVU.MOREX.r3.1HG0021100 | 624.8 | 10.45 | 0.8802 | 11.87 | 1.679e-32 | 1.249e-28 |
| HORVU.MOREX.r3.1HG0092550 | 86.55 | 23.17 | 2.107 | 11 | 3.897e-28 | 2.173e-24 |
| HORVU.MOREX.r3.2HG0216230 | 547.1 | -6.122 | 0.6282 | -9.745 | 1.933e-22 | 8.621e-19 |

### pairwise_k12_vs_k5.tsv

Wald test, k12 vs k5. DESeq2 results for every gene: mean normalised count (`baseMean`), log2 fold change, standard error, test statistic, p and adjusted p. `log2FoldChange` > 0 means higher in k12.

[`pairwise_k12_vs_k5.tsv`](pairwise_k12_vs_k5.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.5HG0482120 | 298.1 | 40.19 | 2.863 | 14.04 | 8.95e-45 | 1.95e-40 |
| HORVU.MOREX.r3.2HG0102640 | 27.27 | 32.42 | 2.346 | 13.82 | 1.906e-43 | 2.076e-39 |
| HORVU.MOREX.r3.2HG0132100 | 213.4 | 34.38 | 2.536 | 13.56 | 7.032e-42 | 5.106e-38 |
| HORVU.MOREX.r3.4HG0404600 | 13.46 | 32.44 | 2.521 | 12.87 | 6.592e-38 | 3.59e-34 |
| HORVU.MOREX.r3.4HG0416320 | 44.62 | 23.06 | 1.856 | 12.43 | 1.847e-35 | 8.046e-32 |

### pairwise_k14_vs_k1.tsv

Wald test, k14 vs k1. DESeq2 results for every gene: mean normalised count (`baseMean`), log2 fold change, standard error, test statistic, p and adjusted p. `log2FoldChange` > 0 means higher in k14.

[`pairwise_k14_vs_k1.tsv`](pairwise_k14_vs_k1.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.2HG0216230 | 547.1 | -8.498 | 0.6978 | -12.18 | 4.097e-34 | 4.569e-30 |
| HORVU.MOREX.r3.5HG0497560 | 116.6 | -23.64 | 1.934 | -12.23 | 2.214e-34 | 4.569e-30 |
| HORVU.MOREX.r3.7HG0666740 | 7208 | 3.177 | 0.2736 | 11.61 | 3.565e-31 | 2.65e-27 |
| HORVU.MOREX.r3.7HG0712330 | 1089 | -2.199 | 0.1996 | -11.02 | 3.169e-28 | 1.767e-24 |
| HORVU.MOREX.r3.1HG0021100 | 624.8 | 11.1 | 1.042 | 10.66 | 1.611e-26 | 7.185e-23 |

### pairwise_k14_vs_k10.tsv

Wald test, k14 vs k10. DESeq2 results for every gene: mean normalised count (`baseMean`), log2 fold change, standard error, test statistic, p and adjusted p. `log2FoldChange` > 0 means higher in k14.

[`pairwise_k14_vs_k10.tsv`](pairwise_k14_vs_k10.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.7HG0712330 | 1089 | -2.594 | 0.1994 | -13.01 | 1.082e-38 | 2.302e-34 |
| HORVU.MOREX.r3.5HG0513750 | 233.6 | 40.51 | 3.16 | 12.82 | 1.301e-37 | 1.384e-33 |
| HORVU.MOREX.r3.5HG0497560 | 116.6 | -23.61 | 1.934 | -12.21 | 2.762e-34 | 1.958e-30 |
| HORVU.MOREX.r3.1HG0021100 | 624.8 | 10.49 | 0.8802 | 11.91 | 1.002e-32 | 5.33e-29 |
| HORVU.MOREX.r3.2HG0216230 | 547.1 | -7.943 | 0.698 | -11.38 | 5.258e-30 | 2.237e-26 |

### pairwise_k14_vs_k12.tsv

Wald test, k14 vs k12. DESeq2 results for every gene: mean normalised count (`baseMean`), log2 fold change, standard error, test statistic, p and adjusted p. `log2FoldChange` > 0 means higher in k14.

[`pairwise_k14_vs_k12.tsv`](pairwise_k14_vs_k12.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.5HG0497560 | 116.6 | -24.71 | 1.933 | -12.78 | 2.096e-37 | 4.46e-33 |
| HORVU.MOREX.r3.7HG0736260 | 34.85 | -22.43 | 2.168 | -10.34 | 4.546e-25 | 4.836e-21 |
| HORVU.MOREX.r3.1HG0077240 | 324.7 | -8.194 | 0.8421 | -9.73 | 2.244e-22 | 1.592e-18 |
| HORVU.MOREX.r3.4HG0408930 | 1.627e+04 | -4.373 | 0.4557 | -9.595 | 8.371e-22 | 4.453e-18 |
| HORVU.MOREX.r3.7HG0712330 | 1089 | -1.911 | 0.1999 | -9.561 | 1.166e-21 | 4.96e-18 |

### pairwise_k14_vs_k5.tsv

Wald test, k14 vs k5. DESeq2 results for every gene: mean normalised count (`baseMean`), log2 fold change, standard error, test statistic, p and adjusted p. `log2FoldChange` > 0 means higher in k14.

[`pairwise_k14_vs_k5.tsv`](pairwise_k14_vs_k5.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.2HG0102640 | 27.27 | 35.5 | 2.338 | 15.19 | 4.402e-52 | 9.588e-48 |
| HORVU.MOREX.r3.2HG0132100 | 213.4 | 34.91 | 2.535 | 13.77 | 3.827e-43 | 4.168e-39 |
| HORVU.MOREX.r3.4HG0404600 | 13.46 | 33.49 | 2.517 | 13.31 | 2.156e-40 | 1.566e-36 |
| HORVU.MOREX.r3.4HG0416320 | 44.62 | 24.01 | 1.854 | 12.95 | 2.381e-38 | 1.297e-34 |
| HORVU.MOREX.r3.5HG0497560 | 116.6 | -23.75 | 1.934 | -12.28 | 1.111e-34 | 4.842e-31 |

### pairwise_k5_vs_k1.tsv

Wald test, k5 vs k1. DESeq2 results for every gene: mean normalised count (`baseMean`), log2 fold change, standard error, test statistic, p and adjusted p. `log2FoldChange` > 0 means higher in k5.

[`pairwise_k5_vs_k1.tsv`](pairwise_k5_vs_k1.tsv): 35,827 rows, 7 columns.

| gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | padj |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.4HG0416320 | 44.62 | -20.39 | 1.867 | -10.92 | 8.931e-28 | 1.945e-23 |
| HORVU.MOREX.r3.1HG0021100 | 624.8 | 10.97 | 1.042 | 10.53 | 6.106e-26 | 6.65e-22 |
| HORVU.MOREX.r3.3HG0254610 | 159.9 | 23.76 | 2.427 | 9.788 | 1.271e-22 | 9.227e-19 |
| HORVU.MOREX.r3.5HG0508220 | 27.25 | -20.56 | 2.156 | -9.536 | 1.485e-21 | 8.086e-18 |
| HORVU.MOREX.r3.7HG0717880 | 132 | 10.36 | 1.131 | 9.163 | 5.032e-20 | 2.192e-16 |

### pairwise_summary.tsv

DEGs per pairwise contrast and direction.

[`pairwise_summary.tsv`](pairwise_summary.tsv): 10 rows, 4 columns.

| contrast | DEGs | up_in_first | down_in_first |
| --- | --- | --- | --- |
| k5_vs_k1 | 730 | 323 | 407 |
| k10_vs_k1 | 6114 | 3742 | 2372 |
| k12_vs_k1 | 950 | 579 | 371 |
| k14_vs_k1 | 4878 | 2532 | 2346 |
| k10_vs_k5 | 7441 | 4410 | 3031 |

### size_factors_genotypes.tsv

DESeq2 size factors (sequencing depth correction) per genotype.

[`size_factors_genotypes.tsv`](size_factors_genotypes.tsv): 20 rows, 2 columns.

| genotype | size_factor |
| --- | --- |
| T1 | 1.16 |
| T10 | 1.107 |
| T11 | 0.9567 |
| T12 | 0.9318 |
| T13 | 1.002 |

### upset_top30_intersections.tsv

The 30 intersections shown in the UpSet plot.

[`upset_top30_intersections.tsv`](upset_top30_intersections.tsv): 30 rows, 3 columns.

| id | pattern | genes |
| --- | --- | --- |
| 1 | k10_vs_k5 | 1185 |
| 2 | k10_vs_k1 & k10_vs_k5 | 862 |
| 3 | k10_vs_k1 & k10_vs_k5 & k12_vs_k10 & k14_vs_k10 | 754 |
| 4 | k10_vs_k1 | 624 |
| 5 | k14_vs_k1 | 602 |

### vst_genotypes.tsv

Variance-stabilised expression per gene and genotype (design-aware, not blind).

[`vst_genotypes.tsv`](vst_genotypes.tsv): 35,827 rows, 21 columns.

All columns: `gene_id`, `T1`, `T10`, `T11`, `T12`, `T13`, `T14`, `T15`, `T16`, `T17`, `T18`, `T19`, `T2`, `T20`, `T3`, `T4`, `T5`, `T6`, `T7`, `T8`, `T9`.

| gene_id | T1 | T10 | T11 | T12 | T13 | T14 | T15 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.1HG0000030 | 7.462 | 8.108 | 8.108 | 9.936 | 8.369 | 9.726 | 8.101 |
| HORVU.MOREX.r3.1HG0000040 | 13.61 | 12.67 | 13.36 | 13.49 | 13.43 | 13.22 | 13.68 |
| HORVU.MOREX.r3.1HG0000050 | 10.57 | 9.567 | 9.702 | 9.426 | 9.886 | 10.24 | 9.693 |
| HORVU.MOREX.r3.1HG0000060 | 12.11 | 12.11 | 10.5 | 11.17 | 11.67 | 12.02 | 11.6 |
| HORVU.MOREX.r3.1HG0000070 | 12.36 | 12.23 | 11.77 | 12.53 | 11.8 | 11.75 | 12.78 |
