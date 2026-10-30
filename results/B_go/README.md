# Part B — Functions of the genes that differ between populations

Subfolders: [`clusters/`](clusters/README.md), [`contrasts/`](contrasts/README.md).

## Index of tables

| Table | Rows | Columns |
|---|---:|---:|
| [`go_coverage.tsv`](#go_coveragetsv) | 6 | 2 |
| [`summary_significant_terms.tsv`](#summary_significant_termstsv) | 22 | 5 |

## Tables

First 5 rows of each table. Values are rounded and long text is cut for display; at most 8 columns are shown. The files hold the full values.

### go_coverage.tsv

GO annotation coverage of the HC genes.

[`go_coverage.tsv`](go_coverage.tsv): 6 rows, 2 columns.

| metric | value |
| --- | --- |
| HC genes | 35,827 |
| HC genes with any GO | 30,196 |
| HC genes without GO | 5631 |
| obsolete GO assignments replaced | 2348 |
| obsolete GO assignments dropped (no replacement) | 4080 |

### summary_significant_terms.tsv

Genes and significant GO terms per namespace, for each cluster and pairwise group.

[`summary_significant_terms.tsv`](summary_significant_terms.tsv): 22 rows, 5 columns.

| list | genes | significant_MF | significant_BP | significant_CC |
| --- | --- | --- | --- | --- |
| cluster_1 | 6058 | 11 | 51 | 16 |
| cluster_2 | 4439 | 9 | 71 | 21 |
| k10_vs_k1_up_in_k10 | 2304 | 11 | 31 | 5 |
| k10_vs_k1_up_in_k1 | 1183 | 4 | 4 | 3 |
| k10_vs_k5_up_in_k10 | 2858 | 39 | 65 | 3 |
