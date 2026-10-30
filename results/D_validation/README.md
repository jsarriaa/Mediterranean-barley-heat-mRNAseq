# Part D — Treatment validation

## Index of figures

| Figure | Content |
|---|---|
| [`hsf_tpm.png`](#hsf_tpmpng) | TPM + 1 (log scale) of the 5 HSFs. |

## Index of tables

| Table | Rows | Columns |
|---|---:|---:|
| [`hsf_tpm.tsv`](#hsf_tpmtsv) | 330 | 7 |
| [`hsf_tpm_means.tsv`](#hsf_tpm_meanstsv) | 5 | 8 |

## Figures

### hsf_tpm.png

![hsf_tpm.png](hsf_tpm.png)

TPM + 1 (log scale) of the 5 HSFs. One row per HSF, with its own y axis. Columns: reference experiment (`Pratx 2025`) and our 5 populations. Each dot is one library; the bar is the geometric mean of TPM + 1.

## Tables

First 5 rows of each table. Values are rounded and long text is cut for display; at most 8 columns are shown. The files hold the full values.

### hsf_tpm.tsv

TPM of the 5 HSFs per library.

[`hsf_tpm.tsv`](hsf_tpm.tsv): 330 rows, 7 columns.

| gene_id | gene | panel | x | sample | library | TPM |
| --- | --- | --- | --- | --- | --- | --- |
| HORVU.MOREX.r3.7HG0717290 | HvHSFA2 | Pratx 2025 | Control | GP_N_R3 | SRR35220183 | 0.02449 |
| HORVU.MOREX.r3.7HG0717290 | HvHSFA2 | Pratx 2025 | Control | GP_N_R2 | SRR35220184 | 0.1234 |
| HORVU.MOREX.r3.7HG0717290 | HvHSFA2 | Pratx 2025 | Control | GP_N_R1 | SRR35220185 | 0 |
| HORVU.MOREX.r3.7HG0717290 | HvHSFA2 | Pratx 2025 | Heat shock | GP_T_R3 | SRR35220177 | 755.3 |
| HORVU.MOREX.r3.7HG0717290 | HvHSFA2 | Pratx 2025 | Heat shock | GP_T_R2 | SRR35220178 | 592.8 |

### hsf_tpm_means.tsv

Arithmetic mean TPM per group: reference control, reference heat shock and each population.

[`hsf_tpm_means.tsv`](hsf_tpm_means.tsv): 5 rows, 8 columns.

| gene | Control | Heat shock | k1 | k5 | k10 | k12 | k14 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| HvHSFA2 | 0.04929 | 721.3 | 6.559 | 8.137 | 111 | 7.894 | 18.55 |
| HvHSFA3 | 0.5428 | 2.786 | 1.077 | 0.463 | 1.462 | 0.5915 | 0.2013 |
| HSFA6/7 (5HG0519000) | 7.987 | 37.22 | 15.16 | 17.21 | 39.23 | 21.07 | 29.09 |
| HSFA6/7 (5HG0509970) | 1.271 | 354.2 | 1.586 | 4.468 | 32.14 | 16.96 | 46.21 |
| putative HSFA6 (4HG0402170) | 14.09 | 175.2 | 218.3 | 209.5 | 193.7 | 254.3 | 568 |
