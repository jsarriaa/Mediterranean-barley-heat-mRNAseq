# TODO: Compare heat-response gene expression

## Goal

Compare expression of the 12 selected heat-response genes among the five
populations, then compare their population-level expression with basal leaf
expression for Morex V3 from BARLEX.

## Steps

1. [x] **Confirm the 12 genes.** Used
   `results/01_deseq2/heat_response_top12_gene_ids.csv`. Keep the Morex V3 gene
   IDs as the keys throughout the analysis.

2. [x] **Collect the population TPM values.** Used the per-sample
   `results/upstream/kallisto/<sample>/abundance.tsv` files. Join transcript
   IDs to gene IDs with `references/MorexV3_tx2gene.tsv`, sum transcript TPMs
   within each gene, and retain the 12 selected genes and samples listed in
   `data/metadata/samples.csv`. Results are in
   `results/04_tpm_leaf_baseline/top12_gene_sample_tpm.csv`.

3. [x] **Summarize by population.** Calculated the mean, standard
   deviation, median, and sample count for each gene and population. The
   summary is in
   `results/04_tpm_leaf_baseline/top12_gene_population_tpm_summary.csv`.
   Sample-level values are retained for reproducibility.

4. **Get the Morex V3 leaf baseline from BARLEX.** Find the leaf-expression
   field, tissue label, units, and source information in the BARLEX resource.
   Confirm that the values refer to Morex V3 and identify whether they are
   TPM, another normalized measure, or a transformed value. Record missing
   values and any ambiguity rather than filling them in.

5. **Match genes exactly.** Join the population summaries to the BARLEX leaf
   values using Morex V3 gene IDs. Check for duplicate IDs, unmatched genes,
   and multiple leaf measurements. Define and document how repeated leaf
   measurements are summarized, if applicable.

6. **Compare and plot.** Make one plot per gene or a faceted plot showing
   population TPM values and the BARLEX leaf baseline. Label the BARLEX value
   as a reference baseline, not as another population or replicate. Use a
   log scale only for display if needed; state the pseudocount and retain
   untransformed values in tables.

7. **Interpret carefully.** The population values are from separate RNA-seq
   samples, while BARLEX is an external Morex leaf reference. Describe the
   comparison as contextual, not as a matched statistical test. Do not test
   population samples against a single reference value as though it were a
   replicated group.

8. **Save auditable outputs.** Store the gene-level sample TPM table,
   population summary, BARLEX leaf values with source fields, matched
   comparison table, plots, and a short note describing the data fields and
   any exclusions under a new `results/04_tpm_leaf_baseline/` directory.

9. **Check the result.** Confirm that the analysis contains the 12 intended
   gene IDs, the expected samples and populations, documented baseline units,
   and explicit records for missing or duplicate matches.
