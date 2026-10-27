# Top-12 population TPM and Morex leaf baseline

## Goal

Compare expression of the 12 selected heat-response genes among the five
populations and against basal Morex V3 leaf expression from the Morex Gene
Atlas. The external Atlas values are a contextual reference, not another
population or a statistical comparison group.

## Completed steps

1. [x] **Confirm the genes.** Use
   `results/01_deseq2/heat_response_top12_gene_ids.csv` and preserve the
   Morex V3 IDs as the analysis keys.

2. [x] **Calculate population TPM.** Sum transcript TPM from each Kallisto
   sample using `references/MorexV3_tx2gene.tsv`. Outputs:
   `results/04_tpm_leaf_baseline/top12_gene_sample_tpm.csv` and
   `results/04_tpm_leaf_baseline/top12_gene_population_tpm_summary.csv`.

3. [x] **Retrieve Atlas inputs and mappings.** Script 5 obtains the Morex Gene
   Atlas TPM archive and experiment-ID metadata, then queries the Morex Gene
   Atlas relator for the top-12 Morex V3 IDs. It uses `HvMxRTD` mappings for
   the Atlas expression join and retains the returned relationships in
   `top12_morex_gene_atlas_relator_mapping.csv`.

4. [x] **Select comparison samples.** Use experiment IDs 1–3 as the
   one-week untreated first-leaf baseline (PRJEB29972). Include two-month
   senescing-leaf experiment IDs 238–243 (PRJEB14349) and the four
   uninfested 5th-leaf controls from 30-day-old plants collected at 2 hours
   (IDs 41–44, PRJNA326683). Keep other control samples in the sample table
   without pooling them into the comparison groups.

5. [x] **Match and aggregate expression.** Sum transcript TPM by matched
   `HvMxRTD` gene ID. Prefer exact intron-chain mappings; use same-strand
   non-exact mappings only when no exact mapping exists. Exclude opposite-
   strand matches. Retain all mapping relationships and flag any fallback,
   missing, or multiple match. The six senescing-leaf experiments are two
   technical batches for each of three biological replicates; average the
   batch TPMs within biological replicate before calculating the group mean
   and standard deviation. Keep the four 30-day control samples as four
   biological replicates.

6. [x] **Compare and plot.** Create the population-versus-baseline table and
   a combined population/reference summary, plus one linear-scale plot per
   gene. Keep genotype groups and replicate labels; plot all three Morex
   comparison groups separately.

7. [x] **Save and validate outputs.** Outputs are under
   `results/04_tpm_leaf_baseline/morex_gene_atlas/`: relator audit,
   300 gene-experiment rows across 25 retained leaf experiments, a 120-row
   biological-replicate table, a 36-row summary across three Morex groups,
   the 60-row primary-baseline comparison, a 96-row combined population and
   reference summary, and 12 plots. Eleven genes have an expression mapping;
   the first-ranked gene has no relator match and remains missing. Three
   mapped genes use non-exact fallbacks; eight use exact matches.

## Interpretation

The two-month Atlas leaf samples are labeled senescing leaf, so they are
age-matched context but not a fully matched physiological baseline. The
30-day uninfested 5th-leaf samples are controls from a spider-mite study,
collected at 2 hours; they provide a separate developmental and experimental
context, not a heat-stress control. Do not treat either Atlas comparison as a
sixth population or perform a matched significance test against the
population samples.

## Re-run

Activate `barley-rnaseq-upstream` and run the population TPM script first,
then the Atlas baseline script:

```bash
python3 scripts/04_population_tpm_heat_response.py
python3 scripts/05_morex_gene_atlas_baseline.py
```

The Atlas script needs internet access and the `unzip` command. It downloads
the public Morex Gene Atlas inputs if they are missing. The large TPM archive
is ignored by Git.
