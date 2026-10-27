# Transcriptomic fingerprints across five Mediterranean barley populations

This repository reproduces the 1–20 population analysis of barley RNA-seq
gene-count data. It compares 60 samples from genotypes T1–T20, grouped into
five Mediterranean populations (`k1`, `k5`, `k10`, `k12`, `k14`), identifies
genes whose expression differs among populations, and investigates the
heat-stress response gene set.

The scope is deliberately limited to samples T1–T20 and those five groups.

## Pipeline scope and data availability

This protocol describes the complete mRNA-seq analysis for the T1–T20
population comparison: paired FASTQ preparation and quality control, Morex V3
transcriptome quantification, gene-count matrix generation, DESeq2 testing,
functional enrichment, heat-response candidate selection, and annotation.

The FASTQ files are intended for release through the European Nucleotide
Archive (ENA). Their public accession and download link are under embargo and
will be added here when release is authorized. The sample manifest and
per-sample QC summary are included now. The reference genome and annotation
are the public Morex V3 release from IPK Gatersleben.

This repository includes individual downstream R scripts and documents the
upstream command protocol. The large FASTQs, reference files, count matrix,
and generated full statistical outputs are excluded from Git. Run each step
separately in the order specified below.

## Repository structure

```text
transcriptomic-fingerprints-five-mediterranean-populations/
├── README.md
├── .gitignore
├── envs/
│   ├── rnaseq.yml
│   ├── deseq2.yml
│   └── conda-linux-64.lock
├── data/
│   ├── raw/                                   # ignored paired FASTQ inputs
│   ├── processed/                             # ignored fastp-cleaned reads
│   ├── input/
│   │   └── kallisto.merged.gene_counts.tsv    # ignored gene-count matrix
│   └── metadata/
│       ├── samples.csv                        # 60 T1–T20 samples
│       ├── fastq_manifest.csv                 # ENA sample-to-read mapping
│       ├── fastp_qc_summary.csv               # per-sample QC metrics
│       └── heat_stress_gene_ids.csv           # 46 GO:0009408 genes
├── references/
│   ├── MorexV3.fa                              # ignored reference genome
│   ├── MorexV3.gff3                            # ignored gene annotation
│   └── BARLEX_MorexV3_annotations.tsv          # ignored functional annotations
├── scripts/
│   ├── 01_deseq2_populations.R
│   ├── 02_go_enrichment.R
│   ├── 03_heat_stress_and_annotation.R
├── results/
│   ├── upstream/
│   │   ├── 00_fastp_qc/                        # ignored per-library reports
│   │   └── 01_kallisto/                       # ignored quantification outputs
│   ├── 01_deseq2/
│   ├── 02_go_enrichment/
│   └── 03_heat_stress/
└── docs/
    └── figures/                               # tracked report figures
```

Generated full differential-expression and enrichment tables, logs, and
intermediate data are written beneath `results/` and ignored by Git. Compact
summary tables, the selected 46/12 gene lists, and small annotation tables are
tracked. The selected report figures under `docs/figures/` are tracked and
embedded below.

## Study design and sample selection

The count matrix contains columns from a broader sequencing set. Only sample
IDs matching `T1M*` through `T20M*` are used. The committed sample sheet lists
the exact 60 retained libraries and their population assignments. There are
12 libraries per population, representing four genotype IDs per population
with three available replicates per genotype.

| Population | Genotypes | Samples |
|---|---|---:|
| k1 | T1–T4 | 12 |
| k5 | T5–T8 | 12 |
| k10 | T9–T12 | 12 |
| k12 | T13–T16 | 12 |
| k14 | T17–T20 | 12 |

The samples and group labels are explicit in
[`data/metadata/samples.csv`](data/metadata/samples.csv); they are never
reconstructed from column order or guessed during analysis. IDs outside T1–T20
are excluded by requiring an exact match between this table and the selected
count-matrix columns.

## Raw FASTQ inputs and sample mapping

The 60 paired-end Illumina libraries are listed in
[`data/metadata/fastq_manifest.csv`](data/metadata/fastq_manifest.csv), with
sample ID, population, library strandedness, and the paired read filenames.
Each mate is a 150-cycle read, and the libraries are reverse stranded. FASTQ
filenames and ENA accessions should be linked in the manifest when the archive
release is public.

These are paired FASTQs, not one file per sample. Do not infer the sample IDs
from filename order: use the committed manifest. The analysis uses 12
libraries per population, 60 libraries in total. The count matrix contains
only these T1–T20 samples.

The FASTQ data are scheduled for ENA release. The accession and download URL
remain under embargo and will be added to the manifest and protocol when
public release is authorized. The raw reads are excluded from version
control.

## FASTQ quality control and preparation

The protocol processes the paired FASTQ files with `fastp`, which writes
cleaned paired reads plus one JSON and one HTML report per library. The reports
record before/after metrics for all 60 T1–T20 libraries. A compact per-sample
table extracted from those reports is committed in
[`data/metadata/fastp_qc_summary.csv`](data/metadata/fastp_qc_summary.csv).
`fastp` version 1.3.6 processed paired reads and reported before/after
read counts, base-quality distributions, adapter trimming, and filtering
counts. The configured options were:

```text
-q 30 -u 15 -y -g -Y 10 -l 100 -b 150 -B 150
```

In practical terms, these settings trim poly-G tails, enable low-complexity
filtering, use the configured complexity threshold, trim/detect adapters using
`fastp`'s paired-read handling, discard reads shorter than 100 bases, cap each
mate at 150 bases, and apply Q30 / unqualified-base filtering. `fastp` HTML and
JSON reports—not just the final read files—were used to inspect the outcome.

Across the 60 selected libraries, the reports record 1,891,663,404 input read
pairs and 1,578,532,879 retained pairs (83.45% overall). Mean per-library Q30
base fractions changed from 94.86% to 98.65% for read 1 and from 93.58% to
98.42% for read 2. These are aggregate QC summaries, not a claim that every
library had identical quality; use the per-sample table for individual
results. The selected libraries were not removed based on these summaries by
the current analysis; T1–T20 scope and the sample sheet define the retained
set.

Run `fastp` once for each row in `data/metadata/fastq_manifest.csv`. For each
sample, set `sample`, `R1`, and `R2` to the sample ID and the two corresponding
manifest filenames. For example, assign them for `T1M1`, then run:

```bash
sample=T1M1
R1=Unknown_CU317-001T0001_1.fq.gz
R2=Unknown_CU317-001T0001_2.fq.gz
mkdir -p data/processed/fastp results/upstream/00_fastp_qc
fastp \
  -i "data/raw/${R1}" -I "data/raw/${R2}" \
  -o "data/processed/fastp/${sample}_R1.fastq.gz" \
  -O "data/processed/fastp/${sample}_R2.fastq.gz" \
  --html "results/upstream/00_fastp_qc/${sample}.fastp.html" \
  --json "results/upstream/00_fastp_qc/${sample}.fastp.json" \
  -q 30 -u 15 -y -g -Y 10 -l 100 -b 150 -B 150
```

Confirm that each report exists, paired input/output counts are consistent,
and per-sample retained-read and quality metrics are acceptable before
quantification. The committed summary is an audit of this QC stage; it is not
a replacement for reviewing individual HTML reports.

## Reference and Kallisto quantification

Use the Morex V3 barley reference genome and gene annotation released by
[IPK Gatersleben](https://doi.org/10.5447/ipk/2021/3). The release,
“Pseudomolecules and annotation of the third version of the reference genome
sequence assembly of barley cv. Morex [Morex V3],” provides genomic
pseudomolecules in FASTA format and structural gene models in GFF3 format.
Download the genome FASTA and gene annotation from this release and record the
release DOI with the analysis:

```text
references/MorexV3.fa
references/MorexV3.gff3
```

Prepare a GTF and transcript FASTA that are both derived from the same
Morex V3 GFF3/FASTA release. If the release does not provide a transcript
FASTA directly, extract transcript sequences and convert the annotation with
`gffread`. Build the Kallisto index once from that transcript FASTA. The index
and transcript-to-gene relationships must come from the same reference
release. Kallisto 0.52.0 pseudoaligns the paired cleaned reads with
reverse-stranded settings (`--rf-stranded`) and 12 threads. The recorded run
used zero bootstrap replicates. This is transcriptome pseudoalignment and
quantification, not splice-aware genomic alignment; a BAM is not used as input
to this analysis.

Example reference preparation and index commands:

```bash
gffread -T references/MorexV3.gff3 \
  -o references/MorexV3.gtf
gffread -w references/MorexV3.transcripts.fa \
  -g references/MorexV3.fa references/MorexV3.gff3
kallisto index -i references/MorexV3.idx \
  references/MorexV3.transcripts.fa
```

Pin the `gffread` version used for this conversion in the record of any
reproduction. Use the resulting GTF to map transcript estimates back to genes.

Run the quantification once per sample using the paired fastp outputs:

```bash
kallisto quant \
  --threads 12 \
  --index references/MorexV3.idx \
  --gtf references/MorexV3.gtf \
  --rf-stranded \
  -o "results/upstream/01_kallisto/${sample}" \
  "data/processed/fastp/${sample}_R1.fastq.gz" \
  "data/processed/fastp/${sample}_R2.fastq.gz"
```

The GTF option links transcript targets back to gene IDs for gene-level
summarization. The `--rf-stranded` setting matters because the sample sheet
records reverse-stranded libraries; using the wrong strand mode would alter
the estimated abundance assigned to transcripts.

The run records 83,990 Kallisto targets and a 31-base k-mer. Across the 60
selected libraries, weighted pseudoalignment was 45.23% (per-library values
ranged from 43.4% to 50.8%). Per-sample Kallisto run metadata record the processed and pseudoaligned read
counts, software version, strandedness, and command.

Kallisto writes per-sample transcript abundance estimates. These were merged
and summarized by gene ID into
`kallisto.merged.gene_counts.tsv`, with `gene_id` and `gene_name` columns and
60 sample columns. That matrix contains 81,687 gene rows. The DESeq2 analysis
rounds estimated counts to integers before modeling; it does not use TPM as
DESeq2 count input.

### Gene-count matrix and annotation resources

The matrix `data/input/kallisto.merged.gene_counts.tsv` is produced by merging
the Kallisto per-sample estimated counts and summarizing transcript estimates
to gene IDs using the matching transcript-to-gene map. Preserve gene IDs and
sample names as columns in the matrix. The committed matrix checksum
identifies the exact input used for the reported downstream results.

The independently required BARLEX table is staged at
`references/BARLEX_MorexV3_annotations.tsv`. It supplies Morex V3 gene
descriptions, Pfam, and InterPro fields to supplement online gene annotation.
Use the Morex V3 genome/GFF3 from the IPK release for read quantification, and
the BARLEX table as a separate functional annotation resource when describing
the candidate genes.

| Input | SHA-256 |
|---|---|
| `kallisto.merged.gene_counts.tsv` (T1–T20 only) | `dfd94f50cfe070009f8971d38d2ed2c8718d14dcf8a075e675f26b5cdf7e2108` |
| `BARLEX_MorexV3_annotations.tsv` | `5f25fcfcb2d69346fdfa6be80dbc0cc790de3cad8f2b58cc53bb46115d9e9275` |
| `samples.csv` | `b1004640c8be6ef79abd7aa0052b8c836e3be1d62ceec2d106aeab6a5f371f35` |
| `fastq_manifest.csv` | `8dd97c18544044b1d03cb8cfc4d4088d17393d1fd9598a4ca743f4d8835f662c` |
| `fastp_qc_summary.csv` | `cf5754c2d5176af4d763f0e3059b785662e84b4afbbf9fbade333b9867dafd9f` |
| `heat_stress_gene_ids.csv` | `42218f0c0d38fafccba2ae40e6c7a3d6ca333864e1f024252bdcec5e123f6a88` |

After generating the matrix as described below, verify the input files:

```bash
echo 'dfd94f50cfe070009f8971d38d2ed2c8718d14dcf8a075e675f26b5cdf7e2108  data/input/kallisto.merged.gene_counts.tsv' | sha256sum --check
echo '5f25fcfcb2d69346fdfa6be80dbc0cc790de3cad8f2b58cc53bb46115d9e9275  references/BARLEX_MorexV3_annotations.tsv' | sha256sum --check
```

The count matrix, reference sequences, GFF3, and BARLEX table are ignored by
Git because of their size or provenance. The FASTQ archive accession will be
added after the embargo is lifted. Do not add credentials or cloud
configuration to Git.

## Software and versions

Versions were taken from the per-sample quantification reports and the active
analysis environment:

| Program/package | Version |
|---|---:|
| fastp (paired-read filtering and QC) | 1.3.6 |
| Kallisto (transcriptome pseudoalignment and quantification) | 0.52.0 |
| R | 4.5.3 |
| DESeq2 | 1.50.2 |
| ggplot2 | 4.0.3 |
| pheatmap | 1.0.13 |
| gprofiler2 | 0.2.4 |
| gffread (transcript FASTA/GTF preparation) | 0.12.9 (protocol environment) |
| seqtk | not required by this protocol |

The R package versions are pinned by the downstream Conda environment files.
`envs/deseq2.yml` specifies direct R dependencies; the accompanying Linux
x86-64 explicit lock records resolved packages for exact recreation on that
platform. `envs/rnaseq.yml` pins fastp, gffread, and Kallisto for upstream read
processing and quantification. fastp and Kallisto versions are recorded for
the analyzed reads; gffread 0.12.9 is the version pinned in this protocol for
reference preparation. The reference files and Kallisto index are derived
from the cited Morex V3 release.

Create the upstream tools environment with:

```bash
conda env create --file envs/rnaseq.yml
conda activate barley-rnaseq-upstream
fastp --version
kallisto version
```

Kallisto indexing requires a transcript FASTA and transcript-to-gene
annotation. Use the transcript sequences and gene identifiers derived from
the same Morex V3 GFF3. The upstream environment pins `gffread` for transcript
extraction. For a different operating system, solve each YAML separately and
retain platform-specific locks.

Create and activate the environment:

```bash
conda create --name five-populations-deseq2 --file envs/conda-linux-64.lock
conda activate five-populations-deseq2
```

The workflow requires R packages `DESeq2`, `ggplot2`, `pheatmap`, and
`gprofiler2`. g:Profiler is an online service: enrichment and gene conversion
require network access, and results may change as its database is updated.
Every returned table and the query timestamp are saved under
`results/02_go_enrichment/` or `results/03_heat_stress/`.

## Step-by-step analysis protocol

Perform the upstream steps in order: stage the FASTQ pairs when the ENA release
is public, run `fastp`, prepare the matching Morex V3 transcript reference,
run Kallisto per sample, and combine estimated gene counts into the matrix.
Then run the three R scripts below, one at a time. The commands in this
protocol assume the filenames and project-relative paths shown in the
repository structure.

### Prepare the gene-count matrix

For each Kallisto output, retain the estimated counts and map each transcript
to its gene using the same Morex V3 GTF used for indexing and quantification.
Sum transcript estimates for each gene and join the 60 sample columns by the
sample IDs in `data/metadata/samples.csv`. Write
`data/input/kallisto.merged.gene_counts.tsv` with `gene_id`, `gene_name`, and
one column per sample. Keep the unrounded estimated counts in this file; the
DESeq2 script validates and rounds them immediately before model fitting.
Verify the completed matrix against the SHA-256 value in this README.

### Run the R scripts separately

Run all commands from the repository root with the `five-populations-deseq2`
environment active. Confirm the required count matrix and BARLEX annotation
are present before starting. There is intentionally no combined shell runner:
inspect each step's outputs before continuing.

#### Script 1 — differential expression and exploratory plots

```bash
Rscript scripts/01_deseq2_populations.R
```

Validates the sample sheet, population labels, count values, and exact match
between the 60 declared sample IDs and matrix columns. It rounds the estimated
counts, fits the population omnibus DESeq2 likelihood-ratio test and all ten
pairwise Wald contrasts, calculates VST values, and produces PCA and heatmaps.
Outputs are written to `results/01_deseq2/`. Check the analysis and pairwise
summary tables and confirm that the significant-gene file exists before the
GO step.

#### Script 2 — GO enrichment

```bash
Rscript scripts/02_go_enrichment.R
```

Uses the significant omnibus genes from Script 1 as the query to g:Profiler
for barley (`hvulgare`). It saves the complete enrichment table, plot,
GO:0009408 intersection, term summary, and query metadata under
`results/02_go_enrichment/`. It requires an internet connection. Check the
term summary and query metadata before running annotation; the online database
can change over time.

#### Script 3 — heat-response gene annotation

```bash
Rscript scripts/03_heat_stress_and_annotation.R
```

Uses the committed 46-gene GO:0009408 list and the top-12 variation ranking
from Script 1. It retrieves barley ID conversions and descriptions from
g:Profiler, joins BARLEX description/Pfam/InterPro by exact gene ID, and
retains all 46 candidates. It requires internet access and
`references/BARLEX_MorexV3_annotations.tsv`. Outputs and annotation coverage
summaries are written to `results/03_heat_stress/`.

## Analysis logic and why each step is used

1. **Keep the intended biological comparison explicit.** The matrix contains
   only the 60 samples in `samples.csv`: T1–T20, with 12 libraries in each of
   five populations. Other sample groups are excluded from this analysis.
   Exact sample-ID matching prevents sample order or cohort membership from
   silently changing the design.
2. **Model population-level expression differences.** Round Kallisto estimated
   counts to integers and fit a DESeq2 negative-binomial model with design
   `~ subpopulation`. The likelihood-ratio test compares that model with
   `~ 1`, testing whether expression differs anywhere among the five
   populations; it does not answer which pair differs. Benjamini–Hochberg
   adjusted p-values control the false discovery rate. A separate Wald model
   supplies all ten population-pair contrasts for follow-up.
3. **Visualize global expression.** Apply DESeq2's variance-stabilizing
   transformation to make expression patterns more comparable across the
   count range. Use those transformed values for PCA and heatmaps. The top-gene
   heatmap uses the most significant omnibus genes and centers each gene across
   samples to emphasize relative patterns, not absolute expression.
4. **Perform functional enrichment on the population-associated genes.**
   Query the omnibus significant-gene IDs against g:Profiler's barley
   (`hvulgare`) annotations. Enrichment asks whether a biological annotation
   appears more often in the significant set than expected; it does not show
   that a treatment caused the expression change.
5. **Focus interpretation on heat response after enrichment.** GO:0009408,
   “response to heat,” was significant in the barley enrichment result
   (adjusted p-value 0.02735), with 46 intersecting genes. Those 46 IDs are
   recorded in `data/metadata/heat_stress_gene_ids.csv`. This is a
   post-enrichment biological focus, not a filter applied to reads or to the
   initial DESeq2 test. It narrows interpretation to a relevant
   process while preserving the all-gene differential-expression and
   enrichment results. The 12 focused genes are ranked reproducibly by the
   standard deviation of each gene's five population-mean VST expression
   values; they are not selected by subjective inspection of the heatmap.
6. **Use two complementary annotation resources, not two orthology tests.**
   g:Profiler supplies online barley ID conversion/descriptions and is the
   source of GO enrichment. The separate BARLEX Morex V3 table supplies
   reference-matched descriptions plus Pfam and InterPro fields. Using both
   provides complementary functional descriptions and domain annotations;
   the code prefers a non-empty g:Profiler description and falls back to
   BARLEX when it is missing or `"None"`. BARLEX is not used to calculate the
   GO enrichment, and the two resources are not being treated as an
   orthology-based differential-expression filter.
7. **Retain and report all candidates.** Join annotations by exact gene ID,
   preserve all 46 genes, and record genes without a description rather than
   silently dropping them. Save the focused top-12 table, figures, complete
   statistical/enrichment outputs, summaries, session information, and logs
   under their result folders.

## Reported analysis results

The DESeq2 LRT identified 20,306 genes with adjusted p-value < 0.05 from
53,244 genes with testable p-values; 34,175 genes had a non-missing adjusted
p-value after independent filtering. The input matrix contained 81,687 gene
rows. Pairwise Wald contrasts produced the following significant-gene counts:

| Contrast | Significant genes |
|---|---:|
| k5 vs k1 | 2,830 |
| k10 vs k1 | 11,351 |
| k12 vs k1 | 4,016 |
| k14 vs k1 | 9,754 |
| k10 vs k5 | 12,963 |
| k12 vs k5 | 3,504 |
| k14 vs k5 | 9,407 |
| k12 vs k10 | 9,102 |
| k14 vs k10 | 7,976 |
| k14 vs k12 | 5,343 |

These are DESeq2 results from the population model, not claims of heat
treatment effects. The live barley g:Profiler query returned 207 significant
terms; GO:0009408 (“response to heat”) had an adjusted p-value of 0.02735 and
an intersection of 46 genes. All 46 exactly match the versioned gene set. The
service reported 7,327 mapped query genes for enrichment from the 20,306 DESeq2
IDs. The query timestamp and full response are saved with the run; the full
response is Git-ignored because of its size.

## Results and figures

The results are partitioned by step:

- `results/01_deseq2/`: omnibus and pairwise DESeq2 tables, contrast summary,
  PCA, heat-response and top-gene heatmaps, the reproducible top-12 ranking,
  model diagnostics, and R session details.
- `results/02_go_enrichment/`: complete barley GO enrichment, highlighted
  response-to-heat term, query details, and heat-gene intersection.
- `results/03_heat_stress/`: 46- and 12-gene annotation tables, annotation
  coverage summaries, and a record of genes without descriptions.
- `docs/figures/`: small PNG report figures embedded in this README.

All full statistical tables and derived run outputs are regenerated during the
protocol and ignored by Git; the selected report images and analysis scripts
are tracked.
The compact tracked result tables are
[`analysis_summary.csv`](results/01_deseq2/analysis_summary.csv),
[`pairwise_contrast_summary.csv`](results/01_deseq2/pairwise_contrast_summary.csv),
[`heat_response_top12_gene_ids.csv`](results/01_deseq2/heat_response_top12_gene_ids.csv),
[`heat_response_46_annotations.csv`](results/03_heat_stress/heat_response_46_annotations.csv),
and [`heat_response_top12_annotations.csv`](results/03_heat_stress/heat_response_top12_annotations.csv).
The fresh annotation step retained all 46 gene IDs and found descriptions for
all 46 by combining g:Profiler and BARLEX annotations.

### Expression heatmap: 46 heat-response genes

![Heat-response gene expression across the five populations](docs/figures/heat_response_46_populations.png)

Columns are samples ordered by population; rows are hierarchically clustered
genes. Color represents per-gene centered DESeq2 VST expression.

### Expression heatmap: 12 reproducibly ranked candidates

![Top 12 heat-response genes ranked by population-mean variation](docs/figures/heat_response_top12_populations.png)

The focused set is selected by descending variation among the five population
mean expression values, not by manually selecting genes from the plot.

### Principal component analysis

![PCA of the 60 samples](docs/figures/pca_five_populations.png)

The figure is regenerated from the VST matrix by
`scripts/01_deseq2_populations.R`.

### Barley GO enrichment

![Barley GO enrichment of population-variable genes](docs/figures/go_enrichment.png)

The full g:Profiler response and a second plot highlighting GO:0009408 are
written under `results/02_go_enrichment/`.

## Output interpretation and limitations

The omnibus LRT answers whether expression varies across the five populations;
it does not identify a single population pair. Use the ten Wald contrast tables
to interpret specific pairwise differences. A significant adjusted p-value is
evidence of association with population under this design, not proof that heat
stress caused the difference. Heat-related genes are an enrichment-guided
candidate set; annotations and enrichment depend on the date/version of
g:Profiler and the supplied BARLEX release.
