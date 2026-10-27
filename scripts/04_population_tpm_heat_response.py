#!/usr/bin/env python3

import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
TOP_GENES = ROOT / "results/01_deseq2/heat_response_top12_gene_ids.csv"
SAMPLES = ROOT / "data/metadata/samples.csv"
TX2GENE = ROOT / "references/MorexV3_tx2gene.tsv"
KALLISTO = ROOT / "results/upstream/kallisto"
OUTPUT = ROOT / "results/04_tpm_leaf_baseline"
PLOT_DIR = OUTPUT / "plots"

POPULATION_COLORS = {
    "k1": "#4477AA",
    "k5": "#EE6677",
    "k10": "#228833",
    "k12": "#CCBB44",
    "k14": "#66CCEE",
}


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    top_rows = read_csv(TOP_GENES)
    if len(top_rows) != 12:
        raise ValueError(f"Expected 12 ranked genes in {TOP_GENES}; found {len(top_rows)}")
    gene_ids = [row["gene_id"] for row in top_rows]
    ranks = [int(row["rank"]) for row in top_rows]
    if len(set(gene_ids)) != len(gene_ids):
        raise ValueError("The top-12 gene list contains duplicate gene IDs")
    if ranks != list(range(1, 13)):
        raise ValueError("The top-12 gene list must have ranks 1 through 12 in order")
    gene_rank = dict(zip(gene_ids, ranks))

    sample_rows = read_csv(SAMPLES)
    samples = [row["sample"] for row in sample_rows]
    if len(samples) != len(set(samples)):
        raise ValueError(f"Duplicate sample IDs in {SAMPLES}")
    sample_population = {row["sample"]: row["subpopulation"] for row in sample_rows}

    tx2gene = {}
    with TX2GENE.open(newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if len(row) != 2:
                raise ValueError(f"Expected 2 columns in {TX2GENE} line {line_number}")
            transcript_id, gene_id = row
            previous = tx2gene.setdefault(transcript_id, gene_id)
            if previous != gene_id:
                raise ValueError(f"Transcript {transcript_id} maps to multiple genes")

    missing_genes = set(gene_ids) - set(tx2gene.values())
    if missing_genes:
        raise ValueError(f"Top genes missing from transcript map: {sorted(missing_genes)}")

    sample_values = {}
    for sample in samples:
        abundance_path = KALLISTO / sample / "abundance.tsv"
        if not abundance_path.is_file():
            raise FileNotFoundError(f"Missing Kallisto abundance file: {abundance_path}")
        gene_tpm = defaultdict(float)
        seen_transcripts = set()
        with abundance_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"target_id", "tpm"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise ValueError(f"Unexpected columns in {abundance_path}")
            for row in reader:
                transcript_id = row["target_id"]
                if transcript_id in seen_transcripts:
                    raise ValueError(f"Duplicate target {transcript_id} in {abundance_path}")
                seen_transcripts.add(transcript_id)
                gene_id = tx2gene.get(transcript_id)
                if gene_id is None:
                    raise ValueError(f"Transcript {transcript_id} is missing from {TX2GENE}")
                if gene_id in gene_rank:
                    tpm = float(row["tpm"])
                    if tpm < 0:
                        raise ValueError(f"Negative TPM for {transcript_id} in {abundance_path}")
                    gene_tpm[gene_id] += tpm
        if not seen_transcripts:
            raise ValueError(f"No transcript records found in {abundance_path}")
        for gene_id in gene_ids:
            sample_values[(sample, gene_id)] = gene_tpm[gene_id]

    OUTPUT.mkdir(parents=True, exist_ok=True)
    sample_output = OUTPUT / "top12_gene_sample_tpm.csv"
    with sample_output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["rank", "gene_id", "sample", "population", "gene_tpm"],
        )
        writer.writeheader()
        for gene_id in gene_ids:
            for sample in samples:
                writer.writerow(
                    {
                        "rank": gene_rank[gene_id],
                        "gene_id": gene_id,
                        "sample": sample,
                        "population": sample_population[sample],
                        "gene_tpm": f"{sample_values[(sample, gene_id)]:.10g}",
                    }
                )

    summary_output = OUTPUT / "top12_gene_population_tpm_summary.csv"
    with summary_output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "gene_id",
                "population",
                "sample_count",
                "mean_gene_tpm",
                "sd_gene_tpm",
                "median_gene_tpm",
            ],
        )
        writer.writeheader()
        populations = sorted(
            set(sample_population.values()), key=lambda population: int(population[1:])
        )
        for gene_id in gene_ids:
            for population in populations:
                values = [
                    sample_values[(sample, gene_id)]
                    for sample in samples
                    if sample_population[sample] == population
                ]
                if not values:
                    raise ValueError(f"No samples found for population {population}")
                writer.writerow(
                    {
                        "rank": gene_rank[gene_id],
                        "gene_id": gene_id,
                        "population": population,
                        "sample_count": len(values),
                        "mean_gene_tpm": f"{statistics.mean(values):.10g}",
                        "sd_gene_tpm": f"{statistics.stdev(values):.10g}"
                        if len(values) > 1
                        else "",
                        "median_gene_tpm": f"{statistics.median(values):.10g}",
                    }
                )

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plot_samples = sorted(
        samples,
        key=lambda sample: (
            int(sample_population[sample][1:]),
            int(sample.split("M", 1)[0][1:]),
            int(sample.split("M", 1)[1]),
        ),
    )
    genotypes = list(
        dict.fromkeys(sample.split("M", 1)[0] for sample in plot_samples)
    )
    genotype_population = {
        genotype: sample_population[
            next(sample for sample in plot_samples if sample.startswith(genotype + "M"))
        ]
        for genotype in genotypes
    }
    genotype_positions = {
        genotype: position for position, genotype in enumerate(genotypes, start=1)
    }
    population_ranges = {}
    cursor = 1
    for population in populations:
        group_genotypes = [
            genotype for genotype in genotypes
            if genotype_population[genotype] == population
        ]
        population_ranges[population] = (
            cursor,
            cursor + len(group_genotypes) - 1,
        )
        cursor += len(group_genotypes)

    for gene_id in gene_ids:
        figure, axis = plt.subplots(figsize=(16, 7))
        for population in populations:
            start, end = population_ranges[population]
            group_genotypes = [
                genotype for genotype in genotypes
                if genotype_population[genotype] == population
            ]
            axis.axvspan(
                start - 0.5,
                end + 0.5,
                color=POPULATION_COLORS[population],
                alpha=0.07,
                zorder=0,
            )
            population_values = []
            for genotype in group_genotypes:
                genotype_samples = [
                    sample for sample in plot_samples
                    if sample.startswith(genotype + "M")
                ]
                replicate_offsets = {
                    sample: (index - (len(genotype_samples) - 1) / 2) * 0.2
                    for index, sample in enumerate(genotype_samples)
                }
                for sample in genotype_samples:
                    value = sample_values[(sample, gene_id)]
                    population_values.append(value)
                    x_position = genotype_positions[genotype] + replicate_offsets[sample]
                    axis.scatter(
                        x_position,
                        value,
                        color=POPULATION_COLORS[population],
                        edgecolor="white",
                        linewidth=0.5,
                        s=48,
                        zorder=3,
                    )
                    axis.annotate(
                        f"M{sample.split('M', 1)[1]}",
                        (x_position, value),
                        xytext=(0, 5),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=7,
                        color=POPULATION_COLORS[population],
                    )

            median = statistics.median(population_values)
            axis.hlines(
                median,
                start - 0.35,
                end + 0.35,
                color=POPULATION_COLORS[population],
                linewidth=2,
                zorder=2,
            )
            if population != populations[-1]:
                axis.axvline(end + 0.5, color="#888888", linewidth=0.8, alpha=0.7)
            axis.text(
                (start + end) / 2,
                1.015,
                population,
                transform=axis.get_xaxis_transform(),
                ha="center",
                va="bottom",
                color=POPULATION_COLORS[population],
                fontsize=10,
                fontweight="bold",
            )

        axis.set_xlim(0.5, len(genotypes) + 0.5)
        axis.set_xticks(
            [genotype_positions[genotype] for genotype in genotypes],
            genotypes,
            fontsize=8,
        )
        for label, genotype in zip(axis.get_xticklabels(), genotypes):
            label.set_color(POPULATION_COLORS[genotype_population[genotype]])
        axis.set_ylim(bottom=0)
        axis.set_xlabel("Genotype (replicate labels above points)")
        axis.set_ylabel("Gene-level TPM")
        axis.grid(axis="y", alpha=0.25)
        figure.suptitle(f"Rank {gene_rank[gene_id]}: {gene_id}", y=0.99)
        figure.tight_layout(rect=(0, 0, 1, 0.92))
        figure.savefig(
            PLOT_DIR / f"top12_rank_{gene_rank[gene_id]:02d}_{gene_id}.png",
            dpi=180,
        )
        plt.close(figure)

    print(
        f"Wrote TPM values for {len(gene_ids)} genes and {len(samples)} samples "
        f"to {sample_output} and {summary_output}; saved plots to {PLOT_DIR}."
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1)
