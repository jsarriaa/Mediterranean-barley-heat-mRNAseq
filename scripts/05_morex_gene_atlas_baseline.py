#!/usr/bin/env python3

import csv
import math
import re
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from html.parser import HTMLParser
from pathlib import Path
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
TOP_GENES = ROOT / "results/01_deseq2/heat_response_top12_gene_ids.csv"
POPULATION_SAMPLE_TPM = ROOT / "results/04_tpm_leaf_baseline/top12_gene_sample_tpm.csv"
POPULATION_SUMMARY = ROOT / "results/04_tpm_leaf_baseline/top12_gene_population_tpm_summary.csv"
EXTERNAL_DATA = ROOT / "data/external/morex_gene_atlas"
TPM_ARCHIVE = EXTERNAL_DATA / "hvmxrtd_all_tpm_files.zip"
EXPERIMENTS = EXTERNAL_DATA / "MorexAtlas_allProjects_paper_23Aug23_exp_ids.txt"
OUTPUT = ROOT / "results/04_tpm_leaf_baseline/morex_gene_atlas"
PLOT_DIR = OUTPUT / "plots"

RELATOR_URL = "https://ics.hutton.ac.uk/morexgeneatlas/relator.cgi"
TPM_ARCHIVE_URL = (
    "https://ics.hutton.ac.uk/morexgeneatlas/downloads/hvmxrtd_all_tpm_files.zip"
)
EXPERIMENTS_URL = (
    "https://ics.hutton.ac.uk/morexgeneatlas/downloads/"
    "MorexAtlas_allProjects_paper_23Aug23_exp_ids.txt"
)
USER_AGENT = "five-mediterranean-populations-rnaseq/1.0"

COMPARISON_LABELS = {
    "morex_1week_control": "1-week first-leaf control",
    "morex_2month_senescing": "2-month senescing leaf",
    "morex_30day_5thleaf_uninfested_2hr": (
        "30-day 5th-leaf uninfested control (2 hr)"
    ),
}

POPULATION_COLORS = {
    "k1": "#4477AA",
    "k5": "#EE6677",
    "k10": "#228833",
    "k12": "#CCBB44",
    "k14": "#66CCEE",
}

REFERENCE_COLORS = {
    "morex_1week_control": "#555555",
    "morex_2month_senescing": "#9673A6",
    "morex_30day_5thleaf_uninfested_2hr": "#D17A22",
}


def read_csv(path, delimiter=","):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def ensure_download(path, url):
    if path.is_file():
        if path.stat().st_size == 0:
            raise ValueError(f"Downloaded source file is empty: {path}")
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    partial_path = path.with_name(path.name + ".part")
    request = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request, timeout=180) as response, partial_path.open("wb") as output:
            shutil.copyfileobj(response, output)
        if partial_path.stat().st_size == 0:
            raise ValueError(f"Downloaded source file is empty: {url}")
        partial_path.replace(path)
    except (OSError, URLError, ValueError):
        partial_path.unlink(missing_ok=True)
        raise


class RelatorTableParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []
        self.current_row = None
        self.current_cell_text = None
        self.current_cell_titles = None

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag == "tr":
            self.current_row = []
        elif tag in {"td", "th"} and self.current_row is not None:
            self.current_cell_text = []
            self.current_cell_titles = []
        elif tag == "a" and self.current_cell_titles is not None:
            title = attributes.get("title", "").strip()
            if title:
                self.current_cell_titles.append(title)

    def handle_data(self, data):
        if self.current_cell_text is not None and data.strip():
            self.current_cell_text.append(data.strip())

    def handle_endtag(self, tag):
        if tag in {"td", "th"} and self.current_cell_text is not None:
            self.current_row.append(
                {
                    "text": " ".join(self.current_cell_text),
                    "titles": list(dict.fromkeys(self.current_cell_titles)),
                }
            )
            self.current_cell_text = None
            self.current_cell_titles = None
        elif tag == "tr" and self.current_row is not None:
            if self.current_row:
                self.rows.append(self.current_row)
            self.current_row = None


def fetch_relatormap(gene_ids):
    payload = urlencode({"list": "\n".join(gene_ids)}).encode("utf-8")
    request = Request(
        RELATOR_URL,
        data=payload,
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": USER_AGENT,
        },
        method="POST",
    )
    with urlopen(request, timeout=120) as response:
        html = response.read().decode("utf-8")

    parser = RelatorTableParser()
    parser.feed(html)
    requested = set(gene_ids)
    relations = []
    for row in parser.rows:
        if len(row) < 4:
            continue
        query_id = row[0]["text"]
        if query_id not in requested or query_id == "Query ID":
            continue
        related_id = row[1]["text"]
        class_text = row[3]["text"]
        relations.append(
            {
                "gene_id": query_id,
                "related_id": "" if related_id == "None" else related_id,
                "dataset": row[2]["text"],
                "class_text": class_text,
                "class_codes": list(dict.fromkeys(re.findall(r"\(([A-Za-z=])\)", class_text))),
                "class_descriptions": "; ".join(row[3]["titles"]),
            }
        )

    if not relations:
        raise ValueError("The Morex Gene Atlas relator returned no rows for the query")
    returned_genes = {row["gene_id"] for row in relations}
    missing_queries = requested - returned_genes
    if missing_queries:
        raise ValueError(
            "The Morex Gene Atlas relator omitted query IDs: "
            + ", ".join(sorted(missing_queries))
        )
    return relations


def select_atlas_gene_ids(gene_ids, relations):
    candidates = defaultdict(dict)
    for relation in relations:
        if relation["dataset"] != "HvMxRTD" or not relation["related_id"]:
            continue
        gene_id = relation["gene_id"]
        related_id = relation["related_id"]
        candidate = candidates[gene_id].setdefault(
            related_id,
            {"codes": set(), "descriptions": set(), "class_texts": set()},
        )
        candidate["codes"].update(relation["class_codes"])
        if relation["class_descriptions"]:
            candidate["descriptions"].add(relation["class_descriptions"])
        if relation["class_text"]:
            candidate["class_texts"].add(relation["class_text"])

    selected = {}
    statuses = {}
    policies = {}
    for gene_id in gene_ids:
        gene_candidates = candidates.get(gene_id, {})
        exact = [
            related_id
            for related_id, candidate in gene_candidates.items()
            if "=" in candidate["codes"]
        ]
        same_strand = [
            related_id
            for related_id, candidate in gene_candidates.items()
            if any(code not in {"x", "-"} for code in candidate["codes"])
        ]
        chosen = exact or same_strand
        selected[gene_id] = sorted(chosen)
        if exact:
            policies[gene_id] = "exact_intron_chain"
            statuses[gene_id] = (
                "single_exact_match" if len(exact) == 1 else "multiple_exact_matches"
            )
        elif same_strand:
            policies[gene_id] = "same_strand_nonexact_fallback"
            statuses[gene_id] = (
                "single_nonexact_match"
                if len(same_strand) == 1
                else "multiple_nonexact_matches"
            )
        elif gene_candidates:
            policies[gene_id] = "opposite_strand_only"
            statuses[gene_id] = "opposite_strand_only"
        else:
            policies[gene_id] = "no_HvMxRTD_relationship"
            statuses[gene_id] = "unmapped"

    return selected, statuses, policies


def read_experiments():
    rows = read_csv(EXPERIMENTS, delimiter="\t")
    required = {
        "expt_id",
        "sample_type",
        "study_accession",
        "run_accession_id",
        "sample_title",
        "cultivar",
        "tissue_type",
        "dev_stage",
        "bio_replicate",
        "sample_description",
        "expt_conditions",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Unexpected columns in Morex Atlas metadata: {EXPERIMENTS}")

    metadata_experiments = [row for row in rows if row.get("expt_id", "").strip()]
    ids = [row["expt_id"].strip() for row in metadata_experiments]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate expt_id values in Morex Atlas metadata")

    experiments = []
    for row in metadata_experiments:
        if row.get("cultivar", "").strip().lower() != "morex":
            continue
        tissue = row.get("tissue_type", "").strip().lower()
        stage = row.get("dev_stage", "").strip().lower()
        study = row.get("study_accession", "").strip()
        conditions = row.get("expt_conditions", "").strip().lower()
        sample_type = row.get("sample_type", "").strip().lower()
        description = row.get("sample_description", "").strip().lower()
        title = row.get("sample_title", "").strip()

        if study == "PRJEB52944":
            continue

        two_month_senescing = (
            study == "PRJEB14349"
            and tissue == "senescing leaf"
            and stage == "2 months"
        )
        thirty_day_uninfested = (
            study == "PRJNA326683"
            and tissue == "5th leaf, uninfested"
            and stage in {"30 days old", "30 days"}
            and sample_type == "control for treatment"
            and conditions == "non-infested"
            and "2hrs" in description
        )
        is_control_leaf = (
            tissue == "seedling leaf"
            and ("control" in sample_type or "control" in conditions)
        )
        if not (two_month_senescing or thirty_day_uninfested or is_control_leaf):
            continue

        primary = (
            study == "PRJEB29972"
            and conditions == "control"
            and "first leaf" in description
            and "no dnase" not in description
            and stage == "1 week"
        )
        comparison_group = ""
        if primary:
            role = "primary_baseline"
            comparison_group = "morex_1week_control"
        elif two_month_senescing:
            role = "two_month_senescing_technical"
            comparison_group = "morex_2month_senescing"
        elif thirty_day_uninfested:
            role = "30_day_5th_leaf_uninfested_2hr_control"
            comparison_group = "morex_30day_5thleaf_uninfested_2hr"
        elif "no dnase" in description:
            role = "no_dnase_control"
        elif "mock" in conditions:
            role = "mock_control"
        else:
            role = "other_control"

        row = dict(row)
        row["baseline_role"] = role
        row["comparison_group"] = comparison_group
        row["comparison_label"] = COMPARISON_LABELS.get(comparison_group, "")
        row["biological_replicate_id"] = row["bio_replicate"].strip()
        technical_match = re.search(r"_B(\d+)$", title)
        row["technical_replicate"] = (
            technical_match.group(1) if technical_match else ""
        )
        row["archive_member"] = (
            f"{row['study_accession']}_{row['run_accession_id']}_load.txt"
        )
        experiments.append(row)

    experiments.sort(key=lambda row: int(row["expt_id"]))
    primary_count = sum(
        row["baseline_role"] == "primary_baseline" for row in experiments
    )
    if primary_count != 3:
        raise ValueError(
            "Expected three 1-week, untreated first-leaf Morex replicates "
            f"from PRJEB29972; found {primary_count}. Review Atlas metadata."
        )
    expected_groups = {
        "morex_1week_control": 3,
        "morex_2month_senescing": 6,
        "morex_30day_5thleaf_uninfested_2hr": 4,
    }
    for group, expected in expected_groups.items():
        count = sum(row["comparison_group"] == group for row in experiments)
        if count != expected:
            raise ValueError(
                f"Expected {expected} Atlas samples in {group}; found {count}"
            )
    senescing_replicates = defaultdict(list)
    for row in experiments:
        if row["comparison_group"] == "morex_2month_senescing":
            senescing_replicates[row["biological_replicate_id"]].append(
                row["technical_replicate"]
            )
    if any(
        sorted(technical_replicates) != ["1", "2"]
        for technical_replicates in senescing_replicates.values()
    ) or len(senescing_replicates) != 3:
        raise ValueError(
            "Expected two technical batches for each of three 2-month "
            "senescing-leaf biological replicates"
        )
    spider_mite_replicates = {
        row["biological_replicate_id"]
        for row in experiments
        if row["comparison_group"] == "morex_30day_5thleaf_uninfested_2hr"
    }
    if spider_mite_replicates != {"1", "2", "3", "4"}:
        raise ValueError(
            "Expected four 30-day, 5th-leaf, uninfested 2-hour control "
            f"replicates; found {sorted(spider_mite_replicates)}"
        )
    if not experiments:
        raise ValueError("No Morex leaf experiments were found")
    return experiments


def list_archive_entries():
    result = subprocess.run(
        ["unzip", "-Z1", str(TPM_ARCHIVE)],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        raise RuntimeError(
            f"Could not list {TPM_ARCHIVE}: {result.stderr.strip() or result.stdout.strip()}"
        )
    return set(result.stdout.splitlines())


def validate_archive():
    result = subprocess.run(
        ["unzip", "-tq", str(TPM_ARCHIVE)],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        raise RuntimeError(
            f"Morex Atlas TPM archive validation failed: "
            f"{result.stderr.strip() or result.stdout.strip()}"
        )


def read_sample_gene_tpm(experiment, selected_ids, archive_entries):
    member = experiment["archive_member"]
    if member not in archive_entries:
        raise ValueError(
            f"Morex Atlas archive is missing {member} for expt_id "
            f"{experiment['expt_id']}"
        )

    result = subprocess.run(
        ["unzip", "-p", str(TPM_ARCHIVE), member],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        raise RuntimeError(
            f"Could not read {member}: {result.stderr.strip() or result.stdout.strip()}"
        )

    selected = set(selected_ids)
    values = defaultdict(float)
    observed_genes = set()
    transcript_ids = set()
    row_count = 0
    expected_expt_id = int(experiment["expt_id"])
    for line_number, line in enumerate(result.stdout.splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) != 6:
            raise ValueError(
                f"Expected 6 tab-separated fields in {member} line {line_number}; "
                f"found {len(fields)}"
            )
        expt_id, transcript_id, _, atlas_gene_id, tpm_text, _ = fields
        if int(expt_id) != expected_expt_id:
            raise ValueError(
                f"Unexpected expt_id {expt_id} in {member}; "
                f"expected {expected_expt_id}"
            )
        tpm = float(tpm_text)
        if not math.isfinite(tpm) or tpm < 0:
            raise ValueError(
                f"Invalid TPM value in {member} line {line_number}: {tpm_text}"
            )
        if atlas_gene_id in selected:
            if transcript_id in transcript_ids:
                raise ValueError(f"Duplicate transcript {transcript_id} in {member}")
            transcript_ids.add(transcript_id)
            observed_genes.add(atlas_gene_id)
            values[atlas_gene_id] += tpm
        row_count += 1

    if row_count == 0:
        raise ValueError(f"No TPM records found in {member}")
    return values, observed_genes, row_count


def write_mapping_audit(gene_ids, relations, selected, statuses, policies):
    relation_rows = []
    for relation in relations:
        gene_id = relation["gene_id"]
        related_id = relation["related_id"]
        used = relation["dataset"] == "HvMxRTD" and related_id in selected[gene_id]
        relation_rows.append(
            {
                "gene_id": gene_id,
                "related_id": related_id,
                "dataset": relation["dataset"],
                "class_text": relation["class_text"],
                "class_codes": "|".join(relation["class_codes"]),
                "class_descriptions": relation["class_descriptions"],
                "selected_for_expression": str(used).lower(),
                "selection_policy": policies[gene_id] if used else "",
                "gene_mapping_status": statuses[gene_id],
            }
        )

    relations_by_gene = {row["gene_id"] for row in relations}
    for gene_id in gene_ids:
        if gene_id not in relations_by_gene:
            relation_rows.append(
                {
                    "gene_id": gene_id,
                    "related_id": "",
                    "dataset": "",
                    "class_text": "",
                    "class_codes": "",
                    "class_descriptions": "",
                    "selected_for_expression": "false",
                    "selection_policy": "",
                    "gene_mapping_status": statuses[gene_id],
                }
            )

    write_csv(
        OUTPUT / "top12_morex_gene_atlas_relator_mapping.csv",
        [
            "gene_id",
            "related_id",
            "dataset",
            "class_text",
            "class_codes",
            "class_descriptions",
            "selected_for_expression",
            "selection_policy",
            "gene_mapping_status",
        ],
        relation_rows,
    )


def write_baseline_outputs(gene_ids, gene_rank, selected, statuses, controls, archive_entries):
    all_sample_rows = []
    for experiment in controls:
        selected_ids = sorted(
            {
                atlas_id
                for gene_id in gene_ids
                for atlas_id in selected[gene_id]
            }
        )
        tpm_values, observed_ids, row_count = read_sample_gene_tpm(
            experiment, selected_ids, archive_entries
        )
        for gene_id in gene_ids:
            target_ids = selected[gene_id]
            missing_ids = [atlas_id for atlas_id in target_ids if atlas_id not in observed_ids]
            if not target_ids:
                gene_tpm = ""
                value_status = statuses[gene_id]
            elif missing_ids:
                gene_tpm = ""
                value_status = "target_gene_missing_in_TPM_file"
            else:
                gene_tpm = f"{sum(tpm_values[atlas_id] for atlas_id in target_ids):.10g}"
                value_status = "complete"

            all_sample_rows.append(
                {
                    "rank": gene_rank[gene_id],
                    "gene_id": gene_id,
                    "expt_id": experiment["expt_id"],
                    "study_accession": experiment["study_accession"],
                    "run_accession_id": experiment["run_accession_id"],
                    "sample_title": experiment["sample_title"],
                    "sample_description": experiment["sample_description"],
                    "sample_type": experiment["sample_type"],
                    "cultivar": experiment["cultivar"],
                    "tissue_type": experiment["tissue_type"],
                    "dev_stage": experiment["dev_stage"],
                    "bio_replicate": experiment["bio_replicate"],
                    "expt_conditions": experiment["expt_conditions"],
                    "baseline_role": experiment["baseline_role"],
                    "comparison_group": experiment["comparison_group"],
                    "comparison_label": experiment["comparison_label"],
                    "biological_replicate_id": experiment[
                        "biological_replicate_id"
                    ],
                    "technical_replicate": experiment["technical_replicate"],
                    "selected_hvmxrtd_gene_ids": "|".join(target_ids),
                    "selected_hvmxrtd_gene_count": len(target_ids),
                    "mapping_status": statuses[gene_id],
                    "gene_tpm": gene_tpm,
                    "value_status": value_status,
                    "atlas_transcript_rows": row_count,
                }
            )

    write_csv(
        OUTPUT / "top12_morex_leaf_tpm_by_sample.csv",
        list(all_sample_rows[0]),
        all_sample_rows,
    )

    biological_rows = []
    comparison_groups = list(COMPARISON_LABELS)
    for gene_id in gene_ids:
        for comparison_group in comparison_groups:
            group_rows = [
                row
                for row in all_sample_rows
                if row["gene_id"] == gene_id
                and row["comparison_group"] == comparison_group
            ]
            replicates = defaultdict(list)
            for row in group_rows:
                replicates[row["biological_replicate_id"]].append(row)
            for replicate_id, replicate_rows in sorted(
                replicates.items(), key=lambda item: int(item[0])
            ):
                observed = [
                    float(row["gene_tpm"])
                    for row in replicate_rows
                    if row["value_status"] == "complete"
                ]
                complete = len(observed) == len(replicate_rows)
                biological_rows.append(
                    {
                        "rank": gene_rank[gene_id],
                        "gene_id": gene_id,
                        "comparison_group": comparison_group,
                        "comparison_label": COMPARISON_LABELS[comparison_group],
                        "biological_replicate_id": replicate_id,
                        "n_technical_replicates": len(replicate_rows),
                        "expt_ids": "|".join(
                            row["expt_id"] for row in replicate_rows
                        ),
                        "sample_titles": "|".join(
                            row["sample_title"] for row in replicate_rows
                        ),
                        "gene_tpm": (
                            f"{statistics.mean(observed):.10g}" if complete else ""
                        ),
                        "value_status": "complete" if complete else "incomplete",
                    }
                )
    write_csv(
        OUTPUT / "top12_morex_biological_replicate_tpm.csv",
        list(biological_rows[0]),
        biological_rows,
    )

    context_summary_rows = []
    for gene_id in gene_ids:
        for comparison_group in comparison_groups:
            group_rows = [
                row
                for row in biological_rows
                if row["gene_id"] == gene_id
                and row["comparison_group"] == comparison_group
            ]
            values = [
                float(row["gene_tpm"])
                for row in group_rows
                if row["value_status"] == "complete"
            ]
            context_summary_rows.append(
                {
                    "rank": gene_rank[gene_id],
                    "gene_id": gene_id,
                    "comparison_group": comparison_group,
                    "comparison_label": COMPARISON_LABELS[comparison_group],
                    "biological_replicate_count": len(values),
                    "mean_gene_tpm": (
                        f"{statistics.mean(values):.10g}" if values else ""
                    ),
                    "sd_gene_tpm": (
                        f"{statistics.stdev(values):.10g}"
                        if len(values) > 1
                        else ""
                    ),
                    "replicate_tpm_values": "|".join(
                        row["gene_tpm"]
                        for row in group_rows
                        if row["value_status"] == "complete"
                    ),
                    "technical_replicate_counts": "|".join(
                        str(row["n_technical_replicates"])
                        for row in group_rows
                        if row["value_status"] == "complete"
                    ),
                    "comparison_status": (
                        "complete" if len(values) == len(group_rows) else "missing"
                    ),
                }
            )
    write_csv(
        OUTPUT / "top12_morex_comparison_group_summary.csv",
        list(context_summary_rows[0]),
        context_summary_rows,
    )

    primary_experiments = [
        row
        for row in controls
        if row["comparison_group"] == "morex_1week_control"
    ]
    baseline_summary = []
    for gene_id in gene_ids:
        gene_values = []
        for row in all_sample_rows:
            if row["gene_id"] == gene_id and row["baseline_role"] == "primary_baseline":
                if row["value_status"] == "complete":
                    gene_values.append(float(row["gene_tpm"]))
        expected = len(primary_experiments)
        if gene_values:
            mean_value = statistics.mean(gene_values)
            sd_value = statistics.stdev(gene_values) if len(gene_values) > 1 else None
        else:
            mean_value = None
            sd_value = None
        baseline_summary.append(
            {
                "rank": gene_rank[gene_id],
                "gene_id": gene_id,
                "selected_hvmxrtd_gene_ids": "|".join(selected[gene_id]),
                "selected_hvmxrtd_gene_count": len(selected[gene_id]),
                "mapping_status": statuses[gene_id],
                "primary_expt_ids": "|".join(
                    row["expt_id"] for row in primary_experiments
                ),
                "baseline_expected_n": expected,
                "baseline_observed_n": len(gene_values),
                "baseline_mean_gene_tpm": (
                    f"{mean_value:.10g}" if mean_value is not None else ""
                ),
                "baseline_sd_gene_tpm": (
                    f"{sd_value:.10g}" if sd_value is not None else ""
                ),
                "baseline_status": (
                    "complete"
                    if len(gene_values) == expected
                    else "missing_mapping"
                    if not selected[gene_id]
                    else "incomplete"
                ),
            }
        )

    write_csv(
        OUTPUT / "top12_morex_primary_leaf_baseline_summary.csv",
        list(baseline_summary[0]),
        baseline_summary,
    )

    summary_by_gene = {row["gene_id"]: row for row in baseline_summary}
    population_summary = read_csv(POPULATION_SUMMARY)
    comparison_rows = []
    for row in population_summary:
        baseline = summary_by_gene[row["gene_id"]]
        comparison_rows.append(
            {
                "rank": row["rank"],
                "gene_id": row["gene_id"],
                "population": row["population"],
                "population_sample_count": row["sample_count"],
                "population_mean_gene_tpm": row["mean_gene_tpm"],
                "population_sd_gene_tpm": row["sd_gene_tpm"],
                "morex_baseline_mean_gene_tpm": baseline["baseline_mean_gene_tpm"],
                "morex_baseline_sd_gene_tpm": baseline["baseline_sd_gene_tpm"],
                "morex_baseline_n": baseline["baseline_observed_n"],
                "morex_baseline_expected_n": baseline["baseline_expected_n"],
                "mapping_status": baseline["mapping_status"],
                "baseline_status": baseline["baseline_status"],
            }
        )
    write_csv(
        OUTPUT / "top12_population_vs_morex_baseline.csv",
        list(comparison_rows[0]),
        comparison_rows,
    )
    extended_comparison_rows = []
    for row in population_summary:
        extended_comparison_rows.append(
            {
                "rank": row["rank"],
                "gene_id": row["gene_id"],
                "comparison_group": f"population_{row['population']}",
                "comparison_label": row["population"],
                "group_type": "population",
                "sample_count": row["sample_count"],
                "mean_gene_tpm": row["mean_gene_tpm"],
                "sd_gene_tpm": row["sd_gene_tpm"],
                "mapping_status": "",
            }
        )
    for row in context_summary_rows:
        extended_comparison_rows.append(
            {
                "rank": row["rank"],
                "gene_id": row["gene_id"],
                "comparison_group": row["comparison_group"],
                "comparison_label": row["comparison_label"],
                "group_type": "Morex Gene Atlas reference",
                "sample_count": row["biological_replicate_count"],
                "mean_gene_tpm": row["mean_gene_tpm"],
                "sd_gene_tpm": row["sd_gene_tpm"],
                "mapping_status": next(
                    item["mapping_status"]
                    for item in baseline_summary
                    if item["gene_id"] == row["gene_id"]
                ),
            }
        )
    write_csv(
        OUTPUT / "top12_population_vs_morex_context_groups.csv",
        list(extended_comparison_rows[0]),
        extended_comparison_rows,
    )
    return all_sample_rows, biological_rows


def make_comparison_plots(gene_ids, gene_rank, atlas_biological_rows):
    population_rows = read_csv(POPULATION_SAMPLE_TPM)
    if len(population_rows) != 12 * 60:
        raise ValueError(
            f"Expected 720 population sample TPM rows; found {len(population_rows)}"
        )
    population_by_gene = defaultdict(list)
    for row in population_rows:
        population_by_gene[row["gene_id"]].append(row)

    references_by_gene = defaultdict(lambda: defaultdict(list))
    for row in atlas_biological_rows:
        if row["value_status"] == "complete":
            references_by_gene[row["gene_id"]][row["comparison_group"]].append(row)

    def sample_sort_key(row):
        match = re.fullmatch(r"T(\d+)M(\d+)", row["sample"])
        if not match:
            raise ValueError(f"Unexpected population sample ID: {row['sample']}")
        return int(match.group(1)), int(match.group(2))

    population_rows.sort(key=sample_sort_key)
    genotypes = list(
        dict.fromkeys(
            re.fullmatch(r"(T\d+)M\d+", row["sample"]).group(1)
            for row in population_rows
        )
    )
    genotype_population = {}
    for row in population_rows:
        genotype = re.fullmatch(r"(T\d+)M\d+", row["sample"]).group(1)
        previous = genotype_population.setdefault(genotype, row["population"])
        if previous != row["population"]:
            raise ValueError(f"Genotype {genotype} is assigned to multiple populations")

    populations = sorted(
        set(genotype_population.values()), key=lambda value: int(value[1:])
    )
    positions = {genotype: index for index, genotype in enumerate(genotypes, start=1)}
    population_ranges = {}
    cursor = 1
    for population in populations:
        group = [g for g in genotypes if genotype_population[g] == population]
        population_ranges[population] = (cursor, cursor + len(group) - 1)
        cursor += len(group)

    reference_groups = list(COMPARISON_LABELS)
    reference_positions = {
        group: len(genotypes) + index
        for index, group in enumerate(reference_groups, start=1)
    }
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    for gene_id in gene_ids:
        figure, axis = plt.subplots(figsize=(21, 7))
        for population in populations:
            start, end = population_ranges[population]
            axis.axvspan(
                start - 0.5,
                end + 0.5,
                color=POPULATION_COLORS[population],
                alpha=0.07,
                zorder=0,
            )
            population_values = []
            for genotype in genotypes:
                if genotype_population[genotype] != population:
                    continue
                genotype_samples = [
                    row for row in population_by_gene[gene_id]
                    if row["sample"].startswith(genotype + "M")
                ]
                genotype_samples.sort(key=sample_sort_key)
                for index, row in enumerate(genotype_samples):
                    offset = (index - (len(genotype_samples) - 1) / 2) * 0.2
                    x_position = positions[genotype] + offset
                    value = float(row["gene_tpm"])
                    population_values.append(value)
                    axis.scatter(
                        x_position,
                        value,
                        color=POPULATION_COLORS[population],
                        edgecolor="white",
                        linewidth=0.5,
                        s=48,
                        zorder=3,
                    )
                    replicate = row["sample"].split("M", 1)[1]
                    axis.annotate(
                        f"M{replicate}",
                        (x_position, value),
                        xytext=(0, 5 + index * 7),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=7,
                        color=POPULATION_COLORS[population],
                    )

            if population_values:
                axis.hlines(
                    statistics.median(population_values),
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

        for comparison_group in reference_groups:
            position = reference_positions[comparison_group]
            color = REFERENCE_COLORS[comparison_group]
            axis.axvspan(
                position - 0.5,
                position + 0.5,
                color=color,
                alpha=0.08,
                zorder=0,
            )
            references = sorted(
                references_by_gene[gene_id].get(comparison_group, []),
                key=lambda row: int(row["biological_replicate_id"]),
            )
            values = [float(row["gene_tpm"]) for row in references]
            if values:
                offsets = {
                    row["biological_replicate_id"]: (
                        index - (len(references) - 1) / 2
                    )
                    * 0.2
                    for index, row in enumerate(references)
                }
                for index, row in enumerate(references):
                    value = float(row["gene_tpm"])
                    x_position = (
                        position + offsets[row["biological_replicate_id"]]
                    )
                    replicate_label = (
                        f"P{row['biological_replicate_id']}"
                        if comparison_group == "morex_1week_control"
                        else f"R{row['biological_replicate_id']}"
                    )
                    axis.scatter(
                        x_position,
                        value,
                        color=color,
                        edgecolor="white",
                        linewidth=0.5,
                        s=52,
                        zorder=4,
                    )
                    axis.annotate(
                        replicate_label,
                        (x_position, value),
                        xytext=(0, 5 + index * 7),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=7,
                        color=color,
                    )
                mean_value = statistics.mean(values)
                sd_value = statistics.stdev(values) if len(values) > 1 else 0
                axis.errorbar(
                    position,
                    mean_value,
                    yerr=sd_value,
                    fmt="D",
                    color="#222222",
                    markersize=5,
                    capsize=3,
                    linewidth=1.2,
                    zorder=5,
                )
            else:
                axis.text(
                    position,
                    0.03,
                    "No mapped\nTPM",
                    transform=axis.get_xaxis_transform(),
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    color="#8B0000",
                )

            if position != reference_positions[reference_groups[-1]]:
                axis.axvline(
                    position + 0.5,
                    color="#888888",
                    linewidth=0.6,
                    alpha=0.5,
                )
        xticks = [positions[genotype] for genotype in genotypes] + [
            reference_positions[group] for group in reference_groups
        ]
        reference_ticklabels = {
            "morex_1week_control": "1-week\nfirst-leaf\ncontrol",
            "morex_2month_senescing": "2-month\nsenescing\nleaf",
            "morex_30day_5thleaf_uninfested_2hr": (
                "30-day\n5th-leaf\nuninfested\ncontrol (2 hr)"
            ),
        }
        ticklabels = genotypes + [
            reference_ticklabels[group] for group in reference_groups
        ]
        axis.set_xlim(0.5, len(genotypes) + len(reference_groups) + 0.5)
        axis.set_xticks(xticks, ticklabels, fontsize=7)
        for label, genotype in zip(axis.get_xticklabels(), genotypes):
            label.set_color(POPULATION_COLORS[genotype_population[genotype]])
        for label, group in zip(
            axis.get_xticklabels()[len(genotypes):], reference_groups
        ):
            label.set_color(REFERENCE_COLORS[group])
        axis.set_ylim(bottom=0)
        axis.set_xlabel(
            "Population genotypes and distinct Morex Gene Atlas reference groups"
        )
        axis.set_ylabel("Gene-level TPM (linear scale)")
        axis.grid(axis="y", alpha=0.25)
        figure.suptitle(f"Rank {gene_rank[gene_id]}: {gene_id}", y=0.99)
        figure.tight_layout(rect=(0, 0, 1, 0.92))
        figure.savefig(
            PLOT_DIR / f"top12_rank_{gene_rank[gene_id]:02d}_{gene_id}.png",
            dpi=180,
        )
        plt.close(figure)


def main():
    top_rows = read_csv(TOP_GENES)
    if len(top_rows) != 12:
        raise ValueError(f"Expected 12 ranked genes in {TOP_GENES}; found {len(top_rows)}")
    gene_ids = [row["gene_id"] for row in top_rows]
    ranks = [int(row["rank"]) for row in top_rows]
    if len(set(gene_ids)) != 12 or ranks != list(range(1, 13)):
        raise ValueError("Top-gene input must contain unique ranks 1 through 12")
    gene_rank = dict(zip(gene_ids, ranks))

    if not POPULATION_SAMPLE_TPM.is_file() or not POPULATION_SUMMARY.is_file():
        raise FileNotFoundError(
            "Run scripts/04_population_tpm_heat_response.py before this script"
        )
    OUTPUT.mkdir(parents=True, exist_ok=True)
    ensure_download(TPM_ARCHIVE, TPM_ARCHIVE_URL)
    ensure_download(EXPERIMENTS, EXPERIMENTS_URL)

    relations = fetch_relatormap(gene_ids)
    selected, statuses, policies = select_atlas_gene_ids(gene_ids, relations)
    write_mapping_audit(
        gene_ids, relations, selected, statuses, policies
    )

    controls = read_experiments()
    validate_archive()
    archive_entries = list_archive_entries()
    missing_members = [
        row["archive_member"]
        for row in controls
        if row["archive_member"] not in archive_entries
    ]
    if missing_members:
        raise ValueError(
            "Morex Atlas archive is missing selected sample files: "
            + ", ".join(missing_members)
        )

    sample_rows, biological_rows = write_baseline_outputs(
        gene_ids,
        gene_rank,
        selected,
        statuses,
        controls,
        archive_entries,
    )
    make_comparison_plots(gene_ids, gene_rank, biological_rows)

    primary_count = sum(
        row["comparison_group"] == "morex_1week_control" for row in controls
    )
    experiment_count = len(controls)
    biological_comparison_groups = len(COMPARISON_LABELS)
    print(
        f"Processed {len(gene_ids)} genes across {experiment_count} Morex leaf "
        f"experiments, including {primary_count} primary baseline replicates "
        f"and {biological_comparison_groups} comparison groups. Mapping audit, "
        f"TPM tables, comparison tables, and plots are in {OUTPUT}."
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, URLError, ValueError, KeyError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1)
