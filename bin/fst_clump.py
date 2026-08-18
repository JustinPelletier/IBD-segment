#!/usr/bin/env python3
"""Iteratively merge graph communities whose recalculated Hudson FST is low."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pgen", required=True, type=Path)
    parser.add_argument("--pvar", required=True, type=Path)
    parser.add_argument("--psam", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--cluster-column", required=True)
    parser.add_argument("--threshold", required=True, type=float)
    parser.add_argument("--min-cluster-size", required=True, type=int)
    parser.add_argument("--plink2", default="plink2")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--memory-mb", type=int, default=2000)
    parser.add_argument("--output-prefix", required=True, type=Path)
    return parser.parse_args()


def read_membership(path: Path, cluster_column: str):
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Membership file has no header.")
        id_column = "ID" if "ID" in reader.fieldnames else reader.fieldnames[0]
        if cluster_column not in reader.fieldnames:
            raise ValueError(
                f"Membership column {cluster_column!r} was not found. "
                f"Available columns: {', '.join(reader.fieldnames)}"
            )
        rows = list(reader)

    if not rows:
        raise ValueError("Membership file contains no samples.")

    ids = [row[id_column] for row in rows]
    if any(not sample_id for sample_id in ids):
        raise ValueError("Membership file contains an empty sample ID.")
    if len(ids) != len(set(ids)):
        raise ValueError("Membership file contains duplicate sample IDs.")

    assignments = {row[id_column]: row[cluster_column] for row in rows}
    if any(not cluster for cluster in assignments.values()):
        raise ValueError(f"Column {cluster_column!r} contains an empty label.")
    return rows, reader.fieldnames, id_column, assignments


def read_psam_ids(path: Path) -> set[str]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("PSAM file has no header.")
        iid_column = "IID" if "IID" in reader.fieldnames else "#IID"
        if iid_column not in reader.fieldnames:
            raise ValueError("PSAM file does not contain IID or #IID.")
        ids = [row[iid_column] for row in reader]
    if len(ids) != len(set(ids)):
        raise ValueError("PSAM file contains duplicate IIDs.")
    return set(ids)


def write_phenotype(path: Path, assignments: dict[str, str]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["#IID", "CLUSTER"])
        for sample_id in sorted(assignments):
            writer.writerow([sample_id, assignments[sample_id]])


def find_fst_output(prefix: Path) -> Path:
    candidates = [
        Path(f"{prefix}.fst.summary"),
        Path(f"{prefix}.fst.summary.zst"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"PLINK 2 did not create {prefix}.fst.summary.")


def read_fst(path: Path) -> list[tuple[str, str, float]]:
    if path.suffix == ".zst":
        command = ["zstdcat", str(path)]
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        lines = completed.stdout.splitlines()
    else:
        lines = path.read_text(encoding="utf-8").splitlines()

    reader = csv.DictReader(lines, delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError(f"FST output {path} has no header.")

    pop1_column = "#POP1" if "#POP1" in reader.fieldnames else "POP1"
    pop2_column = "POP2"
    fst_candidates = [name for name in reader.fieldnames if name.upper().endswith("FST")]
    if pop1_column not in reader.fieldnames or pop2_column not in reader.fieldnames:
        raise ValueError(f"Cannot identify population columns in {path}.")
    if not fst_candidates:
        raise ValueError(f"Cannot identify the FST column in {path}.")
    fst_column = fst_candidates[-1]

    results = []
    for row in reader:
        try:
            fst = float(row[fst_column])
        except (TypeError, ValueError):
            fst = math.nan
        results.append((row[pop1_column], row[pop2_column], fst))
    return results


def run_fst(args: argparse.Namespace, assignments: dict[str, str], work: Path, step: int):
    phenotype = work / f"step_{step}.pheno.tsv"
    output = work / f"step_{step}"
    write_phenotype(phenotype, assignments)

    pfile_prefix = args.pgen.with_suffix("")
    command = [
        args.plink2,
        "--pfile", str(pfile_prefix), "vzs",
        "--pheno", str(phenotype),
        "--pheno-name", "CLUSTER",
        "--fst", "CLUSTER", "method=hudson",
        "--threads", str(args.threads),
        "--memory", str(args.memory_mb),
        "--out", str(output),
    ]
    subprocess.run(command, check=True)
    return read_fst(find_fst_output(output))


def main() -> None:
    args = parse_args()
    if args.threshold < 0:
        raise ValueError("--threshold must be nonnegative.")
    if args.min_cluster_size < 2:
        raise ValueError("--min-cluster-size must be at least two.")
    if args.threads < 1 or args.memory_mb < 1:
        raise ValueError("--threads and --memory-mb must be positive.")
    if shutil.which(args.plink2) is None:
        raise FileNotFoundError(f"PLINK 2 executable not found: {args.plink2}")

    for path in (args.pgen, args.pvar, args.psam, args.membership):
        if not path.is_file():
            raise FileNotFoundError(path)

    rows, fieldnames, id_column, assignments = read_membership(
        args.membership, args.cluster_column
    )
    psam_ids = read_psam_ids(args.psam)
    membership_ids = set(assignments)
    if membership_ids != psam_ids:
        missing_genotypes = sorted(membership_ids - psam_ids)
        missing_membership = sorted(psam_ids - membership_ids)
        raise ValueError(
            "Membership and genotype sample IDs differ: "
            f"{len(missing_genotypes)} membership IDs lack genotypes; "
            f"{len(missing_membership)} genotype IDs lack membership."
        )

    pairwise_rows: list[list[object]] = []
    merge_rows: list[list[object]] = []
    initial_cluster_count = len(set(assignments.values()))
    step = 0

    with tempfile.TemporaryDirectory(prefix="fst_clump_") as temporary_directory:
        work = Path(temporary_directory)
        while True:
            sizes = Counter(assignments.values())
            if len(sizes) < 2:
                break
            estimates = run_fst(args, assignments, work, step)
            eligible_pairs = []
            for cluster1, cluster2, fst in estimates:
                eligible = (
                    sizes[cluster1] >= args.min_cluster_size
                    and sizes[cluster2] >= args.min_cluster_size
                    and math.isfinite(fst)
                )
                if eligible:
                    eligible_pairs.append((fst, cluster1, cluster2))
                pairwise_rows.append([
                    step, cluster1, cluster2, fst, sizes[cluster1], sizes[cluster2],
                    "true" if eligible else "false",
                ])

            if not eligible_pairs:
                break
            fst, cluster1, cluster2 = min(
                eligible_pairs, key=lambda value: (value[0], value[1], value[2])
            )
            if fst >= args.threshold:
                break

            merged_label = min(cluster1, cluster2)
            removed_label = max(cluster1, cluster2)
            for sample_id, cluster in list(assignments.items()):
                if cluster == removed_label:
                    assignments[sample_id] = merged_label
            merge_rows.append([
                step + 1, cluster1, cluster2, fst, merged_label,
                sizes[cluster1], sizes[cluster2], sizes[cluster1] + sizes[cluster2],
            ])
            step += 1

    final_sizes = Counter(assignments.values())
    ordered_clusters = sorted(final_sizes, key=lambda label: (-final_sizes[label], label))
    final_labels = {
        cluster: f"Cluster_{index}" for index, cluster in enumerate(ordered_clusters, 1)
    }

    pairwise_path = Path(f"{args.output_prefix}.pairwise_fst.tsv.gz")
    with gzip.open(pairwise_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "iteration", "cluster_1", "cluster_2", "hudson_fst",
            "cluster_1_size", "cluster_2_size", "eligible_for_clumping",
        ])
        writer.writerows(pairwise_rows)

    history_path = Path(f"{args.output_prefix}.fst_merge_history.tsv")
    with history_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "step", "cluster_1", "cluster_2", "fst_before_merge", "retained_label",
            "cluster_1_size", "cluster_2_size", "merged_size",
        ])
        writer.writerows(merge_rows)

    final_path = Path(f"{args.output_prefix}.final_membership.tsv.gz")
    output_fields = list(fieldnames) + ["pre_fst_cluster", "final_cluster"]
    with gzip.open(final_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            sample_id = row[id_column]
            original_cluster = row[args.cluster_column]
            current_cluster = assignments[sample_id]
            output_row = dict(row)
            output_row["pre_fst_cluster"] = original_cluster
            output_row["final_cluster"] = final_labels[current_cluster]
            writer.writerow(output_row)

    summary_path = Path(f"{args.output_prefix}.fst_clumping_summary.tsv")
    with summary_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "samples", "initial_clusters", "merges", "final_clusters",
            "fst_threshold", "minimum_cluster_size",
        ])
        writer.writerow([
            len(assignments), initial_cluster_count, len(merge_rows), len(final_sizes),
            args.threshold, args.min_cluster_size,
        ])


if __name__ == "__main__":
    main()
