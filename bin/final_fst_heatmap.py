#!/usr/bin/env python3
"""Recompute final pairwise Hudson FST and draw an annotated heatmap."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def natural_key(label: str):
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", label)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pgen", required=True, type=Path)
    parser.add_argument("--pvar", required=True, type=Path)
    parser.add_argument("--psam", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--plink2", default="plink2")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--memory-mb", type=int, default=2000)
    parser.add_argument("--output-prefix", required=True, type=Path)
    return parser.parse_args()


def write_phenotype(membership: Path, phenotype: Path):
    clusters = set()
    with open_text(membership) as source, phenotype.open(
        "wt", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Membership file has no header.")
        id_column = "ID" if "ID" in reader.fieldnames else reader.fieldnames[0]
        if "final_cluster" not in reader.fieldnames:
            raise ValueError("Membership file lacks final_cluster.")
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(["#IID", "FINAL_CLUSTER"])
        for row in reader:
            writer.writerow([row[id_column], row["final_cluster"]])
            clusters.add(row["final_cluster"])
    if len(clusters) < 2:
        raise ValueError("At least two final clusters are required for pairwise FST.")
    return sorted(clusters, key=natural_key)


def locate_fst(prefix: Path):
    for suffix in (".fst.summary", ".fst.summary.zst"):
        candidate = Path(f"{prefix}{suffix}")
        if candidate.exists():
            return candidate
    raise FileNotFoundError("PLINK 2 did not create an FST summary.")


def read_lines(path: Path):
    if path.suffix == ".zst":
        return subprocess.run(
            ["zstdcat", str(path)], check=True, text=True, capture_output=True
        ).stdout.splitlines()
    return path.read_text(encoding="utf-8").splitlines()


def main() -> None:
    args = parse_args()
    if shutil.which(args.plink2) is None:
        raise FileNotFoundError(args.plink2)
    phenotype = Path(f"{args.output_prefix}.final_clusters.pheno.tsv")
    clusters = write_phenotype(args.membership, phenotype)
    plink_prefix = Path(f"{args.output_prefix}.final_fst")
    pfile_prefix = args.pgen.with_suffix("")
    subprocess.run([
        args.plink2, "--pfile", str(pfile_prefix), "vzs",
        "--pheno", str(phenotype), "--pheno-name", "FINAL_CLUSTER",
        "--fst", "FINAL_CLUSTER", "method=hudson",
        "--threads", str(args.threads), "--memory", str(args.memory_mb),
        "--out", str(plink_prefix),
    ], check=True)

    reader = csv.DictReader(read_lines(locate_fst(plink_prefix)), delimiter="\t")
    if not reader.fieldnames:
        raise ValueError("FST output has no header.")
    pop1 = "#POP1" if "#POP1" in reader.fieldnames else "POP1"
    fst_columns = [name for name in reader.fieldnames if name.upper().endswith("FST")]
    if pop1 not in reader.fieldnames or "POP2" not in reader.fieldnames or not fst_columns:
        raise ValueError("Cannot identify population or FST columns in PLINK output.")
    fst_column = fst_columns[-1]
    values = []
    for row in reader:
        try:
            value = float(row[fst_column])
        except (TypeError, ValueError):
            value = math.nan
        values.append((row[pop1], row["POP2"], value))

    long_path = Path(f"{args.output_prefix}.final_fst.tsv")
    with long_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["cluster_1", "cluster_2", "hudson_fst"])
        for first, second, value in values:
            writer.writerow([first, second, f"{value:.10g}"])

    index = {cluster: i for i, cluster in enumerate(clusters)}
    matrix = np.full((len(clusters), len(clusters)), np.nan)
    np.fill_diagonal(matrix, 0.0)
    for first, second, value in values:
        matrix[index[first], index[second]] = value
        matrix[index[second], index[first]] = value

    matrix_path = Path(f"{args.output_prefix}.final_fst_matrix.tsv")
    with matrix_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["cluster", *clusters])
        for cluster, row in zip(clusters, matrix):
            writer.writerow([cluster, *["NA" if math.isnan(x) else f"{x:.10g}" for x in row]])

    size = max(8.0, 0.58 * len(clusters) + 3.0)
    fig, ax = plt.subplots(figsize=(size, size))
    image = ax.imshow(matrix, cmap="YlOrRd", aspect="equal")
    ax.set_xticks(range(len(clusters)), clusters, rotation=90)
    ax.set_yticks(range(len(clusters)), clusters)
    ax.set_title("Pairwise Hudson FST between final clusters")
    finite = matrix[np.isfinite(matrix)]
    midpoint = (finite.min() + finite.max()) / 2 if finite.size else 0 # midpoint = 0.07
    for row in range(len(clusters)):
        for column in range(len(clusters)):
            value = matrix[row, column]
            if math.isfinite(value):
                ax.text(column, row, f"{value:.4f}", ha="center", va="center",
                        fontsize=max(4, 9 - len(clusters) // 4),
                        color="black" if value < midpoint else "white")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.8)
    colorbar.set_label("Hudson FST")
    fig.tight_layout()
    fig.savefig(f"{args.output_prefix}.final_fst_heatmap.png", dpi=600,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
