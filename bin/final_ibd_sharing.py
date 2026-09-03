#!/usr/bin/env python3
"""Summarize and plot IBD sharing for final post-FST clusters."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from array import array
from collections import Counter, defaultdict
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
    parser.add_argument("--edges", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--weight-column", required=True)
    parser.add_argument("--output-prefix", required=True, type=Path)
    return parser.parse_args()


def read_membership(path: Path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Membership file has no header.")
        id_column = "ID" if "ID" in reader.fieldnames else reader.fieldnames[0]
        if "final_cluster" not in reader.fieldnames:
            raise ValueError("Membership file lacks final_cluster.")
        membership = {}
        for row in reader:
            sample_id = row[id_column]
            cluster = row["final_cluster"]
            if not sample_id or not cluster:
                raise ValueError("Empty sample ID or final cluster label.")
            if sample_id in membership:
                raise ValueError(f"Duplicate sample ID: {sample_id}")
            membership[sample_id] = cluster
    if not membership:
        raise ValueError("Membership file contains no samples.")
    return membership


def identify_edge_columns(fieldnames: list[str]):
    candidate_pairs = [
        ("ID1", "ID2"), ("id1", "id2"), ("sample_1", "sample_2"),
        ("sample1", "sample2"), ("iid1", "iid2"),
    ]
    for first, second in candidate_pairs:
        if first in fieldnames and second in fieldnames:
            return first, second
    if len(fieldnames) >= 2:
        return fieldnames[0], fieldnames[1]
    raise ValueError("Cannot identify the two sample-ID columns in the edge table.")


def write_matrix(path: Path, clusters: list[str], matrix: np.ndarray):
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["cluster", *clusters])
        for cluster, values in zip(clusters, matrix):
            writer.writerow([cluster, *[f"{value:.10g}" for value in values]])


def heatmap(path: Path, clusters: list[str], matrix: np.ndarray, title: str):
    size = max(8.0, 0.48 * len(clusters) + 3.0)
    fig, ax = plt.subplots(figsize=(size, size))
    image = ax.imshow(matrix, cmap="YlOrRd", aspect="equal")
    ax.set_xticks(range(len(clusters)), clusters, rotation=90)
    ax.set_yticks(range(len(clusters)), clusters)
    ax.set_title(title)
    finite = matrix[np.isfinite(matrix)]
    midpoint = (finite.min() + finite.max()) / 2 if finite.size else 0
    for row in range(len(clusters)):
        for column in range(len(clusters)):
            value = matrix[row, column]
            if math.isfinite(value):
                ax.text(
                    column, row, f"{value:.3f}", ha="center", va="center",
                    fontsize=max(4, 9 - len(clusters) // 4),
                    color="black" if value < midpoint else "white"
                )
    colorbar = fig.colorbar(image, ax=ax, shrink=0.8)
    colorbar.set_label("Mean total IBD length per possible pair (cM)")
    fig.tight_layout()
    fig.savefig(path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def boxplot(
    path: Path,
    clusters: list[str],
    sizes: Counter,
    values: dict[str, array],
):
    ordered = sorted(
        clusters,
        key=lambda cluster: (
            -float(np.mean(values[cluster])) if values[cluster] else math.inf,
            natural_key(cluster),
        ),
    )
    plotted = [cluster for cluster in ordered if values[cluster]]
    if not plotted:
        raise ValueError("No detected within-cluster IBD pairs were found.")

    height = max(7.0, 0.55 * len(plotted) + 2.5)
    fig, ax = plt.subplots(figsize=(11.0, height))
    data = [np.asarray(values[cluster], dtype=float) for cluster in plotted]
    labels = [f"{cluster} (n={sizes[cluster]})" for cluster in plotted]
    result = ax.boxplot(
        data,
        vert=False,
        tick_labels=labels,
        showfliers=False,
        patch_artist=True,
    )
    for box in result["boxes"]:
        box.set_facecolor("#9e9e9e")
        box.set_edgecolor("#333333")
    means = [float(np.mean(group)) for group in data]
    ax.scatter(means, range(1, len(plotted) + 1), color="black", s=28, zorder=3)
    ax.invert_yaxis()
    ax.set_xlabel("Within-cluster IBD sharing per detected sample pair (cM)")
    ax.set_ylabel("Final IBD clusters")
    ax.set_title("Distribution of within-cluster IBD sharing")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    membership = read_membership(args.membership)
    sizes = Counter(membership.values())
    clusters = sorted(sizes, key=natural_key)

    pair_sum = defaultdict(float)
    pair_detected = Counter()
    within_pair_values = defaultdict(lambda: array("d"))

    with open_text(args.edges) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Edge table has no header.")
        first_column, second_column = identify_edge_columns(reader.fieldnames)
        if args.weight_column not in reader.fieldnames:
            raise ValueError(
                f"Weight column {args.weight_column!r} not found. Available: "
                + ", ".join(reader.fieldnames)
            )
        for row in reader:
            first, second = row[first_column], row[second_column]
            if first not in membership or second not in membership:
                continue
            try:
                weight = float(row[args.weight_column])
            except (TypeError, ValueError):
                raise ValueError(f"Invalid IBD weight for pair {first}, {second}.")
            if not math.isfinite(weight) or weight < 0:
                raise ValueError(f"Non-finite or negative IBD weight for {first}, {second}.")
            first_cluster, second_cluster = membership[first], membership[second]
            key = tuple(sorted((first_cluster, second_cluster), key=natural_key))
            pair_sum[key] += weight
            pair_detected[key] += 1
            if first_cluster == second_cluster:
                within_pair_values[first_cluster].append(weight)

    prefix = args.output_prefix
    long_path = Path(f"{prefix}.ibd_sharing.tsv")
    matrix = np.zeros((len(clusters), len(clusters)), dtype=float)
    cluster_index = {cluster: index for index, cluster in enumerate(clusters)}
    with long_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "cluster_1", "cluster_2", "relationship", "cluster_1_size",
            "cluster_2_size", "possible_pairs", "detected_pairs",
            "detected_pair_fraction", "total_ibd_cm",
            "mean_ibd_cm_all_possible_pairs", "mean_ibd_cm_detected_pairs",
        ])
        for i, first_cluster in enumerate(clusters):
            for second_cluster in clusters[i:]:
                key = (first_cluster, second_cluster)
                if first_cluster == second_cluster:
                    possible = sizes[first_cluster] * (sizes[first_cluster] - 1) // 2
                    relationship = "within"
                else:
                    possible = sizes[first_cluster] * sizes[second_cluster]
                    relationship = "between"
                detected = pair_detected[key]
                total = pair_sum[key]
                all_mean = total / possible if possible else math.nan
                detected_mean = total / detected if detected else 0.0
                fraction = detected / possible if possible else math.nan
                writer.writerow([
                    first_cluster, second_cluster, relationship,
                    sizes[first_cluster], sizes[second_cluster], possible, detected,
                    f"{fraction:.10g}", f"{total:.10g}", f"{all_mean:.10g}",
                    f"{detected_mean:.10g}",
                ])
                first_index = cluster_index[first_cluster]
                second_index = cluster_index[second_cluster]
                matrix[first_index, second_index] = all_mean
                matrix[second_index, first_index] = all_mean

    distribution_path = Path(f"{prefix}.within_cluster_distribution.tsv")
    with distribution_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "final_cluster", "cluster_size", "detected_within_pairs",
            "mean_ibd_cm", "minimum_ibd_cm", "q1_ibd_cm", "median_ibd_cm",
            "q3_ibd_cm", "maximum_ibd_cm",
        ])
        for cluster in clusters:
            values = np.asarray(within_pair_values[cluster], dtype=float)
            if values.size:
                minimum = float(np.min(values))
                q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
                maximum = float(np.max(values))
                mean = float(np.mean(values))
                statistics = [mean, minimum, q1, median, q3, maximum]
            else:
                statistics = [math.nan] * 6
            writer.writerow([
                cluster, sizes[cluster], values.size,
                *[f"{value:.10g}" for value in statistics],
            ])

    write_matrix(Path(f"{prefix}.ibd_mean_matrix.tsv"), clusters, matrix)
    heatmap(Path(f"{prefix}.ibd_heatmap.png"), clusters, matrix,
            "Mean IBD sharing between final clusters")
    boxplot(
        Path(f"{prefix}.within_ibd_boxplot.png"),
        clusters,
        sizes,
        within_pair_values,
    )


if __name__ == "__main__":
    main()