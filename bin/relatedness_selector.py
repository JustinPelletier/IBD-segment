#!/usr/bin/env python3
"""Select an unrelated cohort from KING and Hap-IBD related pairs.

The script builds the union of:

* KING pairs with a kinship coefficient at or above the configured cutoff;
* Hap-IBD pairs with total genome-wide IBD strictly above the configured cM
  cutoff.

Two removal modes are supported. ``remove_all`` removes every sample incident
to a relatedness edge. ``optimal`` finds a minimum vertex cover within small
relatedness components, which is equivalent to retaining a maximum unrelated
set. Large or time-limited components use a deterministic maximum-degree
fallback and are identified in the audit output.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import math
import sys
import time
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


IBD_TOTAL_COLUMNS = (
    "total_IBD_length_cM",
    "total_ibd_length_cm",
    "total_ibd_cm",
)
SAMPLE_1_COLUMNS = ("sample_1", "sample1", "iid1", "id1")
SAMPLE_2_COLUMNS = ("sample_2", "sample2", "iid2", "id2")
KING_1_COLUMNS = ("iid1", "id1", "sample_1", "sample1")
KING_2_COLUMNS = ("iid2", "id2", "sample_2", "sample2")
KING_VALUE_COLUMNS = ("kinship", "king_kinship")


@dataclass
class PairEvidence:
    sample_1: str
    sample_2: str
    king_kinship: float | None = None
    total_ibd_cm: float | None = None

    @property
    def source(self) -> str:
        if self.king_kinship is not None and self.total_ibd_cm is not None:
            return "KING+IBD"
        if self.king_kinship is not None:
            return "KING"
        return "IBD"


class SearchTimeout(RuntimeError):
    """Raised when exact minimum-vertex-cover search exceeds its time limit."""


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine KING and Hap-IBD related pairs and select unrelated samples."
    )
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--king", required=True, type=Path)
    parser.add_argument("--ibd", required=True, type=Path)
    parser.add_argument("--king-cutoff", required=True, type=float)
    parser.add_argument("--ibd-total-cm-cutoff", required=True, type=float)
    parser.add_argument(
        "--removal-mode",
        required=True,
        choices=("optimal", "remove_all"),
    )
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument(
        "--exact-component-size-limit",
        type=int,
        default=60,
        help="Largest component submitted to exact optimization (default: 60).",
    )
    parser.add_argument(
        "--exact-search-seconds",
        type=float,
        default=60.0,
        help="Exact-search time limit per component (default: 60 seconds).",
    )
    args = parser.parse_args()

    if not 0 < args.king_cutoff < 0.5:
        parser.error("--king-cutoff must be greater than 0 and smaller than 0.5")
    if args.ibd_total_cm_cutoff <= 0:
        parser.error("--ibd-total-cm-cutoff must be greater than 0")
    if args.exact_component_size_limit < 2:
        parser.error("--exact-component-size-limit must be at least 2")
    if args.exact_search_seconds <= 0:
        parser.error("--exact-search-seconds must be greater than 0")
    return args


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def open_reproducible_gzip(path: Path) -> TextIO:
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return io.TextIOWrapper(compressed, encoding="utf-8", newline="")


def normalized_header(value: str) -> str:
    return value.strip().lstrip("#").lower()


def resolve_column(fieldnames: Iterable[str], candidates: Iterable[str], label: str) -> str:
    lookup = {normalized_header(name): name for name in fieldnames if name is not None}
    for candidate in candidates:
        match = lookup.get(normalized_header(candidate))
        if match is not None:
            return match
    raise ValueError(
        f"Unable to find {label}. Available columns: {', '.join(fieldnames)}"
    )


def detect_delimiter(header: str) -> str | None:
    return "\t" if "\t" in header else None


def dictionary_rows(path: Path) -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        first = handle.readline()
        if not first:
            raise ValueError(f"Input table is empty: {path}")
        delimiter = detect_delimiter(first)
        # Chain without materializing a potentially very large IBD table.
        def stream() -> Iterator[str]:
            yield first
            yield from handle

        if delimiter == "\t":
            reader = csv.DictReader(stream(), delimiter="\t")
            yield from reader
        else:
            headers = first.strip().split()
            for line_number, line in enumerate(handle, start=2):
                if not line.strip():
                    continue
                values = line.strip().split()
                if len(values) != len(headers):
                    raise ValueError(
                        f"Malformed whitespace-delimited row {line_number} in {path}: "
                        f"expected {len(headers)} fields, observed {len(values)}"
                    )
                yield dict(zip(headers, values))


def table_columns(path: Path) -> list[str]:
    with open_text(path) as handle:
        header = handle.readline().rstrip("\r\n")
    if not header:
        raise ValueError(f"Input table is empty: {path}")
    return header.split("\t") if "\t" in header else header.split()


def parse_float(value: str, path: Path, column: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Non-numeric value in {path}, column {column}: {value!r}") from error
    if not math.isfinite(result):
        raise ValueError(f"Non-finite value in {path}, column {column}: {value!r}")
    return result


def canonical_pair(sample_1: str, sample_2: str) -> tuple[str, str]:
    if not sample_1 or not sample_2:
        raise ValueError("Relatedness tables contain an empty sample ID")
    if sample_1 == sample_2:
        raise ValueError(f"Self-pair encountered for sample {sample_1!r}")
    return tuple(sorted((sample_1, sample_2)))


def read_samples(path: Path) -> list[str]:
    samples: list[str] = []
    seen: set[str] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            sample = line.strip()
            if not sample:
                continue
            if sample in seen:
                raise ValueError(
                    f"Duplicate sample ID {sample!r} in {path} at line {line_number}"
                )
            seen.add(sample)
            samples.append(sample)
    if not samples:
        raise ValueError(f"Sample list is empty: {path}")
    return samples


def add_king_pairs(
    path: Path,
    cutoff: float,
    cohort: set[str],
    evidence: dict[tuple[str, str], PairEvidence],
) -> int:
    columns = table_columns(path)
    sample_1_col = resolve_column(columns, KING_1_COLUMNS, "KING sample-1 column")
    sample_2_col = resolve_column(columns, KING_2_COLUMNS, "KING sample-2 column")
    kinship_col = resolve_column(columns, KING_VALUE_COLUMNS, "KING kinship column")
    retained = 0
    for row in dictionary_rows(path):
        kinship = parse_float(row[kinship_col], path, kinship_col)
        if kinship < cutoff:
            continue
        pair = canonical_pair(row[sample_1_col].strip(), row[sample_2_col].strip())
        validate_pair_samples(pair, cohort, path)
        item = evidence.setdefault(pair, PairEvidence(*pair))
        item.king_kinship = max(
            kinship,
            item.king_kinship if item.king_kinship is not None else -math.inf,
        )
        retained += 1
    return retained


def add_ibd_pairs(
    path: Path,
    cutoff: float,
    cohort: set[str],
    evidence: dict[tuple[str, str], PairEvidence],
) -> int:
    columns = table_columns(path)
    sample_1_col = resolve_column(columns, SAMPLE_1_COLUMNS, "IBD sample-1 column")
    sample_2_col = resolve_column(columns, SAMPLE_2_COLUMNS, "IBD sample-2 column")
    total_col = resolve_column(columns, IBD_TOTAL_COLUMNS, "total IBD length column")
    retained = 0
    for row in dictionary_rows(path):
        total_ibd = parse_float(row[total_col], path, total_col)
        # The paper's rule is strictly greater than 1000 cM.
        if total_ibd <= cutoff:
            continue
        pair = canonical_pair(row[sample_1_col].strip(), row[sample_2_col].strip())
        validate_pair_samples(pair, cohort, path)
        item = evidence.setdefault(pair, PairEvidence(*pair))
        item.total_ibd_cm = max(
            total_ibd,
            item.total_ibd_cm if item.total_ibd_cm is not None else -math.inf,
        )
        retained += 1
    return retained


def validate_pair_samples(pair: tuple[str, str], cohort: set[str], path: Path) -> None:
    absent = [sample for sample in pair if sample not in cohort]
    if absent:
        raise ValueError(
            f"Sample ID(s) in {path} are absent from the cohort list: {', '.join(absent)}"
        )


def adjacency_from_pairs(
    pairs: Iterable[tuple[str, str]],
) -> dict[str, set[str]]:
    adjacency: dict[str, set[str]] = defaultdict(set)
    for sample_1, sample_2 in pairs:
        adjacency[sample_1].add(sample_2)
        adjacency[sample_2].add(sample_1)
    return dict(adjacency)


def connected_components(adjacency: dict[str, set[str]]) -> list[set[str]]:
    remaining = set(adjacency)
    components: list[set[str]] = []
    while remaining:
        start = min(remaining)
        component: set[str] = set()
        queue = deque([start])
        remaining.remove(start)
        while queue:
            node = queue.popleft()
            component.add(node)
            for neighbor in sorted(adjacency[node]):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    queue.append(neighbor)
        components.append(component)
    components.sort(key=lambda component: min(component))
    return components


def greedy_vertex_cover(component: set[str], adjacency: dict[str, set[str]]) -> set[str]:
    active = set(component)
    removed: set[str] = set()
    while True:
        degrees = {
            node: len(adjacency[node] & active)
            for node in active
        }
        maximum = max(degrees.values(), default=0)
        if maximum == 0:
            break
        # Remove the highest-degree vertex; use sample ID for reproducibility.
        selected = min(node for node, degree in degrees.items() if degree == maximum)
        removed.add(selected)
        active.remove(selected)
    return removed


def exact_minimum_vertex_cover(
    component: set[str],
    adjacency: dict[str, set[str]],
    seconds: float,
) -> set[str]:
    deadline = time.monotonic() + seconds
    best = greedy_vertex_cover(component, adjacency)

    def search(active: frozenset[str], removed: frozenset[str]) -> None:
        nonlocal best
        if time.monotonic() > deadline:
            raise SearchTimeout
        if len(removed) >= len(best):
            return

        degrees = {
            node: len(adjacency[node] & active)
            for node in active
        }
        maximum = max(degrees.values(), default=0)
        if maximum == 0:
            best = set(removed)
            return

        # A maximal matching supplies a valid lower bound on the additional
        # number of removals required by any vertex cover.
        unmatched = set(active)
        matching_size = 0
        for node in sorted(active):
            if node not in unmatched:
                continue
            candidates = adjacency[node] & unmatched
            if candidates:
                neighbor = min(candidates)
                unmatched.remove(node)
                unmatched.remove(neighbor)
                matching_size += 1
        if len(removed) + matching_size >= len(best):
            return

        vertex = min(node for node, degree in degrees.items() if degree == maximum)
        neighbors = adjacency[vertex] & active

        # Every cover either contains vertex, or contains all its neighbors.
        search(active - {vertex}, removed | {vertex})
        if len(removed) + len(neighbors) < len(best):
            search(active - neighbors, removed | frozenset(neighbors))

    search(frozenset(component), frozenset())
    return best


def select_removed_samples(
    components: list[set[str]],
    adjacency: dict[str, set[str]],
    mode: str,
    exact_size_limit: int,
    exact_seconds: float,
) -> tuple[set[str], dict[int, str]]:
    removed: set[str] = set()
    algorithms: dict[int, str] = {}
    for component_number, component in enumerate(components, start=1):
        if mode == "remove_all":
            component_removed = set(component)
            algorithm = "remove_all"
        elif len(component) > exact_size_limit:
            component_removed = greedy_vertex_cover(component, adjacency)
            algorithm = "greedy_max_degree_size_fallback"
        else:
            try:
                component_removed = exact_minimum_vertex_cover(
                    component, adjacency, exact_seconds
                )
                algorithm = "exact_minimum_vertex_cover"
            except SearchTimeout:
                component_removed = greedy_vertex_cover(component, adjacency)
                algorithm = "greedy_max_degree_time_fallback"
        removed.update(component_removed)
        algorithms[component_number] = algorithm
    return removed, algorithms


def fmt_float(value: float | None) -> str:
    return "NA" if value is None else format(value, ".12g")


def write_outputs(
    prefix: Path,
    samples: list[str],
    evidence: dict[tuple[str, str], PairEvidence],
    adjacency: dict[str, set[str]],
    components: list[set[str]],
    removed: set[str],
    algorithms: dict[int, str],
    args: argparse.Namespace,
    king_rows: int,
    ibd_rows: int,
) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    unrelated_path = Path(f"{prefix}.unrelated.samples.txt")
    removed_path = Path(f"{prefix}.related_samples_removed.txt")
    pairs_path = Path(f"{prefix}.related_pairs.tsv.gz")
    components_path = Path(f"{prefix}.related_components.tsv.gz")
    summary_path = Path(f"{prefix}.selection_summary.tsv")

    with unrelated_path.open("wt", encoding="utf-8") as handle:
        for sample in samples:
            if sample not in removed:
                handle.write(f"{sample}\n")

    with removed_path.open("wt", encoding="utf-8") as handle:
        for sample in samples:
            if sample in removed:
                handle.write(f"{sample}\n")

    with open_reproducible_gzip(pairs_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("sample_1", "sample_2", "king_kinship", "total_IBD_length_cM", "source")
        )
        for pair in sorted(evidence):
            item = evidence[pair]
            writer.writerow(
                (
                    item.sample_1,
                    item.sample_2,
                    fmt_float(item.king_kinship),
                    fmt_float(item.total_ibd_cm),
                    item.source,
                )
            )

    component_lookup: dict[str, int] = {}
    for number, component in enumerate(components, start=1):
        for sample in component:
            component_lookup[sample] = number

    with open_reproducible_gzip(components_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "component",
                "sample",
                "component_size",
                "degree",
                "action",
                "selection_algorithm",
            )
        )
        for sample in sorted(component_lookup, key=lambda item: (component_lookup[item], item)):
            number = component_lookup[sample]
            writer.writerow(
                (
                    f"R{number}",
                    sample,
                    len(components[number - 1]),
                    len(adjacency[sample]),
                    "remove" if sample in removed else "retain",
                    algorithms[number],
                )
            )

    source_counts = defaultdict(int)
    for item in evidence.values():
        source_counts[item.source] += 1
    fallback_components = sum("fallback" in method for method in algorithms.values())

    with summary_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        metrics = (
            ("removal_mode", args.removal_mode),
            ("king_cutoff", args.king_cutoff),
            ("ibd_total_cm_cutoff", args.ibd_total_cm_cutoff),
            ("cohort_samples", len(samples)),
            ("related_pairs_union", len(evidence)),
            ("king_input_pairs_passing_cutoff", king_rows),
            ("ibd_input_pairs_passing_cutoff", ibd_rows),
            ("king_only_pairs", source_counts["KING"]),
            ("ibd_only_pairs", source_counts["IBD"]),
            ("king_and_ibd_pairs", source_counts["KING+IBD"]),
            ("related_components", len(components)),
            ("samples_in_related_components", len(adjacency)),
            ("samples_removed", len(removed)),
            ("samples_retained", len(samples) - len(removed)),
            ("fallback_components", fallback_components),
            ("exact_component_size_limit", args.exact_component_size_limit),
            ("exact_search_seconds_per_component", args.exact_search_seconds),
        )
        writer.writerows(metrics)


def main() -> int:
    args = parse_arguments()
    try:
        samples = read_samples(args.samples)
        cohort = set(samples)
        evidence: dict[tuple[str, str], PairEvidence] = {}
        king_rows = add_king_pairs(args.king, args.king_cutoff, cohort, evidence)
        ibd_rows = add_ibd_pairs(
            args.ibd, args.ibd_total_cm_cutoff, cohort, evidence
        )
        adjacency = adjacency_from_pairs(evidence)
        components = connected_components(adjacency)
        removed, algorithms = select_removed_samples(
            components,
            adjacency,
            args.removal_mode,
            args.exact_component_size_limit,
            args.exact_search_seconds,
        )
        write_outputs(
            args.output_prefix,
            samples,
            evidence,
            adjacency,
            components,
            removed,
            algorithms,
            args,
            king_rows,
            ibd_rows,
        )
    except (OSError, ValueError, csv.Error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
