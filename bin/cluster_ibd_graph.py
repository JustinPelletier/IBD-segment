#!/usr/bin/env python3

"""
Recursively cluster a weighted IBD-sharing graph using Louvain or Leiden.

The graph is constructed from:

1. A genome-wide pair-level IBD summary containing graph edges.
2. A complete sample list containing all graph vertices.

Every sample is added as a vertex, including participants with no detected
IBD-sharing edges. Such participants remain in the graph as isolated nodes.

Fixed-depth recursive clustering
--------------------------------
Level 0 assigns the complete cohort to one initial community.

At level 1, the selected clustering algorithm is applied to the complete graph.

At each subsequent level, every sufficiently large community from the previous
level is independently submitted to the same clustering algorithm. Communities
smaller than --min-cluster-size remain in the output but are not refined again.

Within every proposed split, child communities smaller than
--min-cluster-size are pooled into one residual child. When the pooled residual
still does not meet the minimum, it is absorbed into the retained large child
with which it shares the greatest total edge weight. The split is rejected only
when fewer than two valid children can remain.

Exactly --max-levels levels are written. When no eligible parent community can
be split at a level, its labels are carried forward unchanged. Global weighted
modularity is reported as a diagnostic but never controls refinement or selects
the final level. The deepest requested level is always the terminal membership
used by downstream FST clumping.
"""

import argparse
import csv
import gzip
import inspect
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import IO, Iterable

import igraph as ig


def parse_boolean(value: str) -> bool:
    """Convert a command-line string to a Boolean value."""

    normalized_value = value.strip().lower()

    if normalized_value == "true":
        return True

    if normalized_value == "false":
        return False

    raise argparse.ArgumentTypeError(
        "Expected either 'true' or 'false'."
    )


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Recursively cluster a weighted IBD-sharing graph using "
            "Louvain or Leiden."
        )
    )

    parser.add_argument(
        "--edges",
        required=True,
        type=Path,
        help=(
            "Genome-wide pair-level IBD summary in .tsv or .tsv.gz format."
        ),
    )

    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
        help=(
            "Complete sample list with one sample ID per line."
        ),
    )

    parser.add_argument(
        "--method",
        required=True,
        choices=("louvain", "leiden"),
        help="Community-detection algorithm.",
    )

    parser.add_argument(
        "--max-levels",
        required=True,
        type=int,
        help=(
            "Maximum recursive refinement depth. Level 1 clusters the "
            "complete graph; subsequent levels recluster communities."
        ),
    )

    parser.add_argument(
        "--weight-column",
        required=True,
        help=(
            "Column from the edge file to use as the graph edge weight."
        ),
    )

    parser.add_argument(
        "--min-cluster-size",
        required=True,
        type=int,
        help=(
            "Communities smaller than this value remain in the results "
            "but are not submitted to another refinement level."
        ),
    )

    parser.add_argument(
        "--min-modularity-gain",
        required=False,
        type=float,
        default=0.0,
        help=(
            "Deprecated compatibility option. Modularity gain is reported "
            "but does not control fixed-depth refinement."
        ),
    )

    parser.add_argument(
        "--resolution",
        required=True,
        type=float,
        help=(
            "Resolution parameter used by Leiden. The standard Louvain "
            "implementation in python-igraph does not use this parameter."
        ),
    )

    parser.add_argument(
        "--seed",
        required=True,
        type=int,
        help="Random seed used for reproducible clustering.",
    )

    parser.add_argument(
        "--auto-select",
        required=False,
        type=parse_boolean,
        default=False,
        help=(
            "Deprecated compatibility option. The deepest requested level "
            "is always used."
        ),
    )

    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix used for all output files.",
    )

    return parser.parse_args()


def open_text_file(
    path: Path,
    mode: str = "rt",
) -> IO[str]:
    """Open a plain-text or gzip-compressed file."""

    if str(path).endswith(".gz"):
        return gzip.open(
            path,
            mode=mode,
            encoding="utf-8",
            newline="",
        )

    return path.open(
        mode=mode,
        encoding="utf-8",
        newline="",
    )


def validate_arguments(args: argparse.Namespace) -> None:
    """Validate file paths and clustering parameters."""

    if not args.edges.is_file():
        raise SystemExit(
            f"ERROR: edge file does not exist: {args.edges}"
        )

    if not args.samples.is_file():
        raise SystemExit(
            f"ERROR: sample file does not exist: {args.samples}"
        )

    if args.max_levels < 1:
        raise SystemExit(
            "ERROR: --max-levels must be greater than or equal to one."
        )

    if args.min_cluster_size < 2:
        raise SystemExit(
            "ERROR: --min-cluster-size must be greater than or equal to two."
        )

    if args.min_modularity_gain < 0:
        raise SystemExit(
            "ERROR: --min-modularity-gain must be greater than or equal to zero."
        )

    if args.resolution <= 0:
        raise SystemExit(
            "ERROR: --resolution must be greater than zero."
        )


def read_samples(sample_path: Path) -> list[str]:
    """
    Read the complete cohort sample list.

    Duplicate sample IDs are not allowed because vertex names must be unique.
    """

    samples: list[str] = []
    observed_samples: set[str] = set()

    with open_text_file(sample_path) as sample_handle:
        for line_number, line in enumerate(
            sample_handle,
            start=1,
        ):
            sample_id = line.strip()

            if not sample_id:
                continue

            if sample_id in observed_samples:
                raise SystemExit(
                    "ERROR: duplicate sample ID in the sample list.\n"
                    f"File: {sample_path}\n"
                    f"Line: {line_number}\n"
                    f"Sample: {sample_id}"
                )

            observed_samples.add(sample_id)
            samples.append(sample_id)

    if not samples:
        raise SystemExit(
            f"ERROR: the sample file is empty: {sample_path}"
        )

    return samples


def read_graph(
    edge_path: Path,
    samples: list[str],
    weight_column: str,
) -> ig.Graph:
    """
    Construct an undirected weighted graph.

    Every sample is added as a vertex, even when the sample has no edge.
    The pair-summary file is expected to contain one row per sample pair.
    """

    sample_to_index = {
        sample_id: index
        for index, sample_id in enumerate(samples)
    }

    graph_edges: list[tuple[int, int]] = []
    graph_weights: list[float] = []

    with open_text_file(edge_path) as edge_handle:
        reader = csv.DictReader(
            edge_handle,
            delimiter="\t",
        )

        required_columns = {
            "ID1",
            "ID2",
            weight_column,
        }

        observed_columns = set(
            reader.fieldnames or []
        )

        missing_columns = (
            required_columns - observed_columns
        )

        if missing_columns:
            raise SystemExit(
                "ERROR: the edge file is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        for line_number, row in enumerate(
            reader,
            start=2,
        ):
            sample1 = row["ID1"]
            sample2 = row["ID2"]

            if sample1 not in sample_to_index:
                raise SystemExit(
                    "ERROR: an edge contains a sample absent from "
                    "the complete sample list.\n"
                    f"File: {edge_path}\n"
                    f"Line: {line_number}\n"
                    f"Sample: {sample1}"
                )

            if sample2 not in sample_to_index:
                raise SystemExit(
                    "ERROR: an edge contains a sample absent from "
                    "the complete sample list.\n"
                    f"File: {edge_path}\n"
                    f"Line: {line_number}\n"
                    f"Sample: {sample2}"
                )

            if sample1 == sample2:
                raise SystemExit(
                    "ERROR: a self-loop was found in the edge file.\n"
                    f"File: {edge_path}\n"
                    f"Line: {line_number}\n"
                    f"Sample: {sample1}"
                )

            try:
                edge_weight = float(
                    row[weight_column]
                )
            except ValueError as error:
                raise SystemExit(
                    "ERROR: invalid edge weight.\n"
                    f"File: {edge_path}\n"
                    f"Line: {line_number}\n"
                    f"Column: {weight_column}\n"
                    f"Value: {row[weight_column]!r}"
                ) from error

            if (
                not math.isfinite(edge_weight)
                or edge_weight < 0
            ):
                raise SystemExit(
                    "ERROR: edge weights must be finite and non-negative.\n"
                    f"File: {edge_path}\n"
                    f"Line: {line_number}\n"
                    f"Value: {edge_weight}"
                )

            # A zero-weight pair is equivalent to no graph edge.
            if edge_weight == 0:
                continue

            graph_edges.append(
                (
                    sample_to_index[sample1],
                    sample_to_index[sample2],
                )
            )

            graph_weights.append(
                edge_weight
            )

    graph = ig.Graph(
        n=len(samples),
        edges=graph_edges,
        directed=False,
    )

    graph.vs["name"] = samples
    graph.es["weight"] = graph_weights

    if graph.ecount() == 0:
        raise SystemExit(
            "ERROR: no positive-weight IBD-sharing edges were found."
        )

    if not graph.is_simple():
        raise SystemExit(
            "ERROR: the pair-summary file contains duplicate sample pairs "
            "or self-loops. It must contain one row per unique sample pair."
        )

    return graph


def run_louvain(
    graph: ig.Graph,
) -> list[int]:
    """Run weighted Louvain clustering."""

    clustering = graph.community_multilevel(
        weights="weight"
    )

    return clustering.membership


def run_leiden(
    graph: ig.Graph,
    resolution: float,
) -> list[int]:
    """
    Run weighted Leiden clustering.

    The resolution argument name changed between python-igraph releases.
    Both current and older supported APIs are handled here.
    """

    base_arguments = {
        "objective_function": "modularity",
        "weights": "weight",
        "n_iterations": -1,
    }

    try:
        signature = inspect.signature(
            graph.community_leiden
        )

        if "resolution" in signature.parameters:
            resolution_arguments = {
                "resolution": resolution
            }
        else:
            resolution_arguments = {
                "resolution_parameter": resolution
            }

        clustering = graph.community_leiden(
            **base_arguments,
            **resolution_arguments,
        )

    except (TypeError, ValueError):
        # Some compiled igraph methods do not expose a complete signature.
        try:
            clustering = graph.community_leiden(
                resolution=resolution,
                **base_arguments,
            )
        except TypeError:
            clustering = graph.community_leiden(
                resolution_parameter=resolution,
                **base_arguments,
            )

    return clustering.membership


def partition_graph(
    graph: ig.Graph,
    method: str,
    resolution: float,
) -> list[int]:
    """Apply the requested community-detection algorithm."""

    if graph.vcount() == 0:
        return []

    if graph.ecount() == 0:
        return [0] * graph.vcount()

    if method == "louvain":
        return run_louvain(
            graph
        )

    return run_leiden(
        graph,
        resolution,
    )


def calculate_global_modularity(
    graph: ig.Graph,
    labels: list[str],
) -> float:
    """Calculate weighted modularity for a complete graph partition."""

    label_to_integer = {
        label: index
        for index, label in enumerate(
            dict.fromkeys(labels)
        )
    }

    numeric_membership = [
        label_to_integer[label]
        for label in labels
    ]

    return graph.modularity(
        membership=numeric_membership,
        weights="weight",
    )


def group_vertices_by_community(
    membership: list[str],
) -> dict[str, list[int]]:
    """Group graph vertex indices according to community membership."""

    grouped_vertices: dict[str, list[int]] = defaultdict(list)

    for vertex_index, community in enumerate(
        membership
    ):
        grouped_vertices[community].append(
            vertex_index
        )

    return dict(grouped_vertices)


def create_stable_local_labels(
    subgraph: ig.Graph,
    local_membership: list[int],
    residual_community: int | None = None,
) -> dict[int, int]:
    """
    Create deterministic local cluster labels.

    Communities are ordered by their lexicographically smallest sample ID.
    A pooled residual community, when present, is placed last. This avoids
    relying on arbitrary internal community numbers while keeping the residual
    child easy to identify.
    """

    local_communities: dict[int, list[str]] = defaultdict(list)

    for vertex_index, community in enumerate(
        local_membership
    ):
        local_communities[community].append(
            subgraph.vs[vertex_index]["name"]
        )

    ordered_communities = sorted(
        (
            community
            for community in local_communities
            if community != residual_community
        ),
        key=lambda community: min(
            local_communities[community]
        ),
    )

    if residual_community is not None:
        ordered_communities.append(
            residual_community
        )

    return {
        community: index + 1
        for index, community in enumerate(
            ordered_communities
        )
    }


def consolidate_small_communities(
    subgraph: ig.Graph,
    local_membership: list[int],
    min_cluster_size: int,
) -> tuple[list[int], int | None, int, bool, bool]:
    """
    Pool undersized child communities into one residual child.

    Returns
    -------
    consolidated_membership
        Local membership after pooling undersized children. The original
        membership is returned when the proposed split must be rejected.

    residual_community
        Identifier of the pooled residual child, or None when no pooling was
        required or when the split was rejected.

    small_communities_merged
        Number of undersized child communities included in the residual.

    reject_split
        True when fewer than two children meeting the minimum can remain.

    undersized_residual_absorbed
        True when a pooled residual smaller than the minimum was assigned to
        the most strongly connected retained large child.
    """

    community_sizes = Counter(
        local_membership
    )

    small_communities = {
        community
        for community, size in community_sizes.items()
        if size < min_cluster_size
    }

    if not small_communities:
        return local_membership, None, 0, False, False

    residual_size = sum(
        community_sizes[community]
        for community in small_communities
    )

    if residual_size < min_cluster_size:
        large_communities = {
            community
            for community, size in community_sizes.items()
            if size >= min_cluster_size
        }

        # Absorbing the residual when only one large child exists would leave
        # the parent unsplit. In that case, retain the original parent label.
        if len(large_communities) < 2:
            return (
                local_membership,
                None,
                len(small_communities),
                True,
                False,
            )

        connection_weights = {
            community: 0.0
            for community in large_communities
        }

        for edge_index, (vertex1, vertex2) in enumerate(
            subgraph.get_edgelist()
        ):
            community1 = local_membership[vertex1]
            community2 = local_membership[vertex2]
            edge_weight = subgraph.es[edge_index]["weight"]

            if (
                community1 in small_communities
                and community2 in large_communities
            ):
                connection_weights[community2] += edge_weight

            elif (
                community2 in small_communities
                and community1 in large_communities
            ):
                connection_weights[community1] += edge_weight

        community_minimum_sample = {}

        for community in large_communities:
            community_minimum_sample[community] = min(
                subgraph.vs[vertex_index]["name"]
                for vertex_index, observed_community in enumerate(
                    local_membership
                )
                if observed_community == community
            )

        absorption_target = sorted(
            large_communities,
            key=lambda community: (
                -connection_weights[community],
                -community_sizes[community],
                community_minimum_sample[community],
            ),
        )[0]

        consolidated_membership = [
            (
                absorption_target
                if community in small_communities
                else community
            )
            for community in local_membership
        ]

        return (
            consolidated_membership,
            None,
            len(small_communities),
            False,
            True,
        )

    residual_community = max(
        community_sizes
    ) + 1

    consolidated_membership = [
        (
            residual_community
            if community in small_communities
            else community
        )
        for community in local_membership
    ]

    if len(set(consolidated_membership)) < 2:
        return (
            local_membership,
            None,
            len(small_communities),
            True,
            False,
        )

    return (
        consolidated_membership,
        residual_community,
        len(small_communities),
        False,
        False,
    )


def refine_membership(
    graph: ig.Graph,
    previous_membership: list[str],
    method: str,
    resolution: float,
    min_cluster_size: int,
) -> tuple[list[str], int, int, int, int, int, int]:
    """
    Recursively cluster communities from the previous refinement level.

    Returns
    -------
    new_membership
        Hierarchical community labels for the new level.

    communities_split
        Number of parent communities divided into at least two communities.

    communities_eligible
        Number of parent communities meeting the minimum size.

    communities_carried_forward
        Number of parent communities whose labels were unchanged.

    small_child_communities_merged
        Number of undersized child communities pooled across accepted splits.

    splits_rejected_small_residual
        Number of proposed splits rejected because pooling could not produce
        at least two children meeting the minimum size.

    undersized_residuals_absorbed
        Number of pooled residuals below the minimum absorbed into their most
        strongly connected retained large child.
    """

    new_membership = previous_membership.copy()
    communities_split = 0
    communities_eligible = 0
    small_child_communities_merged = 0
    splits_rejected_small_residual = 0
    undersized_residuals_absorbed = 0

    grouped_vertices = group_vertices_by_community(
        previous_membership
    )

    graph_sample_to_index = {
        sample_id: vertex_index
        for vertex_index, sample_id in enumerate(
            graph.vs["name"]
        )
    }

    for parent_label, vertex_indices in grouped_vertices.items():
        # Small communities remain in the output but are not refined again.
        if len(vertex_indices) < min_cluster_size:
            continue

        communities_eligible += 1

        subgraph = graph.induced_subgraph(
            vertex_indices
        )

        if subgraph.ecount() == 0:
            continue

        local_membership = partition_graph(
            graph=subgraph,
            method=method,
            resolution=resolution,
        )

        observed_local_communities = set(
            local_membership
        )

        if len(observed_local_communities) < 2:
            continue

        (
            local_membership,
            residual_community,
            small_communities_merged,
            reject_split,
            undersized_residual_absorbed,
        ) = consolidate_small_communities(
            subgraph=subgraph,
            local_membership=local_membership,
            min_cluster_size=min_cluster_size,
        )

        if reject_split:
            splits_rejected_small_residual += 1
            continue

        small_child_communities_merged += (
            small_communities_merged
        )

        undersized_residuals_absorbed += int(
            undersized_residual_absorbed
        )

        stable_labels = create_stable_local_labels(
            subgraph,
            local_membership,
            residual_community,
        )

        for local_vertex_index, local_community in enumerate(
            local_membership
        ):
            sample_id = subgraph.vs[
                local_vertex_index
            ]["name"]

            original_vertex_index = graph_sample_to_index[
                sample_id
            ]

            new_membership[original_vertex_index] = (
                f"{parent_label}."
                f"{stable_labels[local_community]}"
            )

        communities_split += 1

    communities_carried_forward = (
        len(grouped_vertices) - communities_split
    )

    return (
        new_membership,
        communities_split,
        communities_eligible,
        communities_carried_forward,
        small_child_communities_merged,
        splits_rejected_small_residual,
        undersized_residuals_absorbed,
    )


def write_tsv(
    output_path: Path,
    header: list[str],
    rows: Iterable[Iterable],
    compressed: bool = False,
) -> None:
    """Write a tab-separated output file."""

    if compressed:
        output_handle = gzip.open(
            output_path,
            mode="wt",
            encoding="utf-8",
            newline="",
        )
    else:
        output_handle = output_path.open(
            mode="wt",
            encoding="utf-8",
            newline="",
        )

    with output_handle:
        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(header)
        writer.writerows(rows)


def write_outputs(
    graph: ig.Graph,
    method: str,
    output_prefix: str,
    memberships: dict[int, list[str]],
    diagnostics: list[dict],
    selected_level: int,
    min_cluster_size: int,
    max_levels: int,
) -> None:
    """Write clustering memberships and diagnostic outputs."""

    sample_ids = graph.vs["name"]
    levels = sorted(memberships)

    membership_path = Path(
        f"{output_prefix}.membership.tsv.gz"
    )

    selected_membership_path = Path(
        f"{output_prefix}.selected_membership.tsv.gz"
    )

    diagnostics_path = Path(
        f"{output_prefix}.diagnostics.tsv"
    )

    cluster_sizes_path = Path(
        f"{output_prefix}.cluster_sizes.tsv"
    )

    membership_header = [
        "ID",
        *[
            f"level_{level}"
            for level in levels
        ],
    ]

    membership_rows = (
        [
            sample_id,
            *[
                memberships[level][sample_index]
                for level in levels
            ],
        ]
        for sample_index, sample_id in enumerate(
            sample_ids
        )
    )

    write_tsv(
        output_path=membership_path,
        header=membership_header,
        rows=membership_rows,
        compressed=True,
    )

    selected_rows = (
        (
            sample_id,
            memberships[selected_level][sample_index],
        )
        for sample_index, sample_id in enumerate(
            sample_ids
        )
    )

    write_tsv(
        output_path=selected_membership_path,
        header=[
            "ID",
            "cluster",
        ],
        rows=selected_rows,
        compressed=True,
    )

    diagnostic_header = [
        "method",
        "level",
        "number_of_clusters",
        "modularity",
        "modularity_gain",
        "communities_eligible",
        "communities_split",
        "communities_carried_forward",
        "small_child_communities_merged",
        "splits_rejected_small_residual",
        "undersized_residuals_absorbed",
        "selected",
        "stopping_reason",
    ]

    diagnostic_rows = (
        (
            method,
            row["level"],
            row["number_of_clusters"],
            format(row["modularity"], ".10g"),
            (
                ""
                if row["modularity_gain"] is None
                else format(
                    row["modularity_gain"],
                    ".10g",
                )
            ),
            row["communities_eligible"],
            row["communities_split"],
            row["communities_carried_forward"],
            row["small_child_communities_merged"],
            row["splits_rejected_small_residual"],
            row["undersized_residuals_absorbed"],
            row["level"] == selected_level,
            row["stopping_reason"],
        )
        for row in diagnostics
    )

    write_tsv(
        output_path=diagnostics_path,
        header=diagnostic_header,
        rows=diagnostic_rows,
    )

    cluster_size_rows = []

    for level in levels:
        cluster_counts = Counter(
            memberships[level]
        )

        for cluster, cluster_size in sorted(
            cluster_counts.items()
        ):
            cluster_size_rows.append(
                (
                    level,
                    cluster,
                    cluster_size,
                    (
                        level < max_levels
                        and cluster_size >= min_cluster_size
                    ),
                    level == selected_level,
                )
            )

    write_tsv(
        output_path=cluster_sizes_path,
        header=[
            "level",
            "cluster",
            "size",
            "eligible_for_further_refinement",
            "terminal_level",
        ],
        rows=cluster_size_rows,
    )


def main() -> None:
    """Run recursive IBD graph clustering."""

    args = parse_arguments()
    validate_arguments(args)

    random_number_generator = random.Random(
        args.seed
    )

    ig.set_random_number_generator(
        random_number_generator
    )

    print(
        f"Reading complete sample list: {args.samples}",
        file=sys.stderr,
    )

    samples = read_samples(
        args.samples
    )

    print(
        f"Reading weighted IBD graph: {args.edges}",
        file=sys.stderr,
    )

    graph = read_graph(
        edge_path=args.edges,
        samples=samples,
        weight_column=args.weight_column,
    )

    isolated_vertex_count = sum(
        degree == 0
        for degree in graph.degree()
    )

    print(
        (
            f"Graph contains {graph.vcount():,} vertices, "
            f"{graph.ecount():,} edges and "
            f"{isolated_vertex_count:,} isolated vertices."
        ),
        file=sys.stderr,
    )

    memberships: dict[int, list[str]] = {
        0: ["C1"] * graph.vcount()
    }

    initial_modularity = calculate_global_modularity(
        graph,
        memberships[0],
    )

    diagnostics = [
        {
            "level": 0,
            "number_of_clusters": 1,
            "modularity": initial_modularity,
            "modularity_gain": None,
            "communities_eligible": int(
                graph.vcount() >= args.min_cluster_size
            ),
            "communities_split": 0,
            "communities_carried_forward": 0,
            "small_child_communities_merged": 0,
            "splits_rejected_small_residual": 0,
            "undersized_residuals_absorbed": 0,
            "stopping_reason": "",
        }
    ]

    for level in range(
        1,
        args.max_levels + 1,
    ):
        previous_membership = memberships[
            level - 1
        ]

        (
            new_membership,
            communities_split,
            communities_eligible,
            communities_carried_forward,
            small_child_communities_merged,
            splits_rejected_small_residual,
            undersized_residuals_absorbed,
        ) = refine_membership(
            graph=graph,
            previous_membership=previous_membership,
            method=args.method,
            resolution=args.resolution,
            min_cluster_size=args.min_cluster_size,
        )

        memberships[level] = new_membership

        modularity = calculate_global_modularity(
            graph,
            new_membership,
        )

        modularity_gain = (
            modularity
            - diagnostics[-1]["modularity"]
        )

        diagnostics.append(
            {
                "level": level,
                "number_of_clusters": len(
                    set(new_membership)
                ),
                "modularity": modularity,
                "modularity_gain": modularity_gain,
                "communities_eligible": communities_eligible,
                "communities_split": communities_split,
                "communities_carried_forward": communities_carried_forward,
                "small_child_communities_merged": (
                    small_child_communities_merged
                ),
                "splits_rejected_small_residual": (
                    splits_rejected_small_residual
                ),
                "undersized_residuals_absorbed": (
                    undersized_residuals_absorbed
                ),
                "stopping_reason": (
                    "fixed_depth_completed"
                    if level == args.max_levels
                    else ""
                ),
            }
        )

        print(
            (
                f"{args.method.capitalize()} level {level}: "
                f"{len(set(new_membership)):,} clusters; "
                f"modularity={modularity:.6f}; "
                f"gain={modularity_gain:.6f}; "
                f"eligible={communities_eligible:,}; "
                f"split={communities_split:,}; "
                f"carried_forward={communities_carried_forward:,}; "
                f"small_children_merged={small_child_communities_merged:,}; "
                f"splits_rejected_small_residual="
                f"{splits_rejected_small_residual:,}; "
                f"undersized_residuals_absorbed="
                f"{undersized_residuals_absorbed:,}."
            ),
            file=sys.stderr,
        )

    selected_level = args.max_levels

    write_outputs(
        graph=graph,
        method=args.method,
        output_prefix=args.output_prefix,
        memberships=memberships,
        diagnostics=diagnostics,
        selected_level=selected_level,
        min_cluster_size=args.min_cluster_size,
        max_levels=args.max_levels,
    )

    print(
        (
            f"Terminal {args.method.capitalize()} refinement "
            f"level: {selected_level}"
        ),
        file=sys.stderr,
    )

    print(
        "Stopping reason: fixed_depth_completed",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
