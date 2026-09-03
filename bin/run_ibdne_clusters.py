#!/usr/bin/env python3
"""Prepare within-cluster Hap-IBD segments and run IBDNe per final cluster."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shlex
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def natural_key(label: str):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", label)
    ]


def parse_boolean(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise argparse.ArgumentTypeError("Expected true or false.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segments", required=True, nargs="+", type=Path)
    parser.add_argument("--maps", required=True, nargs="+", type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--ibdne-jar", required=True, type=Path)
    parser.add_argument("--method", required=True)
    parser.add_argument("--minimum-cluster-size", required=True, type=int)
    parser.add_argument("--mincm", required=True, type=float)
    parser.add_argument("--nits", required=True, type=int)
    parser.add_argument("--nboots", required=True, type=int)
    parser.add_argument("--filtersamples", required=True, type=parse_boolean)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--parallel-clusters", required=True, type=int)
    parser.add_argument("--java-heap-gb", required=True, type=int)
    parser.add_argument("--temporary-directory", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    return parser.parse_args()


@dataclass(frozen=True)
class ClusterRun:
    cluster: str
    sample_count: int
    segment_count: int
    npairs: int
    input_path: Path


def read_membership(path: Path):
    membership: dict[str, str] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Membership file has no header.")
        id_column = "ID" if "ID" in reader.fieldnames else reader.fieldnames[0]
        if "final_cluster" not in reader.fieldnames:
            raise ValueError("Membership file lacks final_cluster.")
        for row in reader:
            sample_id = row[id_column]
            cluster = row["final_cluster"]
            if not sample_id or not cluster:
                raise ValueError("Empty sample ID or final cluster label.")
            if sample_id in membership:
                raise ValueError(f"Duplicate membership ID: {sample_id}")
            membership[sample_id] = cluster
    if not membership:
        raise ValueError("Membership file contains no samples.")
    return membership, Counter(membership.values())


def write_combined_map(map_paths: list[Path], output: Path) -> None:
    """Concatenate chromosome-specific PLINK maps in supplied chromosome order."""
    with output.open("wt", encoding="utf-8", newline="") as destination:
        for map_path in map_paths:
            with open_text(map_path) as source:
                for line in source:
                    if line.strip():
                        destination.write(line if line.endswith("\n") else line + "\n")
    if output.stat().st_size == 0:
        raise ValueError("Combined genetic map is empty.")


def prepare_cluster_inputs(
    segment_paths: list[Path],
    membership: dict[str, str],
    sizes: Counter,
    minimum_size: int,
    work_directory: Path,
):
    eligible = {
        cluster for cluster, size in sizes.items() if size >= minimum_size
    }
    segment_counts = Counter()
    unknown_ids = set()
    input_paths = {
        cluster: work_directory / f"{cluster}.ibd"
        for cluster in eligible
    }

    with ExitStack() as stack:
        handles = {
            cluster: stack.enter_context(
                path.open("wt", encoding="utf-8", newline="")
            )
            for cluster, path in input_paths.items()
        }
        for segment_path in segment_paths:
            with open_text(segment_path) as source:
                for line_number, line in enumerate(source, 1):
                    if not line.strip():
                        continue
                    fields = line.split()
                    if len(fields) < 8:
                        raise ValueError(
                            f"{segment_path}:{line_number} has fewer than 8 fields."
                        )
                    first, second = fields[0], fields[2]
                    if first not in membership or second not in membership:
                        if len(unknown_ids) < 20:
                            unknown_ids.update(
                                sample for sample in (first, second)
                                if sample not in membership
                            )
                        continue
                    cluster = membership[first]
                    if cluster == membership[second] and cluster in eligible:
                        handles[cluster].write(line)
                        if not line.endswith("\n"):
                            handles[cluster].write("\n")
                        segment_counts[cluster] += 1

    if unknown_ids:
        examples = ", ".join(sorted(unknown_ids)[:20])
        raise ValueError(
            "Hap-IBD segment IDs were absent from final membership. "
            f"Examples: {examples}"
        )

    runs = []
    for cluster in sorted(eligible, key=natural_key):
        sample_count = sizes[cluster]
        haplotypes = 2 * sample_count
        npairs = haplotypes * (haplotypes - 2) // 2
        runs.append(ClusterRun(
            cluster=cluster,
            sample_count=sample_count,
            segment_count=segment_counts[cluster],
            npairs=npairs,
            input_path=input_paths[cluster],
        ))
    return runs


def run_ibdne(
    run: ClusterRun,
    args: argparse.Namespace,
    combined_map: Path,
    threads_per_run: int,
):
    cluster_directory = args.output_directory / run.cluster
    cluster_directory.mkdir(parents=True, exist_ok=True)
    prefix = cluster_directory / run.cluster
    command = [
        "java",
        f"-Xmx{args.java_heap_gb}g",
        "-jar",
        str(args.ibdne_jar),
        f"map={combined_map}",
        f"out={prefix}",
        f"nthreads={threads_per_run}",
        f"filtersamples={str(args.filtersamples).lower()}",
        f"npairs={run.npairs}",
        f"nits={args.nits}",
        f"nboots={args.nboots}",
        f"mincm={args.mincm}",
        f"seed={args.seed}",
    ]
    command_path = cluster_directory / f"{run.cluster}.command.txt"
    command_path.write_text(shlex.join(command) + "\n", encoding="utf-8")
    with run.input_path.open("rb") as standard_input:
        completed = subprocess.run(
            command,
            stdin=standard_input,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    driver_log = cluster_directory / f"{run.cluster}.driver.log"
    driver_log.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(
            completed.returncode, command, output=completed.stdout
        )
    # IBDNe may omit .pair.excl when filtersamples=false and may omit
    # .region.excl when no region is excluded. Only the three core outputs
    # are required to mark the cluster analysis as successful.
    missing_outputs = [
        Path(f"{prefix}{suffix}")
        for suffix in (".ne", ".boot", ".log")
        if not Path(f"{prefix}{suffix}").is_file()
    ]
    if missing_outputs:
        missing = ", ".join(path.name for path in missing_outputs)
        raise FileNotFoundError(f"IBDNe did not create expected outputs: {missing}")
    return run


def write_status(path: Path, method: str, sizes: Counter, minimum_size: int,
                 runs: list[ClusterRun], completed: set[str], failures: dict[str, str]):
    run_by_cluster = {run.cluster: run for run in runs}
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "method", "final_cluster", "sample_count", "minimum_cluster_size",
            "within_cluster_segments", "npairs", "status", "message",
        ])
        for cluster in sorted(sizes, key=natural_key):
            size = sizes[cluster]
            if size < minimum_size:
                row = [method, cluster, size, minimum_size, 0, "NA",
                       "skipped_small_cluster", ""]
            else:
                run = run_by_cluster[cluster]
                if run.segment_count == 0:
                    status, message = "skipped_no_segments", ""
                elif cluster in completed:
                    status, message = "completed", ""
                else:
                    status = "failed"
                    message = failures.get(cluster, "Unknown IBDNe failure")
                row = [method, cluster, size, minimum_size, run.segment_count,
                       run.npairs, status, message]
            writer.writerow(row)


def validate_args(args: argparse.Namespace) -> None:
    for path in [*args.segments, *args.maps, args.membership, args.ibdne_jar]:
        if not path.is_file():
            raise FileNotFoundError(path)
    if shutil.which("java") is None:
        raise FileNotFoundError("java")
    if args.minimum_cluster_size < 2:
        raise ValueError("--minimum-cluster-size must be at least 2.")
    if args.mincm <= 0 or args.nits < 1 or args.nboots < 1:
        raise ValueError("mincm, nits and nboots must be positive.")
    if args.threads < 1 or args.parallel_clusters < 1 or args.java_heap_gb < 1:
        raise ValueError("threads, parallel-clusters and java-heap-gb must be positive.")


def main() -> None:
    args = parse_args()
    validate_args(args)
    args.temporary_directory.mkdir(parents=True, exist_ok=True)
    args.output_directory.mkdir(parents=True, exist_ok=True)

    membership, sizes = read_membership(args.membership)
    combined_map = args.temporary_directory / "ibdne.genomewide.map"
    write_combined_map(args.maps, combined_map)
    runs = prepare_cluster_inputs(
        args.segments,
        membership,
        sizes,
        args.minimum_cluster_size,
        args.temporary_directory,
    )

    runnable = [run for run in runs if run.segment_count > 0]
    maximum_parallel = min(args.parallel_clusters, len(runnable)) if runnable else 1
    threads_per_run = max(1, args.threads // maximum_parallel)
    completed: set[str] = set()
    failures: dict[str, str] = {}

    with ThreadPoolExecutor(max_workers=maximum_parallel) as executor:
        futures = {
            executor.submit(run_ibdne, run, args, combined_map, threads_per_run): run
            for run in runnable
        }
        for future in as_completed(futures):
            run = futures[future]
            try:
                future.result()
                completed.add(run.cluster)
            except Exception as error:  # Preserve every cluster's status before failing.
                failures[run.cluster] = str(error)

    status_path = args.output_directory / f"{args.method}.ibdne_status.tsv"
    write_status(
        status_path, args.method, sizes, args.minimum_cluster_size,
        runs, completed, failures,
    )
    if failures:
        failed = ", ".join(sorted(failures, key=natural_key))
        raise RuntimeError(f"IBDNe failed for: {failed}")


if __name__ == "__main__":
    main()