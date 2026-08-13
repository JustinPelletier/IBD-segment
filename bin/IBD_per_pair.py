#!/usr/bin/env python3

"""
Generate pair-level IBD-sharing summary statistics from Hap-IBD output.

The script accepts one or more headerless Hap-IBD .ibd or .ibd.gz files.
Each Hap-IBD row must contain the following eight tab-separated fields:

    ID1
    HAPLOTYPE1
    ID2
    HAPLOTYPE2
    CHROMOSOME
    START_BP
    END_BP
    LENGTH_CM

The output contains one row for every pair sharing at least one detected
IBD segment. Pairs without detected IBD segments are not written because
their absence from the edge list implicitly represents zero IBD sharing.

A disk-backed SQLite database is used to avoid loading all IBD segments
or all sample pairs into memory.
"""

import argparse
import csv
import gzip
import math
import sqlite3
import sys
from pathlib import Path
from typing import IO, Iterable


EXPECTED_HAPIBD_COLUMNS = 8
BATCH_SIZE = 100_000


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Generate pair-level IBD-sharing summary statistics from one "
            "or more Hap-IBD segment files."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        nargs="+",
        type=Path,
        help=(
            "One or more headerless Hap-IBD .ibd or .ibd.gz files. "
            "Multiple chromosome files can be supplied for a genome-wide summary."
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output pair-summary file in .tsv or .tsv.gz format.",
    )

    parser.add_argument(
        "--threshold-cm",
        required=True,
        type=float,
        help=(
            "Segment-length threshold in centimorgans. The output reports "
            "the number and total length of segments greater than or equal "
            "to this threshold."
        ),
    )

    return parser.parse_args()


def open_text_file(
    path: Path,
    mode: str = "rt",
) -> IO[str]:
    """
    Open a plain-text or gzip-compressed file.

    Compression is determined from the .gz filename suffix.
    """

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
    """Validate input paths and numeric arguments."""

    if args.threshold_cm < 0:
        raise SystemExit(
            "ERROR: --threshold-cm must be greater than or equal to zero."
        )

    missing_files = [
        str(input_path)
        for input_path in args.input
        if not input_path.is_file()
    ]

    if missing_files:
        formatted_files = "\n  ".join(missing_files)

        raise SystemExit(
            "ERROR: the following input files do not exist:\n"
            f"  {formatted_files}"
        )

    if args.output.exists() and args.output.is_dir():
        raise SystemExit(
            f"ERROR: the output path is a directory: {args.output}"
        )

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


def initialize_database(
    database_path: Path,
) -> sqlite3.Connection:
    """Create and initialize the temporary SQLite aggregation database."""

    database_path.unlink(missing_ok=True)

    connection = sqlite3.connect(database_path)

    connection.executescript(
        """
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = OFF;
        PRAGMA temp_store = FILE;

        CREATE TABLE pair_summary (
            id1                  TEXT    NOT NULL,
            id2                  TEXT    NOT NULL,
            segment_count        INTEGER NOT NULL,
            total_length         REAL    NOT NULL,
            total_squared_length REAL    NOT NULL,
            minimum_length       REAL    NOT NULL,
            maximum_length       REAL    NOT NULL,
            threshold_count      INTEGER NOT NULL,
            threshold_length     REAL    NOT NULL,

            PRIMARY KEY (id1, id2)
        ) WITHOUT ROWID;
        """
    )

    return connection


def canonicalize_pair(
    id1: str,
    id2: str,
) -> tuple[str, str]:
    """
    Return sample IDs in a consistent order.

    This ensures that ID1-ID2 and ID2-ID1 are treated as the same pair.
    """

    if id1 <= id2:
        return id1, id2

    return id2, id1


def insert_batch(
    connection: sqlite3.Connection,
    batch: list[tuple],
) -> None:
    """Insert or update a batch of IBD segments."""

    if not batch:
        return

    upsert_statement = """
        INSERT INTO pair_summary (
            id1,
            id2,
            segment_count,
            total_length,
            total_squared_length,
            minimum_length,
            maximum_length,
            threshold_count,
            threshold_length
        )
        VALUES (
            ?,
            ?,
            1,
            ?,
            ?,
            ?,
            ?,
            ?,
            ?
        )

        ON CONFLICT (id1, id2)
        DO UPDATE SET
            segment_count =
                pair_summary.segment_count + 1,

            total_length =
                pair_summary.total_length +
                excluded.total_length,

            total_squared_length =
                pair_summary.total_squared_length +
                excluded.total_squared_length,

            minimum_length =
                MIN(
                    pair_summary.minimum_length,
                    excluded.minimum_length
                ),

            maximum_length =
                MAX(
                    pair_summary.maximum_length,
                    excluded.maximum_length
                ),

            threshold_count =
                pair_summary.threshold_count +
                excluded.threshold_count,

            threshold_length =
                pair_summary.threshold_length +
                excluded.threshold_length
    """

    connection.executemany(
        upsert_statement,
        batch,
    )

    connection.commit()
    batch.clear()


def read_hapibd_segments(
    input_paths: Iterable[Path],
    threshold_cm: float,
    connection: sqlite3.Connection,
) -> tuple[int, int]:
    """
    Read Hap-IBD segments and aggregate them in the SQLite database.

    Returns
    -------
    segment_count
        Total number of segment rows processed.

    input_file_count
        Total number of input files processed.
    """

    batch: list[tuple] = []
    total_segment_count = 0
    input_file_count = 0

    for input_path in input_paths:
        input_file_count += 1

        print(
            f"Reading Hap-IBD file: {input_path}",
            file=sys.stderr,
        )

        with open_text_file(input_path) as input_handle:
            for line_number, line in enumerate(
                input_handle,
                start=1,
            ):
                stripped_line = line.rstrip("\r\n")

                if not stripped_line:
                    continue

                fields = stripped_line.split("\t")

                if len(fields) != EXPECTED_HAPIBD_COLUMNS:
                    raise SystemExit(
                        "ERROR: invalid Hap-IBD row.\n"
                        f"File: {input_path}\n"
                        f"Line: {line_number}\n"
                        f"Expected fields: {EXPECTED_HAPIBD_COLUMNS}\n"
                        f"Observed fields: {len(fields)}"
                    )

                sample1 = fields[0]
                sample2 = fields[2]

                if not sample1 or not sample2:
                    raise SystemExit(
                        "ERROR: an empty sample ID was found.\n"
                        f"File: {input_path}\n"
                        f"Line: {line_number}"
                    )

                if sample1 == sample2:
                    raise SystemExit(
                        "ERROR: an IBD segment connects a sample to itself.\n"
                        f"File: {input_path}\n"
                        f"Line: {line_number}\n"
                        f"Sample: {sample1}"
                    )

                try:
                    segment_length_cm = float(fields[7])
                except ValueError as error:
                    raise SystemExit(
                        "ERROR: invalid segment length.\n"
                        f"File: {input_path}\n"
                        f"Line: {line_number}\n"
                        f"Observed value: {fields[7]!r}"
                    ) from error

                if (
                    not math.isfinite(segment_length_cm)
                    or segment_length_cm < 0
                ):
                    raise SystemExit(
                        "ERROR: segment length must be a finite, "
                        "non-negative number.\n"
                        f"File: {input_path}\n"
                        f"Line: {line_number}\n"
                        f"Observed value: {segment_length_cm}"
                    )

                id1, id2 = canonicalize_pair(
                    sample1,
                    sample2,
                )

                passes_threshold = int(
                    segment_length_cm >= threshold_cm
                )

                threshold_length = (
                    segment_length_cm
                    if passes_threshold
                    else 0.0
                )

                batch.append(
                    (
                        id1,
                        id2,
                        segment_length_cm,
                        segment_length_cm**2,
                        segment_length_cm,
                        segment_length_cm,
                        passes_threshold,
                        threshold_length,
                    )
                )

                total_segment_count += 1

                if len(batch) >= BATCH_SIZE:
                    insert_batch(
                        connection,
                        batch,
                    )

    insert_batch(
        connection,
        batch,
    )

    return total_segment_count, input_file_count


def format_number(value: float) -> str:
    """
    Format floating-point values without unnecessary trailing zeros.

    Ten significant digits preserve substantially more precision than
    Hap-IBD normally reports while keeping the output readable.
    """

    return format(value, ".10g")


def write_pair_summary(
    connection: sqlite3.Connection,
    output_path: Path,
    threshold_cm: float,
) -> int:
    """
    Write pair-level summary statistics.

    Returns
    -------
    pair_count
        Number of sample pairs written.
    """

    threshold_label = format_number(threshold_cm)

    header = [
        "ID1",
        "ID2",
        "number_of_IBD_segments",
        "total_IBD_length_cM",
        "mean_IBD_length_cM",
        "min_IBD_length_cM",
        "max_IBD_length_cM",
        "sd_IBD_length_cM",
        f"number_of_segments_ge_{threshold_label}_cM",
        f"total_length_of_segments_ge_{threshold_label}_cM",
    ]

    select_statement = """
        SELECT
            id1,
            id2,
            segment_count,
            total_length,
            total_squared_length,
            minimum_length,
            maximum_length,
            threshold_count,
            threshold_length
        FROM pair_summary
        ORDER BY id1, id2
    """

    pair_count = 0

    with open_text_file(
        output_path,
        mode="wt",
    ) as output_handle:
        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(header)

        for (
            id1,
            id2,
            segment_count,
            total_length,
            total_squared_length,
            minimum_length,
            maximum_length,
            threshold_count,
            threshold_length,
        ) in connection.execute(select_statement):
            mean_length = total_length / segment_count

            # Population variance across all segments observed for the pair.
            variance = (
                total_squared_length / segment_count
                - mean_length**2
            )

            # Protect against very small negative values caused by
            # floating-point rounding.
            variance = max(
                0.0,
                variance,
            )

            standard_deviation = math.sqrt(variance)

            writer.writerow(
                [
                    id1,
                    id2,
                    segment_count,
                    format_number(total_length),
                    format_number(mean_length),
                    format_number(minimum_length),
                    format_number(maximum_length),
                    format_number(standard_deviation),
                    threshold_count,
                    format_number(threshold_length),
                ]
            )

            pair_count += 1

    return pair_count


def main() -> None:
    """Run the pair-level Hap-IBD summarization."""

    args = parse_arguments()
    validate_arguments(args)

    database_path = Path(
        f"{args.output}.sqlite.tmp"
    )

    connection: sqlite3.Connection | None = None

    try:
        connection = initialize_database(
            database_path
        )

        segment_count, input_file_count = read_hapibd_segments(
            input_paths=args.input,
            threshold_cm=args.threshold_cm,
            connection=connection,
        )

        pair_count = write_pair_summary(
            connection=connection,
            output_path=args.output,
            threshold_cm=args.threshold_cm,
        )

        print(
            (
                f"Processed {segment_count:,} IBD segments from "
                f"{input_file_count:,} input file(s)."
            ),
            file=sys.stderr,
        )

        print(
            f"Wrote summary statistics for {pair_count:,} sample pairs.",
            file=sys.stderr,
        )

        print(
            f"Output: {args.output}",
            file=sys.stderr,
        )

    finally:
        if connection is not None:
            connection.close()

        database_path.unlink(
            missing_ok=True
        )


if __name__ == "__main__":
    main()
