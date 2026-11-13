#!/usr/bin/env python3
"""Compare copy-loss co-occurrence across thresholds on a fixed sample set."""

import argparse
import math
from pathlib import Path

from audit_sample_join import rows
from mutation_summary import write_table

FIELDS = (
    "sample_barcode", "sample_uuid", "gdc_segment_mean", "sequence_calls",
    "protein_or_splice_calls", "progenetix_evidence",
)


def summarize(gene, records, thresholds):
    for threshold in thresholds:
        loss = [row for row in records if float(row["gdc_segment_mean"]) <= -threshold]
        yield (
            gene, threshold, len(records), len(loss),
            sum(int(row["sequence_calls"]) > 0 for row in loss),
            sum(int(row["protein_or_splice_calls"]) > 0 for row in loss),
            sum(float(row["gdc_segment_mean"]) >= threshold for row in records),
            sum(row["progenetix_evidence"] == "Copy loss" for row in records),
            sum(row["progenetix_evidence"] == "Copy loss" for row in loss),
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--locus", nargs=2, action="append", metavar=("GENE", "TSV"), required=True)
    parser.add_argument("--thresholds", nargs="+", type=float, default=[0.15, 0.2, 0.25])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not all(math.isfinite(value) and value > 0 for value in args.thresholds):
        parser.error("Thresholds must be positive and finite")
    result = []
    for gene, source in args.locus:
        records = list(rows(Path(source), FIELDS))
        if not records or len({row["sample_uuid"] for row in records}) != len(records):
            parser.error(f"{source}: expected unique sample UUIDs and at least one record")
        if not all(math.isfinite(float(row["gdc_segment_mean"])) for row in records):
            parser.error(f"{source}: expected finite GDC segment means")
        result.extend(summarize(gene, records, sorted(set(args.thresholds))))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_table(args.output, (
        "gene", "threshold", "evaluable_samples", "gdc_copy_loss",
        "loss_with_any_sequence_call", "loss_with_protein_or_splice_call",
        "gdc_copy_gain", "progenetix_copy_loss", "both_sources_copy_loss",
    ), result)


if __name__ == "__main__":
    main()
