#!/usr/bin/env python3
"""Audit exact TCGA sample UUID matches before combining mutation and CNV data."""

import argparse
from collections import defaultdict
from pathlib import Path
from uuid import UUID
import csv
import gzip


def rows(path, columns):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        missing = set(columns) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {', '.join(sorted(missing))}")
        yield from reader


def audit(variants, gdc, progenetix):
    samples = {row["sample_barcode"] for row in rows(variants, ["sample_barcode"])}
    uuids = defaultdict(set)
    biosamples = defaultdict(set)
    reverse = defaultdict(set)
    for row in rows(gdc, ["sample_barcode", "sample_uuid"]):
        if row["sample_barcode"] in samples:
            uuids[row["sample_barcode"]].add(str(UUID(row["sample_uuid"])))
    for row in rows(progenetix, ["sample_uuid", "biosample_id"]):
        uuid = str(UUID(row["sample_uuid"]))
        biosample = row["biosample_id"].strip()
        if not biosample:
            raise ValueError(f"{progenetix}: missing biosample ID for {uuid}")
        biosamples[uuid].add(biosample)
        reverse[biosample].add(uuid)

    for sample in sorted(samples):
        ids = uuids[sample]
        candidates = set().union(*(biosamples[uuid] for uuid in ids))
        if not ids:
            status = "no_gdc_uuid"
        elif len(ids) > 1:
            status = "ambiguous_gdc_uuid"
        elif not candidates:
            status = "no_progenetix_match"
        elif len(candidates) > 1 or any(len(reverse[item]) > 1 for item in candidates):
            status = "ambiguous_progenetix_match"
        else:
            status = "matched"
        yield sample, ",".join(sorted(ids)), ",".join(sorted(candidates)), status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variants", type=Path)
    parser.add_argument("gdc_samples", type=Path, help="TSV with sample_barcode and sample_uuid")
    parser.add_argument("progenetix_samples", type=Path, help="TSV with sample_uuid and biosample_id")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = list(audit(args.variants, args.gdc_samples, args.progenetix_samples))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as target:
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample_barcode", "gdc_sample_uuids", "progenetix_biosample_ids", "status"])
            writer.writerows(result)
        matched = sum(row[-1] == "matched" for row in result)
        print(f"Matched {matched:,} of {len(result):,} samples. Wrote {args.output}")
    except (OSError, EOFError, ValueError, csv.Error) as error:
        parser.exit(1, f"Sample join audit failed: {error}\n")


if __name__ == "__main__":
    main()
