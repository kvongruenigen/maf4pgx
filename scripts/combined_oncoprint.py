#!/usr/bin/env python3
"""Plot sequence calls and masked DNAcopy segments for a GRCh38 locus."""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import io
import json
import math
from pathlib import Path
import re
import time
from tempfile import NamedTemporaryFile
from urllib.parse import urlencode
from urllib.request import urlopen
from uuid import UUID

from audit_sample_join import rows
from mutation_summary import write_table
from oncoprint import read_calls, write_plot
from public_variants import GDC_API

# GDC Portal GRCh38 TP53 gene span, 1-based inclusive.
TP53 = ("17", 7661779, 7687538)


def store_cache(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=path.parent, delete=False) as target:
        target.write(content)
    Path(target.name).replace(path)


def add_locus_args(parser):
    parser.add_argument("--gene", default="TP53")
    parser.add_argument("--interval", nargs=3, metavar=("CHROM", "START", "END"),
                        help="GRCh38 interval, 1-based inclusive")


def parse_locus(parser, args):
    if not re.fullmatch(r"[A-Z0-9][A-Z0-9-]*", args.gene):
        parser.error("Expected an uppercase gene symbol")
    if args.interval is None and args.gene == "TP53":
        return TP53
    if args.interval is None:
        parser.error("--interval is required for this gene")
    chromosome, start, end = args.interval
    try:
        start, end = int(start), int(end)
        if 1 <= start <= end:
            return chromosome.removeprefix("chr"), start, end
    except ValueError:
        pass
    parser.error("Expected a positive 1-based inclusive interval")


def get(url):
    for attempt in range(3):
        try:
            with urlopen(url, timeout=20) as response:
                return response.read()
        except OSError:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def segment_files(project):
    filters = {
        "cases.project.project_id": project, "access": "open",
        "data_type": "Masked Copy Number Segment", "analysis.workflow_type": "DNAcopy",
    }
    params = {
        "filters": json.dumps({"op": "and", "content": [
            {"op": "in", "content": {"field": field, "value": [value]}}
            for field, value in filters.items()
        ]}),
        "fields": "file_id,file_size,cases.samples.submitter_id,cases.samples.sample_id",
        "size": 500, "sort": "file_id:asc",
    }
    offset = 0
    while True:
        params["from"] = offset
        page = json.loads(get(f"{GDC_API}/files?{urlencode(params)}"))["data"]
        yield from page["hits"]
        offset += len(page["hits"])
        if offset >= page["pagination"]["total"]:
            break
        if not page["hits"]:
            raise ValueError("Incomplete GDC file listing")


def classify(text, locus, threshold):
    chromosome, first, last = locus
    means = []
    partial = False
    segments = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not {"Chromosome", "Start", "End", "Segment_Mean"}.issubset(segments.fieldnames or ()):
        raise ValueError("Invalid GDC segment header")
    for row in segments:
        if row["Chromosome"].removeprefix("chr") != chromosome:
            continue
        start, end = int(row["Start"]), int(row["End"])
        if start <= first and last <= end:
            mean = float(row["Segment_Mean"])
            if not math.isfinite(mean):
                raise ValueError("Nonfinite segment mean")
            means.append(mean)
        elif start <= last and first <= end:
            partial = True
    if len(means) != 1 or partial:
        return "Unresolved coverage", ""
    mean = means[0]
    if mean <= -threshold:
        return "Copy loss", mean
    if mean >= threshold:
        return "Copy gain", mean
    return "Within threshold", mean


def run(project, source, prefix, threshold, gene="TP53", locus=TP53, cache=None, manifest=None):
    if {row["assembly"] for row in rows(source, ["assembly"])} != {"GRCh38"}:
        raise ValueError("The variant export must use GRCh38")
    samples, genes = read_calls(source)
    candidates = defaultdict(dict)
    if manifest:
        selected = {
            row["sample_barcode"]: (str(UUID(row["segment_file_id"])), (row["gdc_sample_uuid"], 0))
            for row in rows(manifest, ["sample_barcode", "segment_file_id", "gdc_sample_uuid"])
            if row["sample_barcode"] in samples and row["segment_file_id"] and row["gdc_sample_uuid"]
        }
        for row in rows(manifest, ["sample_barcode", "copy_state"]):
            if row["copy_state"] == "Ambiguous files":
                candidates[row["sample_barcode"]]["ambiguous"] = None
    else:
        for hit in segment_files(project):
            identities = {(sample["submitter_id"], sample["sample_id"])
                          for case in hit.get("cases", []) for sample in case.get("samples", [])}
            if len(identities) == 1:
                barcode, uuid = next(iter(identities))
                if barcode in samples:
                    candidates[barcode][hit["file_id"]] = (uuid, hit["file_size"])
        selected = {sample: next(iter(files.items())) for sample, files in candidates.items()
                    if len(files) == 1}
    if not selected:
        raise ValueError("No uniquely matched segment files")
    size = sum(size for _, (_, size) in selected.values())
    if size > 32 * 1024**2:
        raise ValueError("Segment files exceed the 32 MiB transfer budget")
    print(f"{len(selected)} unique files" if manifest else
          f"{len(selected)} unique files, {size / 1024**2:.2f} MiB", flush=True)

    def download(item):
        sample, (file_id, (uuid, _)) = item
        try:
            cache_path = cache / f"{file_id}.tsv" if cache else None
            content = cache_path.read_bytes() if cache_path and cache_path.exists() else get(f"{GDC_API}/data/{file_id}")
            text = content.decode("utf-8-sig")
            state, mean = classify(text, locus, threshold)
            if cache_path and not cache_path.exists():
                store_cache(cache_path, content)
        except (OSError, ValueError, KeyError) as error:
            return sample, (uuid, file_id, "Download failed", "", str(error))
        return sample, (uuid, file_id, state, mean, "")

    with ThreadPoolExecutor(max_workers=8) as pool:
        results = dict(pool.map(download, sorted(selected.items())))
    calls = genes.get(gene, {})
    prefix.parent.mkdir(parents=True, exist_ok=True)
    records = []
    for sample in sorted(samples):
        if sample in results:
            uuid, file_id, state, mean, error = results[sample]
        else:
            uuid, file_id, mean, error = "", "", "", ""
            state = "Ambiguous files" if candidates[sample] else "No matched file"
        records.append((sample, uuid, file_id, len(calls.get(sample, {})), state, mean, error))
    write_table(Path(f"{prefix}.{gene.lower()}_combined.tsv"), (
        "sample_barcode", "gdc_sample_uuid", "segment_file_id", "sequence_calls",
        "copy_state", "segment_mean", "download_error",
    ), records)

    eligible = {sample for sample, _, _, _, state, _, _ in records
                if state in {"Copy loss", "Copy gain", "Within threshold"}}
    if not eligible:
        raise ValueError(f"No samples have resolved {gene} segment coverage. See the audit table.")
    plot_data = {
        f"{gene} sequence": {sample: variants for sample, variants in calls.items() if sample in eligible},
        f"{gene} copy": {sample: {"segment": {results[sample][2]}} for sample in eligible
                      if results[sample][2] != "Within threshold"},
    }
    title = f"{project} · {gene} sequence and copy number"
    write_plot(Path(f"{prefix}.{gene.lower()}_combined.svg"), title, list(plot_data), eligible, plot_data,
               subtitle=f"{len(eligible)} samples with resolved coverage · DNAcopy segment mean ±{threshold:g}",
               extra_colors={"Copy loss": "#b54f2c", "Copy gain": "#097570"})
    print(dict(Counter(row[4] for row in records)), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project")
    parser.add_argument("variants", type=Path)
    parser.add_argument("--out-prefix", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.2)
    add_locus_args(parser)
    parser.add_argument("--segment-cache", type=Path)
    parser.add_argument("--segment-manifest", type=Path, help="Prior GDC locus audit for the same cohort")
    args = parser.parse_args()
    locus = parse_locus(parser, args)
    if not math.isfinite(args.threshold) or args.threshold <= 0:
        parser.error("--threshold must be positive and finite")
    try:
        run(args.project, args.variants, args.out_prefix, args.threshold, args.gene, locus,
            args.segment_cache, args.segment_manifest)
    except (OSError, ValueError, KeyError, csv.Error) as error:
        parser.exit(1, f"Combined oncoprint failed: {error}\n")


if __name__ == "__main__":
    main()
