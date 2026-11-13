#!/usr/bin/env python3
"""Export a compact, sample-level variant table from open TCGA MAF files."""

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
import gzip
import io
import json
from pathlib import Path
import re
import sys
import tempfile
import time
from urllib.parse import urlencode
from urllib.request import urlopen
from uuid import UUID

GDC_API = "https://api.gdc.cancer.gov/v0"
MIB = 1024 * 1024
MAF_COLUMNS = {
    "Hugo_Symbol": "gene",
    "NCBI_Build": "assembly",
    "Chromosome": "chromosome",
    "Start_Position": "start",
    "End_Position": "end",
    "Reference_Allele": "ref",
    "Tumor_Seq_Allele2": "alt",
    "Variant_Classification": "consequence",
    "Variant_Type": "variant_type",
    "Tumor_Sample_Barcode": "tumor_barcode",
}
SAMPLE_BARCODE = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}[A-Z])(?:-|$)")


def gdc_files(project, workflow, max_mib):
    """Page through metadata before making any file downloads."""
    filters = {
        "cases.project.project_id": project,
        "access": "open",
        "data_type": "Masked Somatic Mutation",
        "data_format": "MAF",
        "experimental_strategy": "WXS",
        "analysis.workflow_type": workflow,
    }
    params = {
        "filters": json.dumps({
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": field, "value": [value]}}
                for field, value in filters.items()
            ],
        }),
        "fields": "file_id,file_name,file_size",
        "sort": "file_id:asc",
        "size": 500,
    }
    files = []
    total_bytes = 0
    while True:
        params["from"] = len(files)
        with urlopen(f"{GDC_API}/files?{urlencode(params)}", timeout=60) as response:
            page = json.load(response)["data"]
        for hit in page["hits"]:
            hit["file_id"] = str(UUID(hit["file_id"]))
            if not hit["file_name"].endswith((".maf", ".maf.gz")):
                raise ValueError(f"Unexpected MAF filename: {hit['file_name']}")
            if hit["file_size"] <= 0:
                raise ValueError(f"Missing file size for {hit['file_id']}")
            total_bytes += hit["file_size"]
            if total_bytes > max_mib * MIB:
                raise ValueError(
                    f"Selection exceeds the {max_mib} MiB download limit. "
                    "Use a smaller project or increase --max-mib."
                )
            files.append(hit)
        if len(files) >= page["pagination"]["total"]:
            break
        if not page["hits"]:
            raise ValueError("GDC returned an incomplete file listing")
    if not files:
        raise ValueError(f"No open masked MAFs for {project} with workflow {workflow!r}")
    if len({hit["file_id"] for hit in files}) != len(files):
        raise ValueError("GDC returned duplicate file IDs in the listing")
    return files, total_bytes


def write_variants(source, writer, file_id, genes):
    # Skip only the preamble so '#' inside a data field is preserved.
    for line in source:
        if line.strip() and not line.startswith("#"):
            header = next(csv.reader([line], delimiter="\t"))
            break
    else:
        raise ValueError("Missing MAF header")

    missing = MAF_COLUMNS.keys() - set(header)
    if missing:
        raise ValueError(f"Missing MAF columns: {', '.join(sorted(missing))}")

    count = 0
    for row in csv.DictReader(source, fieldnames=header, delimiter="\t"):
        if genes and row["Hugo_Symbol"] not in genes:
            continue
        barcode = row["Tumor_Sample_Barcode"]
        match = SAMPLE_BARCODE.match(barcode or "")
        if not match:
            raise ValueError(f"Invalid TCGA tumor barcode: {barcode!r}")
        writer.writerow([file_id, match[1], *(row[column] for column in MAF_COLUMNS)])
        count += 1
    return count


def fetch_variants(hit, genes):
    file_id = hit["file_id"]
    for attempt in range(3):
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
        try:
            with urlopen(f"{GDC_API}/data/{file_id}", timeout=30) as response:
                compressed = hit["file_name"].endswith(".gz")
                stream = gzip.GzipFile(fileobj=response) if compressed else nullcontext(response)
                with stream as data, io.TextIOWrapper(data, encoding="utf-8-sig") as source:
                    count = write_variants(source, writer, file_id, genes)
            return buffer.getvalue(), count
        except (OSError, EOFError) as error:
            if attempt == 2:
                raise ValueError(f"{file_id}: {error}") from error
            time.sleep(2 ** attempt)
        except (ValueError, csv.Error) as error:
            raise ValueError(f"{file_id}: {error}") from error


def export(files, output, genes, workers=1):
    """Stream the source files and replace the output only after a complete run."""
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with tempfile.TemporaryDirectory(dir=output.parent) as staging:
        staged_output = Path(staging) / "variants.tsv.gz"
        with gzip.open(staged_output, "wt", encoding="utf-8", newline="") as destination:
            writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
            writer.writerow(["gdc_file_id", "sample_barcode", *MAF_COLUMNS.values()])
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for offset in range(0, len(files), workers):
                    batch = files[offset:offset + workers]
                    tasks = [pool.submit(fetch_variants, hit, genes) for hit in batch]
                    for task in tasks:
                        rows, row_count = task.result()
                        destination.write(rows)
                        count += row_count
        staged_output.replace(output)
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project", help="One GDC project, for example TCGA-COAD")
    parser.add_argument(
        "--workflow", default="Aliquot Ensemble Somatic Variant Merging and Masking",
        help="GDC workflow label (default: %(default)s)",
    )
    parser.add_argument("--gene", action="append", default=[], help="Gene symbol, repeatable")
    parser.add_argument("--max-mib", type=int, default=256, help="Download limit (default: %(default)s MiB)")
    parser.add_argument("--workers", type=int, default=1, help="Concurrent downloads (default: %(default)s)")
    parser.add_argument("--plan", action="store_true", help="List files and sizes without downloading")
    parser.add_argument("--output", type=Path, help="Output path (default: data/<PROJECT>.variants.tsv.gz)")
    args = parser.parse_args()
    project = args.project.upper()
    if not re.fullmatch(r"TCGA-[A-Z0-9]+", project):
        parser.error("project must be a TCGA project ID such as TCGA-COAD")
    if args.max_mib <= 0:
        parser.error("--max-mib must be positive")
    if not 1 <= args.workers <= 8:
        parser.error("--workers must be between 1 and 8")
    output = args.output or Path(f"data/{project}.variants.tsv.gz")
    if not args.plan and not output.name.endswith(".tsv.gz"):
        parser.error("--output must end in .tsv.gz")

    try:
        files, total_bytes = gdc_files(project, args.workflow, args.max_mib)
        print(f"{project}: {len(files)} MAFs, {total_bytes / MIB:.1f} MiB", file=sys.stderr)
        if args.plan:
            writer = csv.DictWriter(
                sys.stdout, fieldnames=["file_id", "file_name", "file_size"],
                extrasaction="ignore", delimiter="\t", lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(files)
        else:
            count = export(files, output, {gene.upper() for gene in args.gene}, args.workers)
            print(f"Wrote {count:,} variant rows to {output}", file=sys.stderr)
    except (OSError, EOFError, ValueError, KeyError, csv.Error) as error:
        parser.exit(1, f"GDC export failed: {error}\n")


if __name__ == "__main__":
    main()
