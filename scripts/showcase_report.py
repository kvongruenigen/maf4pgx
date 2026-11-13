#!/usr/bin/env python3
"""Build an offline mutation and copy-number report from cached GRCh38 inputs."""

import argparse
from collections import defaultdict
import csv
from html import escape
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from uuid import UUID

from audit_sample_join import audit, rows
from combined_oncoprint import TP53
from oncoprint import VARIANT_KEY
from progenetix_pilot import PROTEIN_ALTERING

HERE = Path(__file__).resolve().parent
LOCI = {"TP53": TP53, "APC": ("5", 112707498, 112846239)}
GDC_FIELDS = ("sample_barcode", "gdc_sample_uuid", "segment_file_id",
              "sequence_calls", "copy_state", "download_error")
PGX_FIELDS = ("sample_barcode", "sample_uuid", "biosample_id", "gdc_segment_file_id",
              "sequence_calls", "protein_or_splice_calls", "progenetix_evidence",
              "progenetix_overlapping_intervals", "error")


def indexed(path, fields, key):
    result = {}
    for row in rows(path, fields):
        value = row[key]
        if key == "sample_uuid":
            value = str(UUID(value))
        if not value or value in result:
            raise ValueError(f"{path}: expected unique {key}, found {value!r}")
        result[value] = row
    return result


def variant_evidence(path):
    samples = set()
    calls = {gene: defaultdict(dict) for gene in LOCI}
    columns = ("sample_barcode", "gene", "consequence", "gdc_file_id", *VARIANT_KEY)
    for row in rows(path, columns):
        sample = row["sample_barcode"]
        key = tuple(row[field] for field in VARIANT_KEY)
        if not sample or not all(key) or row["assembly"] != "GRCh38":
            raise ValueError("Expected complete sample and variant keys in a GRCh38 export")
        samples.add(sample)
        if row["gene"] not in LOCI:
            continue
        variant = calls[row["gene"]][sample].setdefault(
            key, {"consequences": set(), "files": set()}
        )
        variant["consequences"].add(row["consequence"])
        if row["gdc_file_id"]:
            variant["files"].add(row["gdc_file_id"])
    if not samples:
        raise ValueError("The variant export is empty")
    return sorted(samples), calls


def load_segments(file_ids, cache):
    if not file_ids:
        return {}
    with TemporaryDirectory(prefix="maf4pgx-ranges-") as directory:
        loci, files, output = (Path(directory) / name for name in ("loci.tsv", "files.txt", "output.json"))
        with loci.open("w", encoding="utf-8", newline="") as target:
            writer = csv.writer(target, delimiter="\t")
            writer.writerow(("gene", "chromosome", "start", "end"))
            writer.writerows((gene, *locus) for gene, locus in LOCI.items())
        files.write_text("\n".join(sorted(file_ids)) + "\n", encoding="utf-8")
        subprocess.run(["Rscript", str(HERE / "showcase_ranges.R"), str(loci), str(files),
                        str(cache), str(output)], check=True)
        return json.loads(output.read_text(encoding="utf-8"))


def gdc_evidence(row, gene, segments):
    file_id = row["segment_file_id"]
    sample_uuid = row["gdc_sample_uuid"]
    evidence = {
        "file": file_id, "uuid": sample_uuid, "mean": None, "segments": [],
        "reason": row["download_error"] or row["copy_state"],
        "pilot_state": row["copy_state"],
    }
    if file_id:
        cached = segments[str(UUID(file_id))]
        overlap = cached.get(gene, cached)
        if not sample_uuid:
            raise ValueError(f"{gene}: segment file {file_id} has no GDC sample UUID")
        if "error" in overlap:
            evidence["reason"] = overlap["error"]
        else:
            evidence.update(overlap)
            if overlap["mean"] is not None and row["copy_state"] == "Download failed":
                evidence["reason"] = "Recovered from local cache after the pilot download failed"
    return evidence


def pgx_evidence(comparison, join, gdc, sample, sequence_count, protein_count):
    _, uuid, biosample, status = join
    evidence = {"join": status, "state": "Unknown", "biosample": biosample, "reason": "",
                "intervals": [], "in_comparison": False}
    if status != "matched":
        evidence["reason"] = {
            "no_gdc_uuid": "No GDC sample UUID",
            "ambiguous_gdc_uuid": "Ambiguous GDC sample UUID",
            "no_progenetix_match": "No matching Progenetix biosample",
            "ambiguous_progenetix_match": "Ambiguous Progenetix biosample",
        }[status]
        return evidence
    if gdc["uuid"] and str(UUID(gdc["uuid"])) != uuid:
        raise ValueError(f"{sample}: GDC segment audit disagrees with the sample crosswalk")
    row = comparison.get(uuid)
    if row is None:
        evidence["reason"] = "No cached Progenetix comparison for this locus"
        return evidence
    if (row["sample_barcode"] != sample or row["biosample_id"] != biosample
            or not gdc["uuid"] or not gdc["file"]
            or str(UUID(row["gdc_segment_file_id"])) != str(UUID(gdc["file"]))):
        raise ValueError(f"{sample}: Progenetix comparison disagrees with the GDC segment audit")
    if (int(row["sequence_calls"]) != sequence_count
            or int(row["protein_or_splice_calls"]) != protein_count):
        raise ValueError(f"{sample}: Progenetix comparison disagrees with the variant export")

    intervals = []
    for interval in filter(None, row["progenetix_overlapping_intervals"].split("; ")):
        span, state_id = interval.split(":", 1)
        start, end = map(int, span.split("-"))
        if not 0 <= start < end:
            raise ValueError(f"{sample}: invalid Beacon interval {span}")
        intervals.append({"start": start + 1, "end": end, "state_id": state_id})
    state = row["progenetix_evidence"]
    evidence.update(
        state=state, intervals=intervals, in_comparison=True, reason=row["error"] or state,
    )
    return evidence


def build(variants, prefix, cache, label, *, synthetic=False):
    samples, calls = variant_evidence(variants)
    gdc = {gene: indexed(Path(f"{prefix}.{gene.lower()}_combined.tsv"), GDC_FIELDS, "sample_barcode")
           for gene in LOCI}
    pgx = {gene: indexed(Path(f"{prefix}.progenetix_{gene.lower()}_cohort.tsv"), PGX_FIELDS, "sample_uuid")
           for gene in LOCI}
    for gene, manifest in gdc.items():
        if set(manifest) != set(samples):
            raise ValueError(f"{gene}: segment audit must cover the all-gene variant sample universe")
    joins = {row[0]: row for row in audit(variants, Path(f"{prefix}.gdc_samples.tsv"),
                                         Path(f"{prefix}.progenetix_samples.tsv"))}
    file_ids = {str(UUID(row["segment_file_id"])) for manifest in gdc.values()
                for row in manifest.values() if row["segment_file_id"]}
    segments = load_segments(file_ids, cache)

    records = {}
    for gene in LOCI:
        records[gene] = {}
        used_comparisons = set()
        for sample in samples:
            variants_for_sample = [
                {**dict(zip(VARIANT_KEY, key)), "consequences": sorted(item["consequences"]),
                 "files": sorted(item["files"])}
                for key, item in sorted(calls[gene].get(sample, {}).items())
            ]
            protein_count = sum(bool(PROTEIN_ALTERING.intersection(item["consequences"]))
                                for item in variants_for_sample)
            gdc_row = gdc[gene][sample]
            if int(gdc_row["sequence_calls"]) != len(variants_for_sample):
                raise ValueError(f"{gene} {sample}: segment audit disagrees with the variant export")
            copy = gdc_evidence(gdc_row, gene, segments)
            pgx_copy = pgx_evidence(pgx[gene], joins[sample], copy, sample,
                                    len(variants_for_sample), protein_count)
            if pgx_copy["in_comparison"]:
                used_comparisons.add(joins[sample][1])
            records[gene][sample] = {
                "variants": variants_for_sample, "protein_count": protein_count,
                "gdc": copy, "pgx": pgx_copy,
            }
        if used_comparisons != pgx[gene].keys():
            raise ValueError(f"{gene}: Progenetix comparison has rows outside the audited matches")
    return {"label": label, "samples": samples, "loci": LOCI, "records": records,
            "synthetic": synthetic}


def write_report(report, output):
    template = (HERE / "showcase_report.html").read_text(encoding="utf-8")
    data = json.dumps(report, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c")
    page = template.replace("{{LABEL}}", escape(report["label"])).replace("{{REPORT_DATA}}", data, 1)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", required=True, type=Path, help="All-gene GRCh38 variant TSV or TSV.GZ")
    parser.add_argument("--pilot-prefix", required=True, type=Path, help="Prefix of the pilot TSVs")
    parser.add_argument("--segment-cache", required=True, type=Path)
    parser.add_argument("--label", required=True, help="Cohort label displayed in the report")
    parser.add_argument("--synthetic", action="store_true", help="Disable external links for fictional IDs")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        report = build(args.variants, args.pilot_prefix, args.segment_cache, args.label,
                       synthetic=args.synthetic)
        write_report(report, args.output)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Showcase report failed: {error}\n")
    print(f"Wrote {args.output} for {len(report['samples']):,} samples")


if __name__ == "__main__":
    main()
