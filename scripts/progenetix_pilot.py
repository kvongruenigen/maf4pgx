#!/usr/bin/env python3
"""Compare TP53 copy-number evidence from public Progenetix and GDC records."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import re
from urllib.parse import urlencode
from uuid import UUID

from audit_sample_join import audit, rows
from combined_oncoprint import TP53, add_locus_args, parse_locus, store_cache
from mutation_summary import write_table
from oncoprint import read_calls, write_plot

BASE = "https://staging.progenetix.org/beacon"
STATES = {
    "EFO:0030067": "Copy loss", "EFO:0030068": "Copy loss", "EFO:0020073": "Copy loss",
    "EFO:0030070": "Copy gain", "EFO:0030071": "Copy gain", "EFO:0030072": "Copy gain",
}


def query(route, cache_path=None, **params):
    url = f"{BASE}/{route}/?{urlencode(params)}"
    try:
        response = cache_path.read_bytes() if cache_path and cache_path.exists() else subprocess.check_output([
            "curl", "-fsSL", "--max-time", "25", "--retry", "2", url,
        ])
    except subprocess.CalledProcessError as error:
        raise ValueError(f"Beacon request failed for {route} (curl exit {error.returncode})") from error
    data = json.loads(response)
    if data.get("error"):
        raise ValueError(f"Beacon returned an error for {route}")
    records = [row for result in data.get("response", {}).get("resultSets", [])
               for row in result.get("results", [])]
    if len(records) != data["responseSummary"]["numTotalResults"]:
        raise ValueError(f"Incomplete Beacon response for {route}")
    if cache_path and not cache_path.exists():
        store_cache(cache_path, response)
    return records


def crosswalk(project):
    term = f"pgx:{project}"
    result = []
    seen = set()
    for record in query("biosamples", filters=term, limit=1000):
        if term not in {ref["id"] for ref in record["externalReferences"]}:
            raise ValueError("Beacon returned a sample outside the requested project")
        biosample = record["id"]
        if biosample in seen or not re.fullmatch(r"pgxbs-[a-z0-9]+", biosample):
            raise ValueError(f"Invalid or duplicate biosample ID: {biosample}")
        seen.add(biosample)
        if record.get("biosampleStatus", {}).get("id") == "EFO:0009656":
            result.append((str(UUID(record["name"])), biosample))
    return sorted(result)


def locus_evidence(records, locus=TP53, refseq="NC_000017.11"):
    states, ids, intervals = set(), set(), set()
    partial = False
    _, first, last = locus
    for record in records:
        allele = record["variation"]
        variation = allele["variation"]
        if variation["type"] != "CopyNumberChange":
            continue
        location = variation["location"]
        if location.get("sequenceId") != f"refseq:{refseq}":
            continue
        # Beacon locations are 0-based half-open. TP53 is 1-based inclusive.
        start, end = location["start"], location["end"]
        if start >= last or end <= first - 1:
            continue
        state_id = allele["variantState"]["id"]
        intervals.add((start, end, state_id))
        if start <= first - 1 and end >= last:
            states.add(STATES.get(state_id, "Other state"))
            ids.add(state_id)
        else:
            partial = True
    evidence = ", ".join(sorted(states))
    if not evidence:
        evidence = "Partial overlap only" if partial else "No spanning CNV reported"
    return evidence, ", ".join(sorted(ids)), "; ".join(
        f"{start}-{end}:{state}" for start, end, state in sorted(intervals)
    )


def compare(item, locus, refseq, cache):
    sample, uuid, biosample, _, gdc, protein_calls = item
    try:
        chromosome, start, end = locus
        cached = cache / f"{biosample}.{refseq}.{start}-{end}.json" if cache else None
        evidence, states, intervals = locus_evidence(query(
            f"biosamples/{biosample}/g_variants", assemblyId="GRCh38", limit=5000,
            referenceName=chromosome, start=start - 1, end=end, cache_path=cached,
        ), locus, refseq)
        error = ""
    except (OSError, ValueError, KeyError, TypeError) as exc:
        evidence, states, intervals, error = "Query failed", "", "", str(exc)
    return (
        sample, uuid, biosample, gdc["segment_file_id"], gdc["sequence_calls"],
        gdc["copy_state"], evidence, states, error, protein_calls, gdc["segment_mean"], intervals,
    )


PROTEIN_ALTERING = {
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
    "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del",
    "Splice_Site", "Translation_Start_Site", "Start_Codon_Del",
    "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins",
}


def plot_comparison(prefix, project, comparison, gene):
    samples = {row[0] for row in comparison}
    tracks = {
        f"{gene} sequence": {
            row[0]: {"sequence": {"Protein or splice" if row[9] else "Other sequence call"}}
            for row in comparison if int(row[4])
        },
        "GDC copy": {row[0]: {"gdc": {row[5]}} for row in comparison},
        "Progenetix copy": {row[0]: {"pgx": {row[6]}} for row in comparison},
    }
    colors = {
        "Protein or splice": "#28618a", "Other sequence call": "#8397aa",
        "No sequence call": "#e7ecf0", "Copy loss": "#b54f2c",
        "Copy gain": "#097570", "Within threshold": "#b8c9db",
        "No spanning CNV reported": "#ece1c7", "Partial overlap only": "#906013",
        "Other state": "#705199", "Query failed": "#45454e",
    }
    for row in comparison:
        colors.setdefault(row[6], "#a54d83")
    write_plot(
        Path(f"{prefix}.progenetix_{gene.lower()}_cohort.svg"),
        f"{project} · {gene} sequence and copy number", list(tracks), samples, tracks,
        subtitle=f"{len(samples)} UUID-matched samples with resolved GDC {gene} coverage",
        palette=colors, background="No sequence call", show_counts=False,
    )


def write_summary(prefix, comparison, gdc, protein, gene):
    output = []
    for label, cohort in (
        ("Resolved GDC coverage", gdc.values()),
        ("UUID-matched comparison", (gdc[row[0]] for row in comparison)),
    ):
        eligible = [row for row in cohort if row["copy_state"] in {
            "Copy loss", "Copy gain", "Within threshold",
        }]
        loss = [row for row in eligible if row["copy_state"] == "Copy loss"]
        output.append((label, len(eligible), len(loss),
                       sum(int(row["sequence_calls"]) > 0 for row in loss),
                       sum(protein.get(row["sample_barcode"], 0) > 0 for row in loss)))
    write_table(Path(f"{prefix}.{gene.lower()}_cooccurrence.tsv"), (
        "cohort", "evaluable_samples", "gdc_copy_loss",
        "loss_with_any_sequence_call", "loss_with_protein_or_splice_call",
    ), output)
    counts = Counter((row[5], row[6]) for row in comparison)
    write_table(Path(f"{prefix}.{gene.lower()}_copy_comparison.tsv"), (
        "gdc_copy_state", "progenetix_evidence", "samples",
    ), ((*key, count) for key, count in sorted(counts.items())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project")
    parser.add_argument("variants", type=Path)
    parser.add_argument("gdc_audit", type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--limit", type=int, default=None, help="Optional cap on matched samples")
    add_locus_args(parser)
    parser.add_argument("--refseq", help="Versioned GRCh38 RefSeq chromosome accession")
    parser.add_argument("--crosswalk", type=Path, help="Previously verified Progenetix sample TSV")
    parser.add_argument("--beacon-cache", type=Path)
    args = parser.parse_args()
    gene = args.gene
    locus = parse_locus(parser, args)
    refseq = args.refseq or ("NC_000017.11" if locus[0] == "17" else None)
    if not refseq or not re.fullmatch(r"NC_[0-9]+\.[0-9]+", refseq):
        parser.error("Provide a versioned --refseq accession for this chromosome")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if not re.fullmatch(r"TCGA-[A-Z0-9]+", args.project):
        parser.error("Expected a TCGA project ID")
    try:
        prefix = args.out_prefix
        prefix.parent.mkdir(parents=True, exist_ok=True)
        gdc = {row["sample_barcode"]: row for row in rows(args.gdc_audit, [
            "sample_barcode", "gdc_sample_uuid", "copy_state", "segment_file_id", "sequence_calls", "segment_mean",
        ])}
        gdc_path = Path(f"{prefix}.gdc_samples.tsv")
        pgx_path = args.crosswalk or Path(f"{prefix}.progenetix_samples.tsv")
        write_table(gdc_path, ("sample_barcode", "sample_uuid"), (
            (sample, row["gdc_sample_uuid"]) for sample, row in sorted(gdc.items())
            if row["gdc_sample_uuid"]
        ))
        if not args.crosswalk:
            write_table(pgx_path, ("sample_uuid", "biosample_id"), crosswalk(args.project))
        matches = list(audit(args.variants, gdc_path, pgx_path))
        write_table(Path(f"{prefix}.sample_join.tsv"), (
            "sample_barcode", "gdc_sample_uuids", "progenetix_biosample_ids", "status",
        ), matches)
        print(dict(Counter(row[-1] for row in matches)), flush=True)
        _, genes = read_calls(args.variants, raw_consequences=True)
        gene_calls = genes.get(gene, {})
        protein = {
            sample: sum(bool(categories & PROTEIN_ALTERING) for categories in variants.values())
            for sample, variants in gene_calls.items()
        }
        for sample, row in gdc.items():
            if int(row["sequence_calls"]) != len(gene_calls.get(sample, {})):
                raise ValueError(f"GDC audit sequence count differs for {gene} sample {sample}")
        selected = [(*row, gdc[row[0]], protein.get(row[0], 0)) for row in matches if row[-1] == "matched"
                    and gdc[row[0]]["copy_state"] in {"Copy gain", "Copy loss", "Within threshold"}]
        if not selected:
            raise ValueError(f"No unique matches with resolved GDC {gene} coverage")
        with ThreadPoolExecutor(max_workers=4) as pool:
            comparison = list(pool.map(lambda item: compare(item, locus, refseq, args.beacon_cache), selected[:args.limit]))
        write_table(Path(f"{prefix}.progenetix_{gene.lower()}_cohort.tsv"), (
            "sample_barcode", "sample_uuid", "biosample_id", "gdc_segment_file_id",
            "sequence_calls", "gdc_copy_state", "progenetix_evidence", "progenetix_state_ids", "error",
            "protein_or_splice_calls", "gdc_segment_mean", "progenetix_overlapping_intervals",
        ), comparison)
        plot_comparison(prefix, args.project, comparison, gene)
        write_summary(prefix, comparison, gdc, protein, gene)
        print(dict(Counter(row[6] for row in comparison)))
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Progenetix pilot failed: {error}\n")


if __name__ == "__main__":
    main()
