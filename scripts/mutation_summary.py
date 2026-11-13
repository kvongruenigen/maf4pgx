#!/usr/bin/env python3
"""Summarize distinct calls and recurrent variants in a public TCGA export."""

import argparse
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
import csv

from oncoprint import VARIANT_KEY, read_calls


def write_table(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def summarize(source, prefix, top):
    samples, genes = read_calls(source, include_unnamed=True)
    calls_by_sample = defaultdict(dict)
    samples_by_variant = defaultdict(set)
    genes_by_variant = defaultdict(set)
    for gene, sample_calls in genes.items():
        for sample, variants in sample_calls.items():
            for key, categories in variants.items():
                calls_by_sample[sample].setdefault(key, set()).update(categories)
                samples_by_variant[key].add(sample)
                if gene:
                    genes_by_variant[key].add(gene)

    prefix.parent.mkdir(parents=True, exist_ok=True)
    sample_rows = []
    for sample in sorted(samples):
        variants = calls_by_sample[sample]
        counts = Counter(
            next(iter(categories)) if len(categories) == 1 else "Multiple consequence classes"
            for categories in variants.values()
        )
        sample_rows.extend((sample, len(variants), category, count)
                           for category, count in sorted(counts.items()))
    write_table(Path(f"{prefix}.sample_calls.tsv"), (
        "sample_barcode", "distinct_calls", "consequence_class", "calls_in_class",
    ), sample_rows)

    recurrent = sorted(
        (key for key, observed in samples_by_variant.items() if len(observed) >= 2),
        key=lambda key: (-len(samples_by_variant[key]), key),
    )
    write_table(Path(f"{prefix}.recurrent_variants.tsv"), (
        *VARIANT_KEY, "genes", "samples_with_call", "samples_in_export",
    ), (
        (*key, ",".join(sorted(genes_by_variant[key])), len(samples_by_variant[key]), len(samples))
        for key in recurrent
    ))

    ranked = sorted((gene for gene in genes if gene), key=lambda gene: (-len(genes[gene]), gene))
    write_table(Path(f"{prefix}.gene_pairs.tsv"), (
        "gene_a", "gene_b", "both", "a_only", "b_only", "neither", "samples_in_export",
    ), (
        (a, b, len(both), len(genes[a]) - len(both), len(genes[b]) - len(both),
         len(samples) - len(genes[a]) - len(genes[b]) + len(both), len(samples))
        for a, b in combinations(ranked[:top], 2)
        for both in [genes[a].keys() & genes[b].keys()]
    ))
    return len(samples), len(recurrent)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export", type=Path, help="All-gene .variants.tsv.gz export")
    parser.add_argument("--out-prefix", type=Path, required=True)
    parser.add_argument("--top", type=int, default=20, help="Number of genes for pair counts")
    args = parser.parse_args()
    if args.top < 2:
        parser.error("--top must be at least 2")
    try:
        samples, recurrent = summarize(args.export, args.out_prefix, args.top)
        print(f"Summarized {samples:,} samples and {recurrent:,} recurrent variants")
    except (OSError, EOFError, ValueError, csv.Error) as error:
        parser.exit(1, f"Mutation summary failed: {error}\n")


if __name__ == "__main__":
    main()
