#!/usr/bin/env python3
"""Build an oncoprint and gene counts from an all-gene TCGA variant export."""

import argparse
from collections import defaultdict
import csv
import gzip
from html import escape
from pathlib import Path

VARIANT_KEY = ("assembly", "chromosome", "start", "end", "ref", "alt")
CONSEQUENCES = {
    "Missense_Mutation": "Missense",
    "Nonsense_Mutation": "Stop gained",
    "Nonstop_Mutation": "Stop lost",
    "Frame_Shift_Ins": "Frameshift",
    "Frame_Shift_Del": "Frameshift",
    "In_Frame_Ins": "In-frame indel",
    "In_Frame_Del": "In-frame indel",
    "Splice_Site": "Splice site",
    "Silent": "Synonymous",
}
COLORS = {
    "Missense": "#28618a",
    "Stop gained": "#9b3675",
    "Stop lost": "#b54f2c",
    "Frameshift": "#574497",
    "In-frame indel": "#097570",
    "Splice site": "#906013",
    "Synonymous": "#8397aa",
    "Other": "#5b6370",
    "Multiple consequence classes": "#202531",
    "No reported call": "#e7ecf0",
}


def read_calls(path, *, include_unnamed=False, raw_consequences=False):
    samples = set()
    genes = defaultdict(lambda: defaultdict(dict))
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as source:
        rows = csv.DictReader(source, delimiter="\t")
        required = {"gene", "sample_barcode", "consequence", *VARIANT_KEY}
        missing = required - set(rows.fieldnames or ())
        if missing:
            raise ValueError(f"Missing columns: {', '.join(sorted(missing))}")
        for row in rows:
            sample = row["sample_barcode"]
            key = tuple(row[field] for field in VARIANT_KEY)
            if not sample or not all(key):
                raise ValueError(f"Incomplete sample or variant key at line {rows.line_num}")
            samples.add(sample)
            gene = row["gene"]
            if not gene or gene in {"Unknown", ".", "-"}:
                gene = ""
            if gene or include_unnamed:
                variants = genes[gene][sample]
                consequence = row["consequence"]
                category = consequence if raw_consequences else CONSEQUENCES.get(consequence, "Other")
                variants.setdefault(key, set()).add(category)
    if not samples or not (genes.keys() - {""}):
        raise ValueError("The export needs at least one sample and one named gene")
    return samples, genes


def write_counts(path, ranked, genes, denominator):
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "gene", "samples_with_call", "samples_in_export",
            "percent_of_samples_in_export", "distinct_sample_variant_calls",
        ])
        for gene in ranked:
            calls = genes[gene]
            writer.writerow([
                gene, len(calls), denominator, f"{100 * len(calls) / denominator:.2f}",
                sum(map(len, calls.values())),
            ])


def write_plot(path, title, displayed, samples, genes, *, subtitle=None, extra_colors=None,
               palette=None, background="No reported call", show_counts=True):
    colors = palette or {**COLORS, **(extra_colors or {})}
    # Group similar call patterns, keeping samples with no displayed calls at the right.
    samples = sorted(samples, key=lambda sample: (
        tuple(sample not in genes[gene] for gene in displayed), sample,
    ))
    plot_left, plot_top, row_height = 120, 134 + 22 * ((len(colors) - 1) // 5 - 1), 25
    cell_width = max(3, 780 / len(samples))
    plot_width = cell_width * len(samples)
    width = plot_left + plot_width + 126
    height = plot_top + len(displayed) * row_height + 55
    with path.open("w", encoding="utf-8") as output:
        def label(x, y, value, size=13, weight=400):
            print(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{weight}">'
                  f'{escape(value)}</text>', file=output)

        print(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" '
              f'height="{height}" viewBox="0 0 {width:.2f} {height}" role="img" '
              'aria-labelledby="title desc">', file=output)
        print(f'<title id="title">{escape(title)}</title>'
              '<desc id="desc">Rows are genes and columns are tumor samples. '
              'Colored cells indicate reported mutation consequences.</desc>'
              '<rect width="100%" height="100%" fill="white"/>'
              '<g font-family="system-ui, sans-serif" fill="#202531">', file=output)
        label(24, 34, title, 24, 700)
        label(24, 59, subtitle or f'Top {len(displayed)} genes by samples with a reported call · '
              f'{len(samples):,} samples represented in the export')
        for index, (name, color) in enumerate(colors.items()):
            x, y = 24 + (index % 5) * 188, 76 + (index // 5) * 22
            print(f'<rect x="{x}" y="{y}" width="12" height="12" fill="{color}"/>', file=output)
            label(x + 17, y + 11, name, 11)
        if show_counts:
            label(plot_left + plot_width + 12, plot_top - 8, "Samples (%)", 11, 700)
        for row, gene in enumerate(displayed):
            y = plot_top + row * row_height
            calls = genes[gene]
            label(24, y + 14, gene, 12, 600)
            if show_counts:
                label(plot_left + plot_width + 12, y + 14,
                      f'{len(calls):,} ({100 * len(calls) / len(samples):.1f}%)')
            print(f'<rect x="{plot_left}" y="{y}" width="{plot_width:.2f}" height="18" '
                  f'fill="{colors[background]}"/>', file=output)
            for column, sample in enumerate(samples):
                if sample not in calls:
                    continue
                variants = calls[sample]
                categories = set().union(*variants.values())
                category = (next(iter(categories)) if len(categories) == 1
                            else "Multiple consequence classes")
                tooltip = f'{sample} | {gene} | {len(variants)} distinct calls | {", ".join(sorted(categories))}'
                print(f'<rect x="{plot_left + column * cell_width:.2f}" y="{y}" '
                      f'width="{cell_width - 0.4:.2f}" height="18" fill="{colors[category]}">'
                      f'<title>{escape(tooltip)}</title></rect>', file=output)
        label(plot_left, plot_top + len(displayed) * row_height + 15, "Tumor samples")
        print("</g></svg>", file=output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export", type=Path, help="All-gene .variants.tsv.gz export")
    parser.add_argument("--top", type=int, default=20, help="Genes to plot (default: %(default)s)")
    parser.add_argument("--out-prefix", type=Path, required=True)
    args = parser.parse_args()
    if args.top < 1:
        parser.error("--top must be positive")
    try:
        samples, genes = read_calls(args.export)
        ranked = sorted(genes, key=lambda gene: (-len(genes[gene]), gene))
        args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
        counts = Path(f"{args.out_prefix}.gene_counts.tsv")
        plot = Path(f"{args.out_prefix}.oncoprint.svg")
        write_counts(counts, ranked, genes, len(samples))
        title = args.export.name.removesuffix(".variants.tsv.gz")
        write_plot(plot, f"{title} mutation oncoprint", ranked[:args.top], samples, genes)
        print(f"Wrote {plot} and {counts} for {len(samples):,} samples")
    except (OSError, EOFError, ValueError, csv.Error) as error:
        parser.exit(1, f"Oncoprint failed: {error}\n")


if __name__ == "__main__":
    main()
