# Offline showcase fixture

These tables contain eight fictional samples, made-up gene names, and invented coordinates and identifiers. They do not represent TCGA patients or biological findings.

From the repository root, run `python3 examples/showcase/run_demo.py`. To reproduce the checked-in preview, run the same command with `--output-dir` pointing to a scratch directory and compare its `synthetic.oncoprint.svg` with the preview here.

The repeated GENE_A call for SYNTHETIC-01 counts once in that sample. The same allele also occurs in SYNTHETIC-02, so it is recurrent across two samples. The last two samples have calls only in GENE_C, outside the two plotted genes, but still count toward the denominator of eight. The crosswalks produce four exact matches and one example of each missing or ambiguous join status.

## Mutation and copy-number report

The separate [copy-number fixture](copy_number/) uses the real GRCh38 TP53 and APC
locus bounds with fictional calls, segments, sample IDs, and source IDs. It
includes exact boundaries, a partial segment, competing segments, and an allele
reported by two source files.

After [installing the R dependencies](../../docs/public_workflow.md#build-the-interactive-copy-number-report),
run `python3 examples/showcase/run_report.py` from the repo root. Open the printed
HTML path in a browser. To regenerate the checked-in preview, pass
`--output examples/showcase/synthetic.copy-number.html`.
