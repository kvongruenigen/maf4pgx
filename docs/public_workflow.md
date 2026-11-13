# Public TCGA workflow

## Public TCGA variant export

Export open, masked whole-exome somatic calls for one TCGA project using Python 3.9 or newer and the [GDC API](https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/). The public scripts use the Python standard library and run independently of Progenetix, MongoDB, R, and ClinVar. For a run without downloads, see the [offline demo](../README.md#try-it-offline).

```sh
# See which files would be downloaded and their total size.
python3 scripts/public_variants.py TCGA-COAD --plan

# Keep only TP53 and APC calls in a compressed local table.
python3 scripts/public_variants.py TCGA-COAD --gene TP53 --gene APC \
  --output data/TCGA-COAD.TP53-APC.variants.tsv.gz
```

The export reads one source MAF at a time by default and saves only `data/TCGA-COAD.variants.tsv.gz`. Use `--workers 4` to overlap transfers while buffering the selected rows from at most four files. It keeps the GDC file UUID, original tumor barcode, derived sample barcode, reference build, original MAF coordinates and alleles, gene, and consequence. Each row represents a source call, so repeated observations across files are retained. The sample barcode can be used to look up a GDC sample UUID for future copy-number joins.

The default download budget is 256 MiB per run, checked against GDC file sizes before fetching MAFs. `--max-mib` changes that budget. Gene filtering reduces the output size, and the selected MAFs are still transferred in full. Omit `--gene` to export all genes, or use `snakemake --cores 1 data/TCGA-COAD.variants.tsv.gz` for the all-gene export. Use a distinct `--output` path to keep multiple gene selections for the same project.

### Mutation oncoprint

Run the export **without** `--gene` so the sample denominator includes samples with calls outside the plotted genes. If you already ran a gene-filtered export at the default path, regenerate the all-gene file first.

```sh
python3 scripts/public_variants.py TCGA-COAD
python3 scripts/oncoprint.py data/TCGA-COAD.variants.tsv.gz --out-prefix data/TCGA-COAD
```

Open `data/TCGA-COAD.oncoprint.svg` in a browser. The plot shows the 20 genes with calls in the most samples (`--top` changes the number). Hover over a colored cell to see the sample, consequence classes, and distinct variant count. The companion `data/TCGA-COAD.gene_counts.tsv` ranks all named genes. Calls are deduplicated within each gene and sample by assembly, chromosome, coordinates, and alleles.

Percentages use all distinct tumor sample barcodes represented in the all-gene export, including samples without calls in the plotted genes. Samples with zero reported calls across the entire MAF export are outside that denominator. Counts include synonymous calls.

The equivalent workflow target is `snakemake --cores 1 data/TCGA-COAD.oncoprint.svg`.

### Mutation summaries

```sh
python3 scripts/mutation_summary.py data/TCGA-COAD.variants.tsv.gz --out-prefix data/TCGA-COAD
```

- `sample_calls.tsv` gives each sample's distinct call total and the breakdown by consequence class. The total repeats on each consequence row.
- `recurrent_variants.tsv` lists exact MAF alleles reported in two or more samples. Alleles use the source coordinates and representation.
- `gene_pairs.tsv` gives the both, either alone, and neither sample counts for each pair among the 20 most frequently called genes. Use `--top` to change the number of genes.

These summaries use the all-gene export and count synonymous calls. Pair counts describe observed co-occurrence among samples represented in the export. Run `snakemake --cores 1 data/TCGA-COAD.gene_pairs.tsv` to generate all three tables.

### Audit a public CNV sample join

Prepare two sample crosswalks with columns `sample_barcode, sample_uuid` from GDC and `sample_uuid, biosample_id` from Progenetix. Both must be tab-separated and use **sample-level UUIDs**. Use the Progenetix samples with CNV data for the selected cohort.

```sh
python3 scripts/audit_sample_join.py data/TCGA-COAD.variants.tsv.gz \
  data/gdc_samples.tsv data/progenetix_samples.tsv --output data/TCGA-COAD.sample_join.tsv
```

The audit emits one row per mutation-export sample and records missing or ambiguous mappings. Only rows with status `matched` have a unique UUID and an unambiguous Progenetix biosample. The audit checks identifiers. CNV coverage, assembly, and segment states still need validation before plotting combined events.

### Experimental GDC-only TP53 pilot

This separate pilot uses GDC for both inputs. It does not join to Progenetix and is outside the offline showcase. Run it against the all-gene GDC export.

```sh
python3 scripts/combined_oncoprint.py TCGA-COAD data/TCGA-COAD.variants.tsv.gz --out-prefix data/TCGA-COAD
```

The pilot uses open GDC masked DNAcopy segments and exact GDC sample barcodes. It records the sample UUID and segment file ID for each unique match. The TSV audits all mutation-export samples, including ambiguous or missing files, failed downloads, and unresolved locus coverage. The SVG shows samples where one segment covers the entire GDC GRCh38 TP53 gene span. The default gain and loss cutoffs are segment means of +0.2 and −0.2, inclusive. Change these with `--threshold`. The sequence row includes synonymous calls. Co-occurrence in this plot is a descriptive observation.

### Compare TP53 with Progenetix

After generating the combined GDC TSV, run the public Beacon pilot with Python 3 and curl.

```sh
python3 scripts/progenetix_pilot.py TCGA-COAD data/TCGA-COAD.variants.tsv.gz \
  data/TCGA-COAD.tp53_combined.tsv --out-prefix data/TCGA-COAD
```

The workflow queries Progenetix staging for `pgx:TCGA-COAD`. It verifies the project reference, keeps neoplastic biosamples, and joins their UUID-formatted names to GDC sample UUIDs from the segment audit. It queries the GRCh38 TP53 interval for every unique match with resolved GDC coverage. Use `--limit` for a smaller run.

- `.progenetix_tp53_cohort.tsv` records each sample’s GDC segment mean, Progenetix overlapping intervals and state IDs, and both sequence call counts.
- `.progenetix_tp53_cohort.svg` shows the sequence, GDC copy-number, and Progenetix copy-number tracks.
- `.tp53_copy_comparison.tsv` tabulates source combinations. Partial overlaps, no reported spanning CNV, and query failures each have their own category. Progenetix coverage remains unknown in these categories.
- `.tp53_cooccurrence.tsv` compares GDC TP53 loss with any sequence call and with protein-altering or splice-site calls. These are descriptive same-sample counts.

### Compare two GRCh38 loci

Both locus commands accept a gene symbol and a 1-based inclusive GRCh38 interval. The Beacon command also takes the versioned chromosome accession. For APC, the [NCBI GRCh38.p14 gene record](https://www.ncbi.nlm.nih.gov/gene/324) gives `NC_000005.10:112707498–112846239`. The TP53 default retains the original GDC interval.

Reuse the segment manifest and verified Progenetix crosswalk from the same cohort. The optional cache directories hold original segment files and complete region responses so reruns can reuse the downloaded data.

```sh
python3 scripts/combined_oncoprint.py TCGA-COAD data/TCGA-COAD.variants.tsv.gz \
  --gene APC --interval 5 112707498 112846239 \
  --segment-manifest data/TCGA-COAD.tp53_combined.tsv \
  --segment-cache data/cache/gdc_segments --out-prefix data/TCGA-COAD

python3 scripts/progenetix_pilot.py TCGA-COAD data/TCGA-COAD.variants.tsv.gz \
  data/TCGA-COAD.apc_combined.tsv --gene APC --interval 5 112707498 112846239 \
  --refseq NC_000005.10 --crosswalk data/TCGA-COAD.progenetix_samples.tsv \
  --beacon-cache data/cache/beacon --out-prefix data/TCGA-COAD

python3 scripts/locus_sensitivity.py \
  --locus TP53 data/TCGA-COAD.progenetix_tp53_cohort.tsv \
  --locus APC data/TCGA-COAD.progenetix_apc_cohort.tsv \
  --output data/TCGA-COAD.locus_sensitivity.tsv
```

The sensitivity table keeps the matched, GDC-evaluable sample set fixed for each locus while applying GDC cutoffs of ±0.15, ±0.20, and ±0.25 to the same segment means. Loci can have different evaluable samples.

### Build the interactive copy-number report

Install Python 3.9 or newer and an R version supported by Bioconductor. Install
the R packages once:

```sh
Rscript -e 'install.packages(c("BiocManager", "jsonlite"), repos="https://cloud.r-project.org"); BiocManager::install("GenomicRanges", ask=FALSE, update=FALSE)'
```

The report build itself is offline. It reads the all-gene variant export,
both combined GDC TSVs, both Progenetix cohort TSVs, the two sample crosswalk
TSVs, and the segment cache produced by the pilot commands above. To populate
the TP53 cache on a fresh run, include `--segment-cache data/cache/gdc_segments`
on the TP53 GDC command. Use the same cache when running APC.

```sh
python3 -B scripts/showcase_report.py \
  --variants data/TCGA-COAD.variants.tsv.gz \
  --pilot-prefix data/TCGA-COAD \
  --segment-cache data/cache/gdc_segments \
  --label "TCGA-COAD · Colon adenocarcinoma" \
  --output data/TCGA-COAD.copy-number.html
```

Open the resulting HTML in a browser. The oncoprint displays every sample in
pages of 24. Use the sample selector to restrict the display to samples with a
cached comparison for either gene. The source comparison columns use each
gene's own matched, GDC-evaluable subset. They include an explicit count for
unknown Progenetix evidence. Selecting a cell opens the underlying calls,
overlapping intervals, and source file identifiers.

The report recomputes GDC coverage from the cache with `GenomicRanges` and
retains the original pilot state alongside the recomputed evidence. In the local
440-sample run, the cache resolved ten TP53 pilot download failures. The GDC
co-occurrence counts at ±0.20 were 161/398 for TP53 and 68/392 for APC. In the
original Progenetix comparison subsets they were 156/386 and 67/390,
respectively. The source comparison counts reconciled with the pilot sensitivity
table at all three cutoffs. Counts depend on the source files and caches used.
