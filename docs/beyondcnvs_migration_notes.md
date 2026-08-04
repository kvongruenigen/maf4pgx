# BeyondCNVs migration notes

These notes preserve the useful context from `kvongruenigen/BeyondCNVs` that is not executable workflow code in this repository.

## Repository relationship

`maf4pgx` contains the maintained MAF-to-Progenetix workflow that was developed out of BeyondCNVs. The older BeyondCNVs repository also contains exploratory notebooks, IDE files, Snakemake logs, and notebook checkpoints. Its executable scripts are superseded here:

- `scripts/gdc_maf_download.py` and `scripts/unpack.sh` are represented by `scripts/gdc_downloader.py`.
- `scripts/data_extraction.py` is represented by `scripts/extractor.py`.
- `scripts/mapping_finish.py` is represented by `scripts/maf_curation_pgx.py`.
- `scripts/aliquot_to_sample.R` has a newer implementation here.

## TCGA identifier mapping

The old `dev/Diagnostics.ipynb` notebook records the identifier hierarchy that led to the current barcode conversion step:

- TCGA barcodes follow the hierarchy case, sample, portion, analyte, aliquot.
- MAF `Tumor_Sample_Barcode` values can identify aliquots rather than samples.
- Progenetix matching uses the sample identifier, so aliquot barcodes are shortened to sample barcodes and converted with `TCGAutils::barcodeToUUID`.

The diagnostics notebook also counted Progenetix mapping hits and showed that matching on cases and samples gives different ambiguity profiles. Keep sample-level matching as the default unless the database import model is changed deliberately.

## Unmatched variants

BeyondCNVs wrote variants without a matching Progenetix biosample to `data/varNew.tsv`. This repository now preserves that audit file in `scripts/maf_curation_pgx.py` while keeping `data/varImport.tsv` limited to importable variants.

## HGVSc reliability

The old diagnostics notebooks compared HGVSc-derived reference and alternate bases against sequence fields from MAF data. The working conclusion was that HGVSc strings are not reliable enough to replace sequence columns as the source of truth for import coordinates or alleles.

## Future-work ideas

BeyondCNVs also contained exploratory notes that are not part of the current workflow:

- COSMIC PubMed PMID extraction and linking to Progenetix biosamples.
- A broader workflow blueprint for additional repositories such as cBioPortal, GRCh38 remapping, identifier matching, and MongoDB update/merge steps.
- Early HGVS formatting experiments.

These should be implemented as new explicit workflow steps only if they become active requirements.
