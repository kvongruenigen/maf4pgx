# Project Description

Genomic copy number variations (CNV) are a major contributor to the mutation load in cancer. The generation of reference cancer CNV datasets and large-scale CNV data analyses constitute a focus of the Theoretical Cytogenetics and Oncogenomics group. Such analyses are complementary to studies focussing on sequence modifications (e.g. SNVs) detected through NGS techniques.

The study will generate compound variant data through integration of the Progenetix reference CNV data with SNVs from external datasets, most notably from the Cancer Genome Atlas (TCGA; ~11k samples with -omics data), cancer cell lines (cancercelllines.org) and data portals such as cBioPortal.

- Comprehensive integration of non-CNV genomic variation data with Progenetix
- Creation of a GA4GH Beacon w/ deep representation of TCGA (and cell line) data
- Development of Beacon queries and representation for compound genomic variations
- Analysis of compound or alternative variant events in the TCGA dataset and beyond

# Workflow

## Novel data

### Download

Modules: requests, json, re, os

The workflow is able to search for mutation annotation format (MAF) files for masked somatic mutations from TCGA program in the NIC/GDC database and download all files. The download happens in chunks of 1000 files, so the server does not time out. With this a file is created to keep track of the downloaded files. This file is read within a new run of the download code and will exclude files that are already present.

### Unpacking zip files

Modules: 

Afterwards, all the downloads are unpacked until the MAF files and then stored in the data directory in the maf_files folder and empty directories are deleted for cleanliness.

## MAF files available

### Data extraction

Modules: pandas, os, glob, tqdm

*For future use it should be possible to load MAF files into the same directory and the workflow will extract the information from the files and store them together in a file called maf_data.csv. If there is already data available from previous imports, it will be read in and the new data will be compared to the existing, so only additional data will get loaded into the CSV.* ****(Not there yet. Necessary?)****

The extraction script will load in the MAF files stored in the `data/maf_files/` directory and store the following columns in a data frame, that will be put in a list of data frames, which will be put together in the end. The relevant columns are:

- Tumor_Sample_UUID
    - GDC aliquot UUID for tumor sample
- Matched_Norm_Sample_UUID
    - GDC aliquot UUID for matched normal sample
- case_id
    - GDC UUID for the case
- Chromosome
    - The affected chromosome (e.g., chr1)
- Start_Position
    - Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate
- End_Position
    - Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate
- Variant_Classification
    - Translational effect of variant allele
- Variant_Type
    - Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ONP, INS, DEL, or Consolidated)
- Reference_Allele
    - The plus strand reference allele at this position. Includes the deleted sequence for a deletion or "-" for an insertion.
- Tumor_Seq_Allele2
    - Tumor sequencing (discovery) allele 2 - ***Tumor_Seq_Allele1 is always the same as the Reference_Allele***
- Tumor_Sample_Barcode
    - Aliquot barcode for the tumor sample

During the extraction the variable “Tumor_Sample_Barcode” will be labeled correctly as “aliquot_barcode”, since the given barcode belongs to an aliquot of a sample, and only the first 16 characters are kept as sample barcode. (Hierarchy in GDC Data Portal: Samples > Portions > Analytes > Aliquots)

### Conversion

Modules: tidyverse, TCGAutils

With the package ‘TCGAutils’ it is possible to convert barcodes to UUIDs, which is used to obtain the original sample id. This is needed since in progenetix the sample id is used instead of the aliquot id. Furthermore, the columns are renamed to match the variable names used in the bycon package variant import script. Lastly, the temporary file mapfile.tsv is created in the temp directory. The columns kept are:

- case_id
- sample_id
- aliquot_id
- reference_id
- chromosome
- start
- end
- variant_classification
- variant_type

### Mapping

Modules: os, pandas, bycon, pymongo, tqdm, numpy

During the mapping process, the file created beforehand will be loaded and several conventions from progenetix will be applied:

1. The case and sample ids will get the prefix ‘pgx:TCGA.’
2. ‘variant_type’ is renamed in ‘snv_type’
3. ‘chromosome’ is renamed into ‘reference_name’ and the ‘chr’ prefix is discarded
4. The sequence ontology for sequence alteration (SO:0001059) is used to label all variants as SNVs - [http://www.sequenceontology.org/browser/](http://www.sequenceontology.org/browser/)
5. Each SNV gets the sequence ontology for its specific type:
- Single nucleotide polymorphisms (SNP): SO:0001483
- Multiple nucleotide polymorphisms (MNP): SO:0002007 (include DNPs, TNPs, ONPs (≥ 4))
- Deletions (DEL): SO:0000159
- Insertions (INS): SO:0000667
6. Conversion of the ‘1-start, fully-closed’ into a ‘0-start, half-open’ genomic coordinate system as recommended by the global alliance for genomics and health ([GA4GH](http://ga4gh.org)) - [https://genomestandards.org/standards/genome-coordinates/](https://genomestandards.org/standards/genome-coordinates/). This is achieved by:
- Subtracting the value 1 of the ‘start’ variable for SNPs, MNPs, and DELs
- Subtracting the value 1 of the ‘end’ variable for INSs


Additional explanation for the genomic coordinate system:

[https://www.biostars.org/p/84686/](https://www.biostars.org/p/84686/)

Then the sample ids are taken from a set of unique aliquot ids and mapped to the progenetix data base, where the internal biosample id and individual id is retrieved, if the sample id is in the data base. Additionally, an internal callset id is generated and added for each unique aliquot id. The generated and the retrieved variables are then assigned to the corresponding aliquot id.

In the end, the variants that couldn’t be mapped to a sample id will be labeled as new variants and stored in a separate file. The variants ready to be imported and the new variants will be stored as a TSV in the data directory with the following format:

| biosample_id | variant_id | callset_id | individual_id | reference_name | start | end | reference_bases | alternate_bases | variant_classification | variant_state_id | specific_so | aliquot_id | reference_id | case_id | sample_id | variant_types |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |

With this the format for the database import is given and the data can be import into the progenetix MongoDB.

---

## ClinVar

Modules: pandas, numpy, json, xml.etree.ElementTree, os, requests, gzip, io.BytesIO, shutil

ClinVar information is being retrieved with the goal to be included in the variant description in the data base.

The ClinVar workflow step uses the script clinvar.py and needs the `maf_data.csv` file to be present.
In a first step variant names are generated from this file in the format `(Gene):HGVSc (HGVSp_Short)`.

Before mapping the script also checks if mapped variants are already present.

The additional file needed is the [ClinVar variation release](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/):
`ClinVarVariationRelease_00-latest.xml`
This will be downloaded and unpacked if not present.

With the needed prerequisit the XML file is parsed in a 'per element' fashion with `iterparse()`, the starting tag for the element is `VariationArchive` and the resulting element will be searched for the desired entries, according to the [BeaconV2 schema](https://github.com/ga4gh-beacon/beacon-v2/blob/main/models/src/beacon-v2-default-model/genomicVariations/defaultSchema.yaml). In the end the mapped variants are stored in the output file `data/mapped_variants.json`.
An example for the ClinVar XML file can be found [here](https://github.com/ncbi/clinvar/blob/master/sample_xmls/vcv_01.xml#L105).

The resulting objects have the follwing structure:
- variant_name
- variant_type
- alternate_bases
- reference_bases
- identifiers
    - clinvar_ids
    - genomicHGVS_id
    - transcriptHGVS_ids
    - proteinHGVS_ids
    - variant_alternative_ids
- molecular_attributes
    - gene_ids
    - aminoacid_changes
    - molecular_effects
- clinical_interpretations
    - category
    - condition_id
    - clinical_relevance
    - effect
        - id
        - label
    - effect_ids
- location
    - chromosome
    - start
    - stop
    - sequence_id
    - cytogenetic
