# Beyond CNVs - Exploring compound mutational hits in reference datasets
### Master project in development
Genomic copy number variations (CNV) are a major contributor to the mutation load in cancer. 
The generation of reference cancer CNV datasets and large-scale CNV data analyses constitute a focus of the Theoretical Cytogenetics and Oncogenomics group. 
Such analyses are complementary to studies focussing on sequence modifications (e.g. SNVs) detected through NGS techniques.

The study will generate compound variant data through integration of the [Progenetix](https://www.progenetix.org "Progenetix Homepage") reference CNV data with SNVs from external datasets, most notably from the Cancer Genome Atlas (TCGA; ~11k samples with -omics data), 
cancer cell lines (cancercelllines.org) and data portals such as cBioPortal.

__Expected outcomes:__

* Comprehensive integration of non-CNV genomic variation data with Progenetix
* Creation of a GA4GH Beacon w/ deep representation of TCGA (and cell line) data
* Development of Beacon queries and representation for compound genomic variations
* Analysis of compound or alternative variant events in the TCGA dataset and beyond



<a name="workflow"></a>
# Workflow
### Contents
- [ Novel data ](#novelData)
    - [ Data mining](#dataMining)
- [ MAF files available ](#mafFiles)
    - [ Data extraction ](#dataExtraction)
    - [ Conversion ](#conversion)
    - [ Curation ](#curation)
    - [ Clinvar annotation mining ](#clinvar)
 ```mermaid
%%{
  init: {
    'flowchart': {'defaultRenderer': 'elk'},
    'theme': 'base',
    'themeVariables': {
      'primaryColor': 'white',
      'primaryTextColor': 'black',
      'primaryBorderColor': 'lightblue',
      'fontFamily': 'Helvetica',
      'lineColor': 'gray',
      'tertiaryColor': 'white'
    }
  }
}%%
graph TD
	%% Objects
	gdc("Data acquisition from GDC API")
	data("Data extraction & Collection")
	duplicates("Removing duplicates (BROAD Filter)")
	barcode("Barcode conversion")
	mapping("Mapping to sample UUID in Progenetix")
	clinvar("ClinVar annotation")
	maf("MAF annotation")
	import("Import file")
	progenetix("Progenetix")
	sampleuuid("Sample UUID")
	samplebarcode("Sample Barcode")
	aliquotbarcode("Aliquot Barcode")

	%% Relationships
	gdc ==> data
    MAF:::sub
	subgraph MAF
	data ==> duplicates ==> barcode ==> mapping
	clinvar ==> maf	
	mapping ==> maf
	mapping ==>|biosample_id\nindividual_id|import

    subgraph Conversion
		aliquotbarcode -.-> samplebarcode -.-> sampleuuid
	end
		Conversion -.-> barcode
		barcode -.-> Conversion
	end
	maf ==>|Variation annotation|progenetix	
	import ==>|Variation import|progenetix
	
	


    classDef node stroke-width: 4px
    classDef border stroke-width: 4px
    classDef a padding-left: 5em, stroke:#add8e6, stroke-width: 4px
    classDef sub stroke:#add8e6, stroke-width: 4px
```

---

<a name="novelData"></a>
## Novel data

<a name="dataMining"></a>
### Data mining

__Script__: gdc_download.py

__Modules__: requests, json, re, os, tarfile, shutil, gzip

The script is able to search for mutation annotation format (MAF) files for masked somatic mutations from TCGA program in the NCI-GDC database and download all files matching the filters. The filters can also be changed to other programs and data types. The download happens in chunks of 1000 files, so the server does not time out. With this a file `existing_file_ids.txt` is created to keep track of the downloaded files. This file is read within a new run of the download code and will exclude files that are already present.
Afterwards, all the downloads are unpacked until the MAF files and then stored in the `data/maf/` directory and empty directories are deleted for cleanliness.

<a name="mafFiles"></a>
## MAF files available

<a name="dataExtraction"></a>
### Data extraction

__Script__: extractor.py

__Modules__: pandas, os, glob, tqdm

The extraction script loads in the MAF files stored in the `data/maf_files/` directory and puts all files in a list of data frames, which will be concatenated in the end. The combined data will be stored in `data/maf_data.csv`.

---

<a name="conversion"></a>
### Conversion (only for TCGA data)

__Script__: aliquot_to_sample.R

__Libraries__: tidyverse, TCGAutils

With the package ‘TCGAutils’ it is possible to convert barcodes to UUIDs, which is used to obtain the original sample id. This is needed since in Progenetix the sample id is used instead of the aliquot id. During the conversion the variable `Tumor_Sample_Barcode` will be labeled correctly as `aliquot_barcode`, since the given barcode belongs to an aliquot of a sample, and only the first 16 characters are kept as sample barcode. (Hierarchy in GDC Data Portal: Samples > Portions > Analytes > Aliquots). The needed columns for the conversion are:
- Tumor_Sample_UUID
    - GDC aliquot UUID for tumor sample
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
    - Tumor sequencing (discovery) allele 2. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases. ***Tumor_Seq_Allele1 is always the same as the Reference_Allele***
- Tumor_Sample_Barcode
    - Aliquot barcode for the tumor sample

Furthermore, the columns are renamed to match the variable names used in the bycon package variant import script. Lastly, the data is stored in `data/pgx_import.tsv`. The columns kept are:

- case_id
- sample_id
- chromosome
- start
- end
- variant_classification
- snv_type
- reference_bases
- alternate_bases

---

<a name="curation"></a>
### Curation

__Script__: maf_curation_pgx.py

__Modules__: os, pandas, bycon, pymongo, tqdm, numpy

The curation script will load the file created beforehand `pgx_import.tsv` and during the curation process several conventions from Progenetix will be applied:

1. The case and sample ids will get the prefix ‘pgx:TCGA.’
2. ‘chromosome’ is renamed into ‘reference_name’ and the ‘chr’ prefix is discarded
3. The [sequence ontology](http://www.sequenceontology.org/) term for sequence alteration ([SO:0001059](http://www.sequenceontology.org/browser/current_release/term/SO:0001059)) is used to label all variants as SNVs.
   Each SNV gets the sequence ontology term for its specific type:
    - [SO:0001483](http://www.sequenceontology.org/browser/current_release/term/SO:0001483): Single nucleotide variants (SNV)
    - [SO:0002007](http://www.sequenceontology.org/browser/current_release/term/SO:0002007): Multiple nucleotide variants (MNV), which includes DNPs, TNPs, ONPs (≥ 4))
    - [SO:0000159](http://www.sequenceontology.org/browser/current_release/term/SO:0000159): Deletions (DEL)
    - [SO:0000667](http://www.sequenceontology.org/browser/current_release/term/SO:0000667): Insertions (INS)
5. Conversion of the ‘1-start, fully-closed’, also called interresidue-based, into a ‘0-start, half-open’, also called residue-based, genomic coordinate system as recommended by the global alliance for genomics and health ([GA4GH](http://ga4gh.org)) - [https://genomestandards.org/standards/genome-coordinates/](https://genomestandards.org/standards/genome-coordinates/).

    This is achieved by:

    - Subtracting the value 1 of the ‘start’ variable for SNPs, MNPs, and DELs
    - Subtracting the value 1 of the ‘end’ variable for INSs


Additional explanation for the genomic coordinate system: [https://www.biostars.org/p/84686/](https://www.biostars.org/p/84686/)

Then the set of unique sample identifiers is mapped to the Progenetix data base, where the internal positions biosample_id, corresponding to the sample identifier, and individual_id, corresponding to the individual the tumor sample belongs to, is retrieved, if the sample identifier is matched in the data base. Additionally, an internal position callset_id is generated and added for each unique sample identifier. The generated and the retrieved variables are then assigned to the corresponding sample identifier.

In the end, the variants that could not be mapped to a sample identifier will be discarded. The variants ready to be imported will be stored as varImport.tsv in the data directory with the following format:

| biosample_id | variant_id | callset_id | individual_id | reference_name | start | end | reference_sequence | sequence | variant_classification | variant_state_id | specific_so | case_id | sample_id | snv_type |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |

With this the format for the database import is given and the data can be import into the Progenetix MongoDB.

---

<a name="clinvar"></a>
### ClinVar annotation mining

__Script__: clinvar.py

__Modules__: pandas, numpy, json, xml.etree.ElementTree, os, requests, gzip, io.BytesIO, shutil

ClinVar information is being retrieved with the goal to be included in the variant description in the data base.

The ClinVar workflow step uses the script clinvar.py and needs the `maf_data.csv` file to be present.
In a first step variant names are generated from this file in the format `(Hugo_Symbol):HGVSc (HGVSp_Short)` and stored under `data/maf_variant_names.txt`.
If the file is already present the variant names will be loaded instead of generated.

Before mapping the script also checks if mapped variants are already present.

The additional file needed is the [ClinVar variation release](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/):
`ClinVarVariationRelease_00-latest.xml`
This will be downloaded and unpacked if not present.

With the needed prerequisit the XML file is parsed in a 'per element' fashion with `iterparse()`, the starting tag for the element is `VariationArchive` and the resulting element will be searched for the desired entries, according to the [BeaconV2 schema](https://github.com/ga4gh-beacon/beacon-v2/blob/main/models/src/beacon-v2-default-model/genomicVariations/defaultSchema.yaml). In the end the mapped variants are stored in the output file `data/mapped_variants.json`.
An example for the ClinVar XML file can be found [here](https://github.com/ncbi/clinvar/blob/master/sample_xmls/vcv_01.xml).

The resulting objects have the following structure:
- variant_name
- variant_type
- sequence
- reference_sequence
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
    - effect_alternative_ids
- location
    - chromosome
    - start
    - stop
    - sequence_id
    - cytogenetic
