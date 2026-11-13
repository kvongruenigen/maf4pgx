##
## Making files
##

# Basic import file
rule import_file:
	input:
		"data/maf_data.csv",
		"data/maf_data_duplicates_removed.csv",
		"data/pgx_import.tsv",
		"data/varImport.tsv",
		"data/varNew.tsv",

# Compact public TCGA export, e.g. data/TCGA-COAD.variants.tsv.gz
rule public_variants:
	input:
		"scripts/public_variants.py"
	output:
		"data/{project}.variants.tsv.gz"
	wildcard_constraints:
		project="TCGA-[A-Z0-9]+"
	shell:
		"python {input:q} {wildcards.project:q} --output {output:q}"

rule public_oncoprint:
	input:
		variants="data/{project}.variants.tsv.gz",
		script="scripts/oncoprint.py"
	output:
		"data/{project}.oncoprint.svg",
		"data/{project}.gene_counts.tsv"
	params:
		prefix="data/{project}"
	shell:
		"python {input.script:q} {input.variants:q} --out-prefix {params.prefix:q}"

rule public_mutation_summary:
	input:
		variants="data/{project}.variants.tsv.gz",
		script="scripts/mutation_summary.py",
		reader="scripts/oncoprint.py"
	output:
		"data/{project}.sample_calls.tsv",
		"data/{project}.recurrent_variants.tsv",
		"data/{project}.gene_pairs.tsv"
	params:
		prefix="data/{project}"
	shell:
		"python {input.script:q} {input.variants:q} --out-prefix {params.prefix:q}"

# All files with annotations
rule make_files:
	input:
		"data/maf_data.csv",
		"data/maf_data_duplicates_removed.csv",
		"data/pgx_import.tsv",
		"data/varImport.tsv",
		"data/varNew.tsv",
		"data/clinvar_variants.json",
		"data/matching_maf_data_curated.csv",
		"data/maf_master.csv"

##
## Imports
##

# For testing purposes
rule test:
	input:
		"data/varImport.tsv"
	shell:
		"~/switchdrive/baudisgroup/dbtools/byconaut/bin/variantsInserter.py --test --i data/varImport.tsv"
# For definite import into the database
rule pgx_import:
	input:
		"data/varImport.tsv"
	shell:
		"~/switchdrive/baudisgroup/dbtools/byconaut/bin/variantsInserter.py --i data/varImport.tsv"

# Importing annotation data into the database
rule annotation_import:
	input:
		"data/maf_master.csv"
	script:
		"scripts/inserter.py"

##
## Novel Data Download
##

# Download the MAF files
rule download:
	script:
		"scripts/gdc_downloader.py"

##
## If MAF files available:
##

# Concatenate the MAF data into a csv file
rule data_extraction:
	output:
		"data/maf_data.csv"
	script:
		"scripts/extractor.py"

# Remove duplicates
rule remove_duplicates:
	input:
		"data/maf_data.csv"
	output:
		"data/maf_data_duplicates_removed.csv"
	script:
		"scripts/duplicates.py"

# Convert the sample barcodes to the sample ids for mapping
rule barcode_conversion:
	input: 
		"data/maf_data_duplicates_removed.csv"
	output: 
		"data/pgx_import.tsv"
	script:
		"scripts/aliquot_to_sample.R"

# Match the sample ids to the database to retrieve biosample ids and individual ids
# Create callset and variant ids
rule mapping:
	input:
		"data/pgx_import.tsv"
	output:
		"data/varImport.tsv",
		"data/varNew.tsv",
		"data/matching_maf_data_curated.csv"
	script:
		"scripts/maf_curation_pgx.py"


# ClinVar annotation
rule clinvar_extraction:
	output:
		"data/clinvar_variants.json"
	script:
		"scripts/clinvar_xml_extractor.py"

# Enhancing MAF data frame
rule maf_annotation:
	input:
		"data/matching_maf_data_curated.csv",
		"data/clinvar_variants.json"
	output:
		"data/maf_master.csv"
	script:
		"scripts/annotator.py"
