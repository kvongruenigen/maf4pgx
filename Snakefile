# This rule executes all the code except the Novel Data Download
rule targets:
	input:
		"data/maf_data.csv",
		"data/pgx_import.tsv",
		"data/varImport.tsv",
		"data/mapped_variants.json"

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
rule extract:
	output:
		"data/maf_data.csv"
	script:
		"scripts/extractor.py"

# Convert the sample barcodes to the sample ids for mapping
rule barcodes:
	input: 
		"data/maf_data.csv"
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
	script:
		"scripts/maf_curation_pgx.py"

# For definite import into the database
rule import:
	input:
		"data/varImport.tsv"
	script:
		"~/switchdrive/baudisgroup/dbtools/byconaut/bin/variantsInserter.py"

# For testing purposes
rule test:
	input:
		"data/varImport.tsv"
	script:
		"~/switchdrive/baudisgroup/dbtools/byconaut/bin/variantsInserter.py --test"

rule clinvar:
	input:
		"data/maf_data.csv"
	output:
		"data/mapped_variants.json"
	script:
		"scripts/clinvar.py"