#!/usr/bin/env Rscript

###############################################################################
# Title: aliquot_to_sample.R
# Aliquot to sample conversion
# This script converts the aliquot barcodes to sample UUIDs
# ONLY FOR TCGA DATA!
###############################################################################

# Load necessary libraries
suppressPackageStartupMessages(library(tidyverse))
library(TCGAutils)
# Import data frame and extract barcodes
cat("Loading data...\n")
data <- read_csv("data/maf_data_duplicates_removed.csv", show_col_types = FALSE)

if ("Tumor_Sample_UUID" %in% colnames(data) == TRUE) {

  # If barcode longer than 16 characters then it is an aliquot barcode
  if (nchar(data$Tumor_Sample_Barcode[1]) > 16) {
    data$aliquot_barcode <- data$Tumor_Sample_Barcode

    sample_barcodes <- unique(data["sample_barcode"])

    # Create an empty list for the ids
    sample_id <- list()

    # Convert barcodes into ids
    cat("Converting barcode to sample UUID...\n")
    for (id in sample_barcodes){
      sam <- barcodeToUUID(id)
      sam <- sam$sample_id
      sample_id <- c(sample_id, sam)
    }

    # Make a data frame for mapping
    mapping_df <- data.frame(
      unlist(as.list(sample_barcodes)),
      unlist(sample_id)
    )
    colnames(mapping_df) <- c("sample_barcode", "sample_id")

    # Join the two data frames based on matching Barcodes
    mapfile <- left_join(
      data, mapping_df,
      by = c("sample_barcode" = "sample_barcode")
    )

    cat("Converting completed.\n")

    # Renaming for further processing
    #####################################################################

    # Rename the columns
    mapfile <- mapfile %>%
      rename(
        aliquot_id = Tumor_Sample_UUID,
        chromosome = Chromosome,
        start = Start_Position,
        end = End_Position,
        variant_classification = Variant_Classification,
        variant_type = Variant_Type,
        reference_sequence = Reference_Allele,
        sequence = Tumor_Seq_Allele2,
      )

    cat("Writing output file...\n")
    # Write file
    write_tsv(mapfile, "data/pgx_import.tsv")

    cat("Done\n")
  } else {
    cat("No conversion necessary.\n")
    mapfile <- data
    # Rename the columns
    mapfile <- mapfile %>%
      rename(
        aliquot_id = Tumor_Sample_UUID,
        chromosome = Chromosome,
        start = Start_Position,
        end = End_Position,
        variant_classification = Variant_Classification,
        variant_type = Variant_Type,
        reference_sequence = Reference_Allele,
        sequence = Tumor_Seq_Allele2,
        sample_barcode = Tumor_Sample_Barcode,
      )
    cat("Writing output file...\n")
    # Write file
    write_tsv(mapfile, "data/pgx_import.tsv")
  }
} else {
  cat("No Tumor_Sample_UUID found. Please check your input file.\n")
}