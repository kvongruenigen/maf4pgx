#! /usr/bin/env python3

###############################################################################
# Preparation for variant import to progenetix
# - Rename columns
# - Create variant_id placeholder
# - Change chromosome naming convention
# - Change reference sequence and sequence "-" to "__None__"
# - Change variant type naming convention (SNV, MNV, DEL, INS)
# - Add sequence ontologies
# - Convert coordinate system from 1-based to 0-based
# - Generate callset_ids per aliquot for import
# - Map sample_id to biosample_id and individual_id
# - Create legacy ids / external references

# Outputs: 
# - data/varImport.tsv -> Variants ready for import
# - data/matching_maf_data_curated.csv -> MAF data with matching biosample
###############################################################################

import os
import pandas as pd
from bycon import *
from pymongo import MongoClient
from tqdm import tqdm
import numpy as np

# Connect to MongoDB
client = MongoClient()
db = client.progenetix
bs = db.biosamples

# Read the data frame
maf_data = pd.read_csv("data/pgx_import.tsv", sep = "\t",
    header = 0, low_memory = False)

print("Preparation for mapping...")

# Create placeholder for variant id, gets created while import
maf_data["variant_id"] = [" "] * len(maf_data)

# Rename columns
maf_data.rename(columns={"Variant_Type": "variant_type",
                         "Start_Position": "start",
                         "End_Position": "end",
                         "Reference_Allele": "reference_sequence",
                         "Tumor_Seq_Allele2": "sequence",
                         "Chromosome": "chromosome",}, inplace=True)

# Naming convention from progenetix
maf_data["chromosome"] = maf_data["chromosome"].str.slice(start=3)
maf_data.loc[maf_data["reference_sequence"] == "-", "reference_sequence"] = "__None__"
maf_data.loc[maf_data["sequence"] == "-", "sequence"] = "__None__"



# DNP, TNP, and ONP are registered as MNV - multiple nucleotide variants
maf_data["variant_type"].replace({"SNP": "SNV",
                              "DNP": "MNV",
                              "TNP": "MNV",
                              "ONP": "MNV",
}, inplace=True)

# Adding sequence ontologies - http://www.sequenceontology.org/browser/
maf_data["variant_state_id"] = ["SO:0001059"] * len(maf_data)
maf_data.loc[maf_data["variant_type"] == "SNV", "specific_so"] = "SO:0001483"
maf_data.loc[maf_data["variant_type"] == "MNV", "specific_so"] = "SO:0002007"
maf_data.loc[maf_data["variant_type"] == "DEL", "specific_so"] = "SO:0000159"
maf_data.loc[maf_data["variant_type"] == "INS", "specific_so"] = "SO:0000667"


# Convert 1-based MAF files to 0-based
# Explanation @ https://www.biostars.org/p/84686/
maf_data.loc[maf_data["variant_type"].isin(["SNV", "MNV", "DEL"]), "start"] -= 1
maf_data.loc[maf_data["variant_type"] == "INS", "end"] -= 1

# Generate callset_ids per aliquot for import and map sample_id to biosample_id
print("Starting progenetix mapping...")

biosample_mapping = {}

for sample_id in tqdm(set(maf_data["sample_id"])):
    if sample_id not in biosample_mapping:
        hit = bs.find({"external_references.id": {"$regex": sample_id},
                        "biosample_status.id": "EFO:0009656"})

        for entry in hit:
            biosample_mapping[sample_id] = (entry["id"], entry["individual_id"])

            break

    biosample_id, individual_id = biosample_mapping.get(sample_id, ("", ""))

    cs_id = generate_id("pgxcs")

    maf_data.loc[maf_data["sample_id"] == sample_id,
        ["callset_id", "biosample_id", "individual_id"]] = cs_id, biosample_id, individual_id

# Create legacy ids / external references
maf_data["case_id"] = "pgx:TCGA." + maf_data["case_id"]
maf_data["sample_id"] = "pgx:TCGA." + maf_data["sample_id"]

print("Mapping completed\nWriting files...")

# Replace empty strings with NaN
maf_data = maf_data.replace("", np.nan)

# Drop unnecessary columns
drop = ["Tumor_Seq_Allele1", "Matched_Norm_Sample_Barcode", "Matched_Norm_Sample_UUID",
        "Allele", "normal_bam_uuid", "tumor_bam_uuid",]
maf_data.drop(drop, axis=1, inplace=True)

# Select only matched variants
matching_maf_data = maf_data.dropna(subset = ["biosample_id"])

# Create import file
import_variants = matching_maf_data[[
    "biosample_id", "variant_id", "callset_id", "individual_id",
    "chromosome", "start", "end", "reference_sequence", "sequence",
    "variant_classification", "specific_so", "case_id", "sample_id",
    "variant_type"]]

import_variants.rename(columns={
    "chromosome": "reference_name",
    "specific_so": "variant_state_id",
    }, inplace=True)

# Write finished mapping file
os.makedirs("data/", exist_ok = True) # Check for the directory
import_variants.to_csv("data/varImport.tsv", sep = "\t", index = False)  # and create .tsv file in the directory
matching_maf_data.to_csv("data/matching_maf_data_curated.csv", index = False)  # and create .tsv file in the directory
print("Removed " + str(len(maf_data) - len(matching_maf_data)) + " variants without matching biosample.")

print("Done.\n- Variants ready for import: data/varImport.tsv")