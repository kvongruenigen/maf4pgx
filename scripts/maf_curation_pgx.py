# Preparation for variant import to progenetix

#####################################################################

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
maf_data.columns = maf_data.columns.str.lower()

print("Preparation for mapping...")

# Create placeholder for variant id, gets created while import
maf_data["variant_id"] = [" "] * len(maf_data)

# Naming convention from progenetix
maf_data["reference_name"] = maf_data["chromosome"].str.slice(start=3)
maf_data.loc[maf_data["reference_sequence"] == "-", "reference_sequence"] = "None"
maf_data.loc[maf_data["sequence"] == "-", "sequence"] = "None"



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

maf_data = maf_data.replace("", np.nan)

matching_maf_data = maf_data.dropna(subset = ["biosample_id"])
# Clean up
import_variants = matching_maf_data[["biosample_id", "variant_id", "callset_id", "individual_id",
    "reference_name", "start", "end", "reference_sequence",
    "sequence", "variant_classification", "variant_state_id",
    "specific_so", "case_id", "sample_id", "variant_type"]]

# Write finished mapping file
os.makedirs("data/", exist_ok = True) # Check for the directory
import_variants.to_csv("data/varImport.tsv", sep = "\t", index = False)  # and create .tsv file in the directory
matching_maf_data.to_csv("data/matching_maf_data_curated.csv")  # and create .tsv file in the directory

print("Done.\n- Variants ready for import: data/varImport.tsv")