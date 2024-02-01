#! /usr/bin/env python3

###############################################################################
# Insert variant annotation data to existing variants in Progenetix
###############################################################################

from pymongo import MongoClient
import pandas as pd

# Connect to MongoDB
client = MongoClient()
db = client.progenetix
variantsCollection = db.variants

maf_master = pd.read_csv("../data/maf_master.csv")

for _, variant in maf_master.iterrows():
    chromosome = variant["chromosome"]
    start = variant["start"]
    stop = variant["end"]
    ref = variant["reference_sequence"]
    alt = variant["sequence"]
    biosample_id = variant["biosample_id"]

    if any(pd.isna(value) for value in [chromosome, start, stop, ref, alt]):
        continue

    query = {
          "location.chromosome": chromosome,
          "location.start": start,
          "location.end": stop,
          "reference_sequence": ref,
          "sequence": alt,
          "biosample_id": biosample_id,
          }

    # Build the variant object
    variant_obj = {
        "variant_type": {"id": variant["specific_so"], "label": variant["variant_type"]},
        "identifiers": {
            key: value for key, value in variant[
                  ["clinvar_ids",
                   "genomicHGVS_id",
                   "proteinHGVS_ids",
                   "transcriptHGVS_ids",
                   "variant_alternative_ids"]
                   ].to_dict().items()
            if value is not None and pd.notna(value)
        },
        "molecular_attributes": {
            key: value for key, value in variant[
                ["aminoacid_changes",
                 "gene_ids",
                 "molecular_effects"]
                 ].to_dict().items()
            if value is not None and pd.notna(value)
        },
        "variant_name": variant["variant_name"] if not pd.isna(variant["variant_name"]) else None,
        "clinical_interpretations": variant["clinical_interpretations"] if pd.notna(variant["clinical_interpretations"]) else None,
        "frequency_in_populations": variant["frequency_in_populations"] if pd.notna(variant["frequency_in_populations"]) else None,
        "allele_origin": {"id": "SO:0001777","label": "somatic variant"},
    }

    # Remove fields with `None` values
    variant_obj = {k: v for k, v in variant_obj.items() if v is not None}

    # Update the variant
    variantsCollection.update_many(query, {"$set": variant_obj}, upsert=False)