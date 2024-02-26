#! /usr/bin/env python3

###############################################################################
# Insert variant annotation data to existing variants in Progenetix
###############################################################################

from pymongo import MongoClient, UpdateOne
import pandas as pd
from tqdm import tqdm
from pprint import pprint
import ast

# Connect to MongoDB
client = MongoClient()
db = client.progenetix
variants_collection = db.variants

print("Loading data...")
maf_master = pd.read_csv("data/maf_master.csv", low_memory=False, dtype={
    "clinvar_ids": object,
    "proteinHGVS_ids": object,
    "transcriptHGVS_ids": object,
    "variant_alternative_ids": object,
    "aminoacid_changes": object,
    "gene_ids": object,
    "molecular_effects": object,
    "frequency_in_populations": object,
    "clinical_interpretations": object,
})


# Iterate through the MAF data and insert the annotations
print("Creating annotations for bulk write operation...")
bulk_operations = []

for _, variant in tqdm(maf_master.iterrows(), total=maf_master.shape[0]):
    chromosome = variant["chromosome"]
    start = variant["start"]
    stop = variant["end"]
    ref = variant["reference_sequence"]
    alt = variant["sequence"]
    biosample_id = variant["biosample_id"]
    geo = eval(variant["geolocation"])

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
    variant_obj = {}
    variant_obj = {
        "variant_type": {"id": "SO:0001059", "label": "sequence_alteration"},
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
        "geo_location": {
                "type": "Feature",
                "geometry": { 
                    "type": "Point",
                    "coordinates": [
                        geo.get("longitude"), geo.get("latitude")
                        ]
                },
                "properties": {
                    "label": f"{geo.get('city')}, {geo.get('country')}",
                    "city": geo.get("city"),
                    "country": geo.get("country"),
                    "continent": geo.get("continent"),
                    "latitude": geo.get("latitude"),
                    "longitude": geo.get("longitude"),
                    "ISO3166alpha3": geo.get("ISO3166alpha3"),
                    "precision": "city"
                }
            }
    }

    # Convert string representations of lists to actual lists
    for key in ["proteinHGVS_ids", "transcriptHGVS_ids", "variant_alternative_ids",
                "clinvar_ids", "aminoacid_changes", "gene_ids", "molecular_effects"]:
        if key in variant_obj["identifiers"]:
            variant_obj["identifiers"][key] = ast.literal_eval(variant_obj["identifiers"][key])

        if key in variant_obj["molecular_attributes"]:
            variant_obj["molecular_attributes"][key] = ast.literal_eval(variant_obj["molecular_attributes"][key])
    
    if 'clinical_interpretations' in variant_obj and variant_obj['clinical_interpretations'] is not None:
        variant_obj['clinical_interpretations'] = ast.literal_eval(variant_obj['clinical_interpretations'])

    # Remove fields with `None` values
    variant_obj = {k: v for k, v in variant_obj.items() if v is not None}

    # Update the variant
    bulk_operations.append(
        UpdateOne(query, {"$set": variant_obj}, upsert=False)
    )

# Execute bulk write operations
print("Inserting data...")
if bulk_operations:
    result = variants_collection.bulk_write(bulk_operations)
    print("Annotations inserted successfully")
    print("Matched", result.bulk_api_result["nMatched"], "variants")
    print("Modified", result.bulk_api_result["nModified"], "variants")
    pprint(result.bulk_api_result)