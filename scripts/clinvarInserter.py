#! /usr/bin/env python3
from pymongo import MongoClient
import json

# Connect to MongoDB
client = MongoClient()
db = client.progenetix
variantsCollection = db.variants

with open('data/mapped_variants.json') as f:
    mapped_variants = [json.loads(line) for line in f]

for variant in mapped_variants:
    location = variant.get("location")
    if location == None:
        continue
    else:
        chromosome = variant.get("location").get("chromosome")
        start = variant.get("location").get("start")
        stop = variant.get("location").get("stop")
        ref = variant.get("reference_bases")
        alt = variant.get("alternate_bases")

        if chromosome == None or start == None or stop == None or ref == None or alt == None:
            continue

    variant_type = variant.get("variant_type") 
    start = int(start)
    stop = int(stop)
    # {'Indel', 'Duplication', 'single nucleotide variant', 'Insertion', 'Microsatellite', 'Deletion'}
    if variant_type == "Insertion":
        stop -= 1
    else:
        start -= 1

    query = {"location.chromosome": chromosome,
             "location.start": start,
             "location.end": stop,
             "reference_sequence": ref,
             "sequence": alt}

    # # Removing fields that are already in the database
    variant.pop("location")
    variant.pop("alternate_bases")
    variant.pop("reference_bases")
    variant.pop("variant_type")

    # Update the variant
    variantsCollection.update_many(query, {"$set": variant}, upsert=False)