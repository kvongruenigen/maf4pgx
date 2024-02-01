#!/usr/bin/env python3

###############################################################################
# Script to extract ClinVar variant data from XML to JSON
###############################################################################

# Import packages #############################################################
import pandas as pd
import numpy as np
import json
import xml.etree.ElementTree as ET
import os
import requests
import gzip
from io import BytesIO
import shutil
# import pprint

# Functions ###################################################################
# ID generator from Rahel #####################################################
def conditionid(id):
    if len(id) == 7:
        conditionID = "ds" + id
    else:
        length = 7 - len(id)
        prefix = length * "0"
        conditionId = "ds" + prefix + id
        return conditionId

# Function to extract clinical information for a given VariationArchive element
def process_variation(variation_archive):
    variation_name = variation_archive.get("VariationName")

    # Needed format: (Gene):HGVSc (HGVSp)
    if "(" in variation_name:
        # Remove transcript information from the variation name
        short_variation_name = variation_name.split("(")[1:]
        # Rejoin with brackets
        short_variation_name = "(" + "(".join(short_variation_name)

    # If there is no bracket in the variation name, it is the wrong format
    else:
        return None
    
    # Run mapping
    
    # Create variant dictionary
    variant_dictionary = {}

    ###########################################################################
    # Location
    ###########################################################################

    # Location ----------------------------------------------------------------
    for location in variation_archive.findall(".//InterpretedRecord/SimpleAllele/Location/SequenceLocation"):
        if location.get("Assembly") == "GRCh38":
            variant_dictionary["location"] = {
                "chromosome": location.get("Chr"),
                "start": location.get("start"),
                "stop": location.get("stop"),
                "sequence_id": location.get("Accession"),
            }

            # Cytogenetic location --------------------------------------------
            cytogenetic_location = variation_archive.find(".//InterpretedRecord/SimpleAllele/Location/CytogeneticLocation")
            if cytogenetic_location != None and cytogenetic_location.text != None:
                variant_dictionary["location"]["cytogenetic"] = cytogenetic_location.text
        
            # Reference and alternate bases -----------------------------------
            variant_dictionary["alternate_bases"] = location.get("alternateAlleleVCF")
            variant_dictionary["reference_bases"] = location.get("referenceAlleleVCF")
    if variant_dictionary.get("location") == None:
        return None
    else:
        chromosome = variant_dictionary.get("location").get("chromosome", None)
        start = variant_dictionary.get("location").get("start", None)
        stop = variant_dictionary.get("location").get("stop", None)
        if chromosome == None or start == None or stop == None:
            return None
    ###########################################################################
    # Variant name and type
    variant_dictionary["variant_name"] = variation_name
    variant_type = variation_archive.get("VariationType")
    variant_dictionary["variant_type"] = variant_type


    ###########################################################################
    # Identifiers
    ###########################################################################

    # ClinVar IDs -------------------------------------------------------------
    variant_dictionary["identifiers"] = {
        "clinvar_ids": [variation_archive.get("Accession"),
                        variation_archive.get("VariationID")],
        "transcriptHGVS_ids": [],
        "proteinHGVS_ids": [],
        "variant_alternative_ids": [],
    }

    ## HGVS -------------------------------------------------------------------
    for HGVS in variation_archive.findall(".//InterpretedRecord/SimpleAllele/HGVSlist/HGVS"):
        # Exclude empty HGVS and GRCh37 variants
        if HGVS.get("Type") != None or HGVS.get("Assembly") != "GRCh37":

            # Genomic HGVS ids ------------------------------------------------
            if "genomic, top-level" in HGVS.get("Type"):
                genomic_hgvs = HGVS.find(".//NucleotideExpression/Expression")
                if genomic_hgvs != None:
                    variant_dictionary["identifiers"]["genomicHGVS_id"] = genomic_hgvs.text

            if "coding" in HGVS.get("Type"):
                # Transcript HGVS ids -----------------------------------------
                transcript_hgvs = HGVS.find(".//NucleotideExpression/Expression")
                if transcript_hgvs != None:
                    variant_dictionary["identifiers"]["transcriptHGVS_ids"].append(transcript_hgvs.text)
                
                # Protein HGVS ids --------------------------------------------
                protein_hgvs = HGVS.find(".//ProteinExpression/Expression")
                if protein_hgvs != None:
                    variant_dictionary["identifiers"]["proteinHGVS_ids"].append(protein_hgvs.text)

    # Alternative ids ---------------------------------------------------------
    for ref in variation_archive.findall(".//InterpretedRecord/SimpleAllele/XRefList/XRef"):
        if ref.attrib != {}:
            variant_dictionary["identifiers"]["variant_alternative_ids"].append(str(ref.get("DB")) + ":" + ref.get("ID"))

    ###########################################################################
    # Molecular attributes
    ###########################################################################

    variant_dictionary["molecular_attributes"] = {
        "gene_ids": [],
        "aminoacid_changes": [],
        "molecular_effects": [],
    }

    # gene_ids ----------------------------------------------------------------
    for gene in variation_archive.findall(".//InterpretedRecord/SimpleAllele/GeneList/Gene"):
        variant_dictionary["molecular_attributes"]["gene_ids"].append(gene.get("Symbol"))

    # aminoacid_changes -------------------------------------------------------
    for aa_change in variation_archive.findall(".//InterpretedRecord/SimpleAllele/ProteinChange"):
        variant_dictionary["molecular_attributes"]["aminoacid_changes"].append(aa_change.text)

    # molecular_effects -------------------------------------------------------
    for molecular_consequence in variation_archive.findall(".//InterpretedRecord/SimpleAllele/HGVSlist/HGVS/MolecularConsequence"):
        if molecular_consequence.attrib != {}:
            molecular_effect = {
                "id": molecular_consequence.get("ID"),
                "label": molecular_consequence.get("Type")
            }
            if molecular_effect not in variant_dictionary["molecular_attributes"]["molecular_effects"]:
                variant_dictionary["molecular_attributes"]["molecular_effects"].append(molecular_effect)

    ###########################################################################
    # Clinical interpretations (Variant level data)
    ###########################################################################

    variant_dictionary["clinical_interpretations"] = []

    # Get all traits
    for trait in variation_archive.findall(".//Interpretation/ConditionList/TraitSet/Trait"):
        if trait.get("Type") != None:
            preferred_label = None
            labels = []
            main_id = None

            # Find effect and label
            for name in trait.findall("Name"):
                if name.find("ElementValue") != None:
                    if name.find("ElementValue").get("Type") == "Preferred":
                        preferred_label = str(name.find("ElementValue").text)

                if name.find("XRef") != None:
                    if name.find("XRef").get("DB") == "MONDO":
                        main_id = str(name.find("XRef").get("ID"))

            # Find all labels
            for xref in trait.findall("XRef"):
                if xref.get("DB") == "MONDO":
                    main_id = str(xref.get("ID"))
                elif ":" in xref.get("ID"):
                    labels.append(str(xref.get("ID")))
                else:
                    labels.append(str(xref.get("DB")) + ":" + str(xref.get("ID")))

            # Hierarchy of preferred ids: MONDO > HPO > OMIM > first label
            if main_id != None:
                labels.append(main_id)
            else:
                for label in labels:
                    if "HP" in label:
                        main_id = label
                        break
                    elif "OMIM" in label:
                        main_id = label
                    else:
                        main_id = labels[0]

            # Create clinical interpretation dictionary and append to variant dictionary
            clinical_interpretation_dictionary = {
                "category": { "id": "MONDO:0000001", "label": "disease or disorder" },
                "condition_id": conditionid(variation_archive.get("VariationID")),
                "clinical_relevance": variation_archive.find(".//Interpretation/Description").text,
                "effect": {
                    "id": main_id,
                    "label": preferred_label
                },
                "effect_ids": labels,
            }
            variant_dictionary["clinical_interpretations"].append(clinical_interpretation_dictionary)

    # pprint.pprint(variant_dictionary)
    return variant_dictionary


# Loading variant data for ClinVar querying ###################################

# Load XML data with iterparse ################################################

# Check if the file exists ----------------------------------------------------


# Output directory
output_directory = "data/"
file_path = os.path.join(output_directory, "ClinVarVariationRelease_00-latest.xml")

if not os.path.isfile(file_path):
    print("Downloading file...")
    # Link to ClinVar XML data
    url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/ClinVarVariationRelease_00-latest.xml.gz"
    response = requests.get(url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        print("Unzipping file...")
        with gzip.open(BytesIO(response.content), 'rt', encoding='utf-8') as gz_file:
            with open(file_path, 'w', encoding='utf-8') as output_file:
                shutil.copyfileobj(gz_file, output_file)
        print("File downloaded and unzipped successfully.")

    else:
        print(f"Error downloading file. Status code: {response.status_code}")

# If the file exists, start mapping -------------------------------------------
print("Start mapping...")

variation_xml = "data/ClinVarVariationRelease_00-latest.xml"

for event, elem in ET.iterparse(variation_xml, events=("start", "end")):

    # VariationArchive is root tag
    if event == "start" and elem.tag == "VariationArchive":

        # Process the VariationArchive element
        found_variant = process_variation(elem)

        # Clear the element to free up memory
        elem.clear()

        # If a variant was found, append it to the list file
        if found_variant != None:
            with open("data/clinvar_variants.json", "a") as f:
                f.write(json.dumps(found_variant) + "\n")