#! /usr/bin/env python3

###############################################################################
# Concatenate MAF files with ClinVar data
# Create missing values for empty fields
# Add allele frequency and geolocation data
###############################################################################

import pandas as pd
import numpy as np

# Functions

# Function to extract nested data from json
def extract_nested_data(row):
    for column in ["identifiers", "location", "molecular_attributes"]:
        nested_data = row[column]
        if pd.notna(nested_data):
            for key, value in nested_data.items():
                if not value:
                    value = None
                if value is not None:
                    row[key] = value
    return row

# Function to combine allele frequencies under frequency_in_populations
def collect_af(df):
    population_dict = {
        "AA": "African American",
        "AFR": "African",
        "AMR": "American",
        "ASJ": "Ashkenazi Jewish",
        "EA": "European American",
        "EAS": "East Asian",
        "EUR": "European",
        "FIN": "Finnish",
        "MID": "Middle Eastern",
        "NFE": "Non-Finnish European",
        "OTH": "Other",
        "SAS": "South Asian",
    }

    df['frequency_in_populations'] = np.nan
    for i, row in df.iterrows():
        afs = {}
        for column in row.keys():
            if 'AF' in column:
                if 'MAX' in column:
                    continue
                if pd.isna(row[column]) or int(row[column]) == 0:
                    continue
                frequency_name = column.split('_')
                source = frequency_name[0]
                if len(frequency_name) == 2:
                    population = "Overall"
                elif len(frequency_name) == 3:
                    population = population_dict[frequency_name[1]]
                elif len(frequency_name) == 4:
                    population = "Non cancer"
                
                if source not in afs:
                    if source == "gnomAD":
                        afs[source] = {
                                "frequencies": [
                                ],
                                "source": "The Genome Aggregation Database (gnomAD)",
                                "sourceReference": "https://gnomad.broadinstitute.org/",
                            }

                    if source == "1000G":
                        afs[source] = {
                            "frequencies": [
                                ],
                            "source": "1000 Genomes (1000G)",
                            "sourceReference": "https://www.internationalgenome.org/",
                        }

                    if source == "ESP":
                        afs[source] = {
                                "frequencies": [
                                {
                                    "allele_frequency": row[column],
                                    "population": population
                                }
                                ],
                            "source": "NHLBI GO Exome Sequencing Project (ESP)",
                            "sourceReference": "https://esp.gs.washington.edu/drupal/",
                        }

                afs[source]["frequencies"].append({
                    "allele_frequency": row[column],
                    "population": population
                })
        # Convert dictionary to list
        if afs == {}:
            continue

        # Convert dictionary to list
        afs_list = list(afs.values())

        if not afs_list:
            afs_list = np.nan
            continue

        df.at[i, 'frequency_in_populations'] = afs_list
        return df

def geolocation_center(df):
    center_locations = {
        "BI": {"city":"Boston",
            "state": "Massachusetts",
            "country": "USA",
            "continent": "North America",
            "latitude": 42.37,
            "longitude": -71.09},
        "WUGSC": {"city":"St. Louis",
                "state": "Missouri",
                "country": "USA",
                "continent": "North America",
                "latitude": 38.65,
                "longitude": -90.31,},
        "BCM": {"city":"Houston",
                "state": "Texas",
                "country": "USA",
                "continent": "North America",
                "latitude": 29.71,
                "longitude": -95.40},
        "SANGER": {"city":"Cambridge",
                "state": "Cambridgeshire",
                "country": "UK",
                "continent": "Europe",
                "latitude": 52.08,
                "longitude": 0.18},
    }
    df['geolocation'] = np.nan
    for i, row in df.iterrows():
        if pd.notna(row['Center']):
            center = row['Center'].split(';')[0]
            if center in center_locations:
                df.at[i, 'geolocation'] = str(center_locations[center])
    return df

# Load variants matching with data base
print("Loading data...")
maf_matches = pd.read_csv("data/matching_maf_data_curated.csv", low_memory=False)

# Load matched clinvar variants
clinvar_data = pd.read_json('data/clinvar_variants.json', lines=True)

# Drop row with no location data
print("Processing ClinVar data...")
clinvar_data.dropna(subset=['location'], inplace=True)

# Initialize new columns with NaN values
columns_to_extract = ['chromosome', 'start', 'stop', 'sequence_id', 'cytogenetic', 'aminoacid_changes', 'gene_ids',
                       'molecular_effects', 'clinvar_ids', 'transcriptHGVS_ids', 'proteinHGVS_ids', 'genomicHGVS_id',
                       'variant_alternative_ids']
clinvar_data[columns_to_extract] = np.nan

clinvar_data = clinvar_data.apply(extract_nested_data, axis=1)

# Drop nested columns
clinvar_data.drop(columns=["identifiers", "location", "molecular_attributes"], inplace=True)

# Use gene name from variant name
clinvar_data['gene'] = clinvar_data['variant_name'].str.extract(r'\((.*?)\)')

# Start and stop positions as integers
clinvar_data[['start', 'stop']] = clinvar_data[['start', 'stop']].astype(int)

# Change genome coordinate system from 1-based to 0-based
clinvar_data.loc[clinvar_data['variant_type'] == "Insertion", 'stop'] -= 1
clinvar_data.loc[clinvar_data['variant_type'] != "Insertion", 'start'] -= 1
clinvar_data.drop(columns=["variant_type"], axis=1, inplace=True)

# Rename columns
clinvar_data.rename(columns={"stop": "end",
                             "reference_bases": "reference_sequence",
                            "alternate_bases": "sequence",
                            "gene": "SYMBOL",}, inplace=True)


# Save clinvar data
clinvar_data.to_csv("data/clinvar_data.csv", index=False)

# Merge matching variants with clinvar data
print("Merging with MAF data...")
merged_data = pd.merge(maf_matches, clinvar_data, how="left", on=["start", "end", "chromosome", "reference_sequence", "sequence", "SYMBOL"])
merged_data["dbSNP_RS"] = 'dbSNP:' + merged_data["dbSNP_RS"].str.strip('rs')
merged_data["CCDS"] = 'CCDS:' + merged_data["CCDS"].str.strip('CCDS')

# Drop misleading 'VARIANT_CLASS' column
merged_data.drop(columns = ['VARIANT_CLASS'], axis=1, inplace=True)

# Sequence ontology mapping for variant classifications
so_dict = {
    'missense_variant': 'SO:0001583',
    'synonymous_variant': 'SO:0001819',
    'splice_donor_variant': 'SO:0001575',
    'intron_variant': 'SO:0001627',
    'stop_gained': 'SO:0001587',
    'splice_region_variant': 'SO:0001630',
    'splice_acceptor_variant': 'SO:0001574',
    'stop_retained_variant': 'SO:0001567',
    '5_prime_UTR_variant': 'SO:0001623',
    '3_prime_UTR_variant': 'SO:0001624',
    'protein_altering_variant': 'SO:0001818',
    'non_coding_transcript_exon_variant': 'SO:0001792',
    'NMD_transcript_variant': 'SO:0001621',
    'inframe_deletion': 'SO:0001822',
    'frameshift_variant': 'SO:0001589',
    'start_lost': 'SO:0002012',
    'downstream_gene_variant': 'SO:0001632',
    'upstream_gene_variant': 'SO:0001631',
    'inframe_insertion': 'SO:0001821',
    'mature_miRNA_variant': 'SO:0001620',
    'stop_lost': 'SO:0001578',
    'coding_sequence_variant': 'SO:0001580',
    'incomplete_terminal_codon_variant': 'SO:0001626',
    'regulatory_region_variant': 'SO:0001566',
    'non_coding_transcript_variant': 'SO:0001619',
    'start_retained_variant': 'SO:0002019',
}

print("Adding missing information in merged data...")
# Enhance merged data with missing data
for i, row in merged_data.iterrows():

    # Create transcriptHGVS_ids by combining MANE and HGVSc
    if np.all(pd.isna(row['transcriptHGVS_ids'])):
        if pd.notna(row['MANE']) and pd.notna(row['HGVSc']):
            merged_data.at[i, 'transcriptHGVS_ids'] = [str(row['MANE'] + ':' + row['HGVSc'])]
        if pd.notna(row['RefSeq']) and pd.notna(row['HGVSc']):
            merged_data.at[i, 'transcriptHGVS_ids'] = [str(row['RefSeq'].split(';')[0] + ':' + row['HGVSc'])]
        if pd.notna(row['HGVSc']):
            merged_data.at[i, 'transcriptHGVS_ids'] = [str(row['HGVSc'])]


    # Create proteinHGVS_ids by combining ENSP and HGVSp
    if np.all(pd.isna(row['proteinHGVS_ids'])):
        if pd.notna(row['ENSP']) and pd.notna(row['HGVSp']):
            merged_data.at[i, 'proteinHGVS_ids'] = [str(row['ENSP'] + ':' + row['HGVSp'])]
        if pd.notna(row['HGVSp']):
            merged_data.at[i, 'proteinHGVS_ids'] = [str(row['HGVSp'])]


    # Create genomicHGVS_id by combining chromosome, start, reference_sequence, and sequence
    if pd.isna(row['genomicHGVS_id']):
        if pd.notna(row['start']) and pd.notna(row['reference_sequence']) and pd.notna(row['sequence']):
            merged_data.at[i, 'genomicHGVS_id'] = 'g.' + str(row['start']+1) + row['reference_sequence'] + '>' + row['sequence']


    # Create variant_alternative_ids from dbSNP_RS and ClinGen (CCDS)
    if np.all(pd.isna(row['variant_alternative_ids'])):

        if pd.notna(row['dbSNP_RS']) and pd.notna(row['CCDS']):
            merged_data.at[i, 'variant_alternative_ids'] = [row['dbSNP_RS'], row['CCDS']]

        elif pd.notna(row['dbSNP_RS']):
            merged_data.at[i, 'variant_alternative_ids'] = [row['dbSNP_RS']]

        elif pd.notna(row['CCDS']):
            merged_data.at[i, 'variant_alternative_ids'] = [row['CCDS']]

    # If variant_alternative_ids is not empty, append dbSNP_RS and/or CCDS
    else:
        alternative_ids = eval(str(row['variant_alternative_ids']))

        if pd.notna(row['CCDS']) and row['CCDS'] not in row['variant_alternative_ids']:
            alternative_ids.append(row['CCDS'])

        if pd.notna(row['dbSNP_RS']) and row['dbSNP_RS'] not in row['variant_alternative_ids']:
            alternative_ids.append(row['dbSNP_RS'])

        merged_data.at[i, 'variant_alternative_ids'] = alternative_ids


    # Create gene_ids by combining SYMBOL and ENSG
    if np.all(pd.isna(row['gene_ids'])):
        if pd.notna(row['SYMBOL']):
            merged_data.at[i, 'gene_ids'] = [str(row['SYMBOL'])]


    # Create aminoacid_changes by stripping 'p.' from HGVSp_Short
    if np.all(pd.isna(row['aminoacid_changes'])):
        if pd.notna(row['HGVSp_Short']) and not '=' in str(row['HGVSp_Short']):
            merged_data.at[i, 'aminoacid_changes'] = [str(row['HGVSp_Short'].split('.')[1])]


    # Create molecular_effects by combining Variant_Classification and BIOTYPE
    if np.all(pd.isna(row['molecular_effects'])):

        if pd.notna(row['Consequence']):
            consequences = row['Consequence'].split(';')
            merged_data.at[i, 'molecular_effects'] = [{'id': so_dict[consequences[0]], 'label': consequences[0]}]

            if len(consequences) > 1:
                for consequence in consequences[1:]:
                    merged_data.at[i, 'molecular_effects'].append({'id': so_dict[consequence], 'label': consequence})
        
        elif pd.notna(row['One_Consequence']):
            merged_data.at[i, 'molecular_effects'] = [{'id': so_dict[row['One_Consequence']], 'label': row['One_Consequence']}]

    # If molecular_effects is not empty, append Consequence or One_Consequence
    else:
        existing_effects = eval(str(row['molecular_effects']))
        existing_effects_labels = [effect['label'] for effect in existing_effects]
        
        if pd.notna(row['Consequence']):
            for consequence in row['Consequence'].split(';'):
                # make sure consequence is not already present
                if consequence not in existing_effects_labels:
                    existing_effects.append({'id': so_dict[consequence], 'label': consequence})

            merged_data.at[i, 'molecular_effects'] = existing_effects

        elif pd.notna(row['One_Consequence']):
            if row['One_Consequence'] not in existing_effects_labels:
                existing_effects.append({'id': so_dict[row['One_Consequence']], 'label': row['One_Consequence']})

            merged_data.at[i, 'molecular_effects'] = existing_effects


    # Create clinical_interpretations by combining CLIN_SIG, PolyPhen, and SIFT
    if np.all(pd.isna(row['clinical_interpretations'])):
        if pd.notna(row['CLIN_SIG']):
            if pd.notna(row['PolyPhen']) and pd.notna(row['SIFT']):
                merged_data.at[i, 'clinical_interpretations'] = [({
                    'category': {'id': 'MONDO:0000001', 'label': 'disease or disorder'},
                    'condition_id': 'ds0517165',
                    'clinical_relevance': row['CLIN_SIG'],
                    'effect': [{'id': 'PolyPhen', 'label': row['PolyPhen']},
                            {'id': 'SIFT', 'label': row['SIFT']}],
                    'effect_ids': []
                })]

            elif pd.notna(row['PolyPhen']):
                merged_data.at[i, 'clinical_interpretations'] = [({
                    'category': {'id': 'MONDO:0000001', 'label': 'disease or disorder'},
                    'condition_id': 'ds0517165',
                    'clinical_relevance': row['CLIN_SIG'],
                    'effect': {'id': 'PolyPhen', 'label': row['PolyPhen']},
                    'effect_ids': []
                })]

            elif pd.notna(row['SIFT']):
                merged_data.at[i, 'clinical_interpretations'] = [({
                    'category': {'id': 'MONDO:0000001', 'label': 'disease or disorder'},
                    'condition_id': 'ds0517165',
                    'clinical_relevance': row['CLIN_SIG'],
                    'effect': {'id': 'SIFT', 'label': row['SIFT']},
                    'effect_ids': []
                })]

            else:
                merged_data.at[i, 'clinical_interpretations'] = [({
                    'category': {'id': 'MONDO:0000001', 'label': 'disease or disorder'},
                    'condition_id': 'ds0517165',
                    'clinical_relevance': row['CLIN_SIG'],
                    'effect': {'id': 'MedGen:CN517202', 'label': 'not provided'},
                    'effect_ids': ['MedGen:CN517202']
                })]

    # If clinical_interpretations is not empty, append CLIN_SIG
    else:
        interpretations = eval(str(row['clinical_interpretations']))
        interpretations_labels = [interpretation['clinical_relevance'] for interpretation in interpretations]
        if pd.notna(row['CLIN_SIG']) and row['CLIN_SIG'] not in interpretations_labels:
            if pd.notna(row['PolyPhen']) and pd.notna(row['SIFT']):
                interpretations.append({
                    'category': {'id': 'MONDO:0000001', 'label': 'disease or disorder'},
                    'condition_id': 'ds0517165',
                    'clinical_relevance': row['CLIN_SIG'],
                    'effect': [{'id': 'PolyPhen', 'label': row['PolyPhen']},
                            {'id': 'SIFT', 'label': row['SIFT']}],
                    'effect_ids': []
                })

            elif pd.notna(row['PolyPhen']):
                interpretations.append({
                    'category': {'id': 'MONDO:0000001', 'label': 'disease or disorder'},
                    'condition_id': 'ds0517165',
                    'clinical_relevance': row['CLIN_SIG'],
                    'effect': {'id': 'PolyPhen', 'label': row['PolyPhen']},
                    'effect_ids': []
                })

            elif pd.notna(row['SIFT']):
                interpretations.append({
                    'category': {'id': 'MONDO:0000001', 'label': 'disease or disorder'},
                    'condition_id': 'ds0517165',
                    'clinical_relevance': row['CLIN_SIG'],
                    'effect': {'id': 'SIFT', 'label': row['SIFT']},
                    'effect_ids': []
                })

            else:
                interpretations.append({
                    'category': {'id': 'MONDO:0000001', 'label': 'disease or disorder'},
                    'condition_id': 'ds0517165',
                    'clinical_relevance': row['CLIN_SIG'],
                    'effect': {'id': 'MedGen:CN517202', 'label': 'not provided'},
                    'effect_ids': ['MedGen:CN517202']
                })

        merged_data.at[i, 'clinical_interpretations'] = interpretations


# Add geolocation data
print("Adding geolocation data...")
geolocation_center(merged_data)


# Add allele frequency data
print("Adding allele frequency data...")
collect_af(merged_data)

# Save merged data
print("Saving merged data...")
merged_data.to_csv("data/maf_master.csv", index=False)
print("Done!")
merged_data.head()