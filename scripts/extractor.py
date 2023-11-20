# Data extraction from MAF files
# All information will be stored at data/maf_data.csv

##############################################################################
import os
import pandas as pd
import numpy as np
import glob
from tqdm import tqdm

# Extract data from MAF files ################################################
# Create a dataframe and a list for the MAF data
data = pd.DataFrame()
df_list = []

# Iterate through directory with downloaded maf files and
# load the relevant information in "data"
print("Starting data extraction...")

for file in tqdm(glob.glob("data/maf/*.maf"), desc="Extraction progress"):

    # Read the file and check the first few lines for comments
    with open(file, "r") as f:
        lines = f.readlines()
    
    # Check for the presence of comments in the first 10 lines
    has_comments = any(line.startswith("#") for line in lines[:10])
    
    if has_comments:
        # Skip commented rows
        combined_data = pd.read_csv(file, sep="\t", skiprows=lambda x: x < 7 or lines[x].startswith("#"), header=0, low_memory=False)
    else:
        # No comments, directly read the file
        combined_data = pd.read_csv(file, sep="\t", header=0, low_memory=False)
        
    df_list.append(combined_data)

# Create one data frame from the list
data = pd.concat(df_list).reset_index(drop=True)
print("Data extraction completed.")

# Removing empty columns #####################################################

print("Removing empty columns...")

original_columns = data.columns

nan_value = float("NaN")

data.replace("", nan_value, inplace=True)

data.dropna(how='all', axis=1, inplace=True)

columns_with_values = data.columns

dropped_columns = set(original_columns) - set(columns_with_values)

print('Columns removed:')
for column in dropped_columns:
    print(column)

# Remove duplicates ##########################################################
print("Checking for duplicated variants...")

# Generate needed columns

#  Barcode explanation: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
data['sample_barcode'] = data['Tumor_Sample_Barcode'].str[0:16]

#  Analyte codes: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/portion-analyte-codes
data['analyte_code'] = data['Tumor_Sample_Barcode'].str[19]

# Plate numbers:
data['plate_number'] = data['Tumor_Sample_Barcode'].str[21:25]

# Identifier for duplicates based on genomic location, mutation and sample
data['duplicate_identifier'] = (data['Chromosome'].astype(str) + ':'
                                    + data['Start_Position'].astype(str) + '-'
                                    + data['End_Position'].astype(str) + ';'
                                    + data['Reference_Allele'] + '>'
                                    + data['Tumor_Seq_Allele2']
                                    + '(' + data['sample_barcode'] + ')')

number_duplications = data.duplicated(subset=['duplicate_identifier']).sum()

duplicated_variants = data[data.duplicated(subset=['duplicate_identifier'], keep = False)]
print(f"Found {number_duplications} duplicated variants in the dataset.")

if number_duplications > 0:
    print("Removing duplicated variants...")

    # Create list for the data frame with the duplicated variants
    variants_to_drop = []

    # Loop through each duplication, so that sorting works
    for identifier in duplicated_variants['duplicate_identifier'].unique():

        # D > W, G, X, unless higher plate number (here higher means higher number)
        df = duplicated_variants.loc[duplicated_variants['duplicate_identifier'] == identifier].sort_values(by=['plate_number', 'analyte_code'], ascending=[False, True])

        # df.iloc[0] will be kept, the rest dropped
        variants_to_drop.append(pd.DataFrame(df.iloc[1:]))

    # Concatenate the dataframes to one
    variants_to_drop = pd.concat(variants_to_drop)

    # Get the indexes for the variants to drop
    indexes_for_dropping = list(variants_to_drop.index)
    data.drop(index=indexes_for_dropping)
    data.drop(columns=['sample_barcode', 'analyte_code', 'plate_number', 'duplicate_identifier'])
    print("Done.")

# Add anchor base for DEL and INS
data.loc[data['Variant_Type'].isin(['DEL', 'INS']), 'vcf_anchor_base'] = data['CONTEXT'].str[4]

# Check for the directory
os.makedirs("data/", exist_ok=True)
print("Writing output file...")
# and create .csv file in the directory
data.to_csv("data/maf_data.csv", index=False)
print("Done.")
