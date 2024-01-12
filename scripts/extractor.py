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
combined_data = pd.DataFrame()
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
        maf_data = pd.read_csv(file, sep="\t", skiprows=lambda x: x < 7 or lines[x].startswith("#"), header=0, low_memory=False)
    else:
        # No comments, directly read the file
        maf_data = pd.read_csv(file, sep="\t", header=0, low_memory=False)
        
    df_list.append(maf_data)

# Create one data frame from the list
combined_data = pd.concat(df_list).reset_index(drop=True)
print("Data extraction completed.")

# Removing empty columns #####################################################

print("Removing empty columns...")

original_columns = combined_data.columns

nan_value = float("NaN")

combined_data.replace("", nan_value, inplace=True)

combined_data.dropna(how='all', axis=1, inplace=True)

columns_with_values = combined_data.columns

dropped_columns = set(original_columns) - set(columns_with_values)

print('Columns removed:')
for column in dropped_columns:
    print(column)

# Remove duplicates ##########################################################
print("Checking for duplicated variants...")

# Generate needed columns

#  Barcode explanation: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
combined_data['sample_barcode'] = combined_data['Tumor_Sample_Barcode'].str[0:16]

#  Analyte codes: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/portion-analyte-codes
combined_data['analyte_code'] = combined_data['Tumor_Sample_Barcode'].str[19]

# Plate numbers:
combined_data['plate_number'] = combined_data['Tumor_Sample_Barcode'].str[21:25]

# Identifier for duplicates based on genomic location, mutation and sample
combined_data['duplicate_identifier'] = (combined_data['Chromosome'].astype(str) + ':'
                                    + combined_data['Start_Position'].astype(str) + '-'
                                    + combined_data['End_Position'].astype(str) + ';'
                                    + combined_data['Reference_Allele'] + '>'
                                    + combined_data['Tumor_Seq_Allele2']
                                    + '(' + combined_data['sample_barcode'] + ')')

number_duplications = combined_data.duplicated(subset=['duplicate_identifier']).sum()

duplicated_variants = combined_data[combined_data.duplicated(subset=['duplicate_identifier'], keep = False)]
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
    combined_data.drop(index=indexes_for_dropping)
    combined_data.drop(columns=['sample_barcode', 'analyte_code', 'plate_number', 'duplicate_identifier'])
    print("Done.")


# Check for the directory
os.makedirs("data/", exist_ok=True)
print("Writing output file...")
# and create .csv file in the directory
combined_data.to_csv("data/maf_data.csv", index=False)
print("Done.")
