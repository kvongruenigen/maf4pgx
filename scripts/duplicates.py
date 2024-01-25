import pandas as pd
# Load data
maf_data = pd.read_csv('data/maf_data.csv')

#  Barcode explanation: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
maf_data['sample_barcode'] = maf_data['Tumor_Sample_Barcode'].str[0:16]

#  Analyte codes: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/portion-analyte-codes
maf_data['analyte_code'] = maf_data['Tumor_Sample_Barcode'].str[19]

# Plate numbers:
maf_data['plate_number'] = maf_data['Tumor_Sample_Barcode'].str[21:25]

# Identifier for duplicates based on genomic location, mutation and sample
maf_data['duplicate_identifier'] = (maf_data['Chromosome'].astype(str) + ':'
                                    + maf_data['Start_Position'].astype(str) + '-'
                                    + maf_data['End_Position'].astype(str) + ';'
                                    + maf_data['Reference_Allele'] + '>'
                                    + maf_data['Tumor_Seq_Allele2']
                                    + '(' + maf_data['sample_barcode'] + ')')

number_duplications = maf_data.duplicated(subset=['duplicate_identifier']).sum() # 72'710

duplicated_variants = maf_data[maf_data.duplicated(subset=['duplicate_identifier'], keep = False)]
print(f"Found {number_duplications} duplicated variants")

print("Applying filters to remove duplicates")
print("Filter analyte code: D > W, G, X, unless higher plate number")
variants_to_drop = []
for identifier in duplicated_variants['duplicate_identifier'].unique():
    
    ####  BROAD filtering ####
    # D > W, G, X, unless higher plate number (here higher means higher number)
    df = duplicated_variants.loc[duplicated_variants['duplicate_identifier'] == identifier].sort_values(by=['plate_number', 'analyte_code'], ascending=[False, True])
    # variants_to_keep.append(pd.DataFrame(df.iloc[0]).transpose())
    variants_to_drop.append(pd.DataFrame(df.iloc[1:]))

    # D    63692
    # W     3887

variants_to_drop = pd.concat(variants_to_drop)

indexes_for_dropping = list(variants_to_drop.index)
maf_data.drop(index=indexes_for_dropping)

#  Save data
maf_data.to_csv('data/maf_data_duplicates_removed.csv', index=False)