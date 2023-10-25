# Data extraction from MAF files
# All information will be stored at data/maf_data.csv

#####################################################################
import os
import pandas as pd
import glob
from tqdm import tqdm

path = input("Enter path to project directory (default = current directory): ") # "/Users/kayvongrunigen/Projects/snpettes/final-workflow"
if path == "":
    path = os.getcwd()

# Set working directory
os.chdir(path)

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

# Check for the directory
os.makedirs("data/", exist_ok=True)
print("Writing output file...")
# and create .csv file in the directory
data.to_csv("data/maf_data.csv", index=False)
print("Done.")
