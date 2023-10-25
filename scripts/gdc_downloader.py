###############################################################################

# Downloading masked somatic mutation files from the TCGA program of GDC
# Results in the download of all unrestricted MAF files 

###############################################################################
import requests
import json
import re
import os
import tarfile
import shutil
import gzip
import pandas as pd
import glob
from tqdm import tqdm

# Access the file endpoint from GDC for id retrieval
files_endpt = "https://api.gdc.cancer.gov/files"


# Filtering for TCGA Masked Somatic Mutation
# - results in open access maf files.
# This set of filters is nested under an "and" operator.
filters = {
    "op": "and",
    "content": [
        {
            "op": "in",
            "content": {
                "field": "cases.project.program.name",
                "value": [input("Enter program name (e.g., TCGA): ")]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "files.data_type",
                "value": [input("Enter data type (e.g., Masked Somatic Mutation): ")]
            }
        }
    ]
}

# Here a GET is used, so the filter parameters should be passed as a
# JSON string.
params = {
    "filters": json.dumps(filters),
    "fields": "file_id",
    "format": "JSON",
    "size": "20000"
}

# Download the ids
response = requests.get(files_endpt, params = params)

# Create a list for the ids
file_uuid_list = []

# This step populates the download list with the file_ids from the 
# previous query
for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
    file_uuid_list.append(file_entry["file_id"])
print("Found ", len(file_uuid_list), "files to download.")

# Import record of existing files
existing_ids = []
if os.path.isfile("data/existing_file_ids.txt"):
    with open("data/existing_file_ids.txt") as f:
        lines = f.readlines()
    for i in lines:
        existing_ids.append(i)

# Remove existing file UUIDs from download queue
removed_ids = []
for i in existing_ids:
    i = i.strip("\n")
    if i in file_uuid_list:
        file_uuid_list.remove(i)
        removed_ids.append(i)
print("Removed", len(removed_ids), "existing ids from download queue.")

# Create an empty list to store UUID chunks
ls = []

# Append sublists to "ls" for 1000 ids each (server limits)
for i in range(0, len(file_uuid_list), 1000):
    ls.append(file_uuid_list[i:i + 1000])

# Download the files
print("Starting", len(ls),"downloads...")
downloaded = 1
for idls in ls:

    # Exclude existing files
    for ids in idls:
        if ids in existing_ids:
            idls.remove(ids)
            print("Removed existing file id:", ids)
    
    data_endpt = "https://api.gdc.cancer.gov/data"

    params = {"ids": idls}

    response = requests.post(data_endpt,
                            data = json.dumps(params),
                            headers = {"Content-Type": "application/json"})

    response_head_cd = response.headers["Content-Disposition"]

    file_name = re.findall("filename=(.+)", response_head_cd)[0]

     # Check for the directory
    os.makedirs("data/maf/", exist_ok = True)
    save_path = "data/maf/"

    complete_name = os.path.join(save_path, file_name)

    with open(complete_name, "wb") as output_file:
        output_file.write(response.content)

    if len(ls) - downloaded != 0:
        print("Downloaded completed. Remaining:", len(ls) - downloaded, "of", len(ls))
        downloaded += 1
    else:
        print("Download finished")

# Create record if imported file UUIDs
for i in file_uuid_list:
    i = i + "\n"
    existing_ids.append(i)

# Save existing file UUIDs
os.makedirs("data/", exist_ok = True)
with open("data/existing_file_ids.txt", "w") as file_ids:
    for i in existing_ids:
        file_ids.write(i)

# Unzip all files
def unpack():
    path = os.getcwd()
    # Move to directory with the downloaded data
    os.chdir(path + '/data/maf')

    # Look for nested tar.gz files and extract them
    for root, _, files in os.walk('.'):
        for file in files:
            if file.endswith('.tar.gz'):
                with tarfile.open(os.path.join(root, file), 'r:gz') as tar:
                    tar.extractall(root)

    # Delete all tar.gz files
    for root, _, files in os.walk('.'):
        for file in files:
            if file.endswith('.tar.gz'):
                os.remove(os.path.join(root, file))

    # Find all directories in the current directory
    for dir_name in os.listdir('.'):
        if os.path.isdir(dir_name):
            # Go into the directory
            os.chdir(dir_name)

            # Look for tar.gz files and extract them
            for file in os.listdir('.'):
                if file.endswith('.maf.gz'):
                    with open(file[:-3], 'wb') as f_out, gzip.open(file, 'rb') as f_in:
                        shutil.copyfileobj(f_in, f_out)

                    os.remove(file)

            # Move all files in the directory up one level
            for file in os.listdir('.'):
                shutil.move(file, '..')

            # Go back up to the parent directory
            os.chdir('..')

            # Remove empty folder
            os.rmdir(dir_name)

    # Go back to the original directory
    os.chdir(path)

print("Unpacking files...")
unpack()