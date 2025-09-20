# Supporting Code: 'Diverse defence systems and prophages in human-associated Bifidobacterium species reveal coevolutionary "arms race" dynamics'.

Created by: James Docherty

# Scripts for Bioinformatics Analysis

Defence Systems & Prophages in Human-associated Bifidobacterium strains

# 01 | Selecting for human samples after downloading genomes from NCBI and collating metadata from assembly data reports

## # Download all genomes from NCBI and accompanying assembly data report then follow these scripts step by step to get genomes from human isolates / MAGs

### **Script for selecting strains isolated from humans. "*Homo sapiens*" is used to select strains from the Assembly Data Report (Python Script)**

```python
# To save output to file, we need to import sys
import sys

# When you run job, it will show 'Enter the text: ' in the prompt
user_input =input("Enter the text: ")

# File of interest
file =open("assembly_data_report.jsonl", 'r')

# Program script | "bif_human_list.txt" is the output file and "a" means to create new file without overwriting anything
for line in file:
    if user_input in line:
        print(line)
        stdoutOrigin = sys.stdout
        sys.stdout = open("bif_human_list.txt", "a")
```

### **In order to transfer all relevant folders containing genomes, we must first use this script to remove all of the text except the accession numbers in the assembly data report (original file is left untouched)**

```bash
cat bif_human_list.txt | awk -F "," '{print $1}' | sed 's/{"accession":"//g' | sed 's/"//g' | sort | uniq > bif_human_list_accessions.txt
```

### **Move relevant genome folders to dedicated directory**

```bash
# Make dedicated directory
mkdir selected_genomes

# Move relevant genome folders to dedicated directory
for X in `cat bif_human_list_accessions.txt`; do mv $X selected_genomes/${SX} ; done

# Move all selected genomes with a certain extension from multiple subdirectories into one directory
mv **/*.fna /path/to/selected_genomes

# Remove all not-needed sub-directories from within a directory
rm -R -- ./*
```

### **To put species name next to accession numbers in filenames (in a separate txt file)**

```python
# To get all the accession numbers paired with the species name
import json

filename = "bif_human_list.txt"

def get_accession_species(json_str):
    data = json.loads(json_str)
    accession = data["accession"]
    species = data["organism"]["organismName"]
    return accession, species

with open(filename, "r") as file:
    for line in file:
        accession, species = get_accession_species(line)
        print("Accession:", accession)
        print("Species:", species)
        print("--------------------")
```

## # Step 2: Extract metadata from assembly data reports

### Extract attributes from assembly data report

```python
import json

def parse_attributes(json_data):
    try:
        attributes = json_data.get('assemblyInfo', {}).get('biosample', {}).get('attributes', [])
        attribute_names = ["strain", "host", "isolation_source", "collection_date", "geo_loc_name", "sample_type"]
        parsed_attributes = {}

        for attribute in attributes:
            name = attribute.get('name', '')
            value = attribute.get('value', '')
            if name in attribute_names:
                parsed_attributes[name] = value

        return parsed_attributes

    except json.JSONDecodeError:
        print("Invalid JSON format.")
        return None

if __name__ == "__main__":
    file_path = "bif_human_list_for_metadata.txt"
    output_file = "output_attributes.tsv"

    try:
        with open(file_path, 'r') as infile, open(output_file, 'w') as outfile:
            outfile.write("accession\torganismName\tstrain\thost\tisolation_source\tcollection_date\tgeo_loc_name\tsample_type\n")
            
            for line in infile:
                try:
                    json_data = json.loads(line)
                    accession = json_data.get('accession', '')
                    organism_name = json_data.get('assemblyInfo', {}).get('biosample', {}).get('description', {}).get('organism', {}).get('organismName', '')
                    attributes_data = parse_attributes(json_data)

                    if attributes_data:
                        # Write the attributes as tab-separated values to the output file
                        outfile.write(f"{accession}\t{organism_name}\t")
                        outfile.write("\t".join(attributes_data.get(attr, "") for attr in ["strain", "host", "isolation_source", "collection_date", "geo_loc_name", "sample_type"]) + "\n")
                except json.JSONDecodeError:
                    print("Invalid JSON format in line:", line)

        print("Output has been written to", output_file)

    except FileNotFoundError:
        print(f"File not found: {file_path}")
```

### Combine metadata with output statistics from CheckM (see Section 02)

```python
import pandas as pd

# Read the input files into pandas dataframes
output_df = pd.read_csv('CheckM_statistics.tsv', sep='\t')
attributes_df = pd.read_csv('output_attributes_from_metadata.tsv', sep='\t')

# Merge the dataframes based on the "accession" column
merged_df = pd.merge(output_df, attributes_df, left_on='Strain', right_on='accession', how='inner')

# Drop the redundant "accession" column from the attributes dataframe
merged_df.drop('accession', axis=1, inplace=True)

# Write the merged dataframe to a new file
merged_df.to_csv('merged_output.tsv', sep='\t', index=False)
```

### Reorder the combined output to distinguish attributes from CheckM output statistics (see Section 02)

```python
import pandas as pd

# Read the input files into pandas dataframes
output_df = pd.read_csv('CheckM_data.tsv', sep='\t')
attributes_df = pd.read_csv('metadata_attributes.tsv', sep='\t')

# Merge the dataframes based on the "accession" column
merged_df = pd.merge(output_df, attributes_df, left_on='Strain', right_on='accession', how='inner')

# Drop the redundant "accession" column from the attributes dataframe
merged_df.drop('accession', axis=1, inplace=True)

# Rearrange the columns as you specified
column_order = ['Strain', 'organismName', 'strain', 'host', 'isolation_source', 'collection_date', 'geo_loc_name', 'sample_type',
                '#_genomes', '#_markers', '#_marker_sets', 'GC', 'GC_std', 'Genome_size', '#_ambiguous_bases', '#_scaffolds',
                '#_contigs', 'Longest_scaffold', 'Longest_contig', 'N50_(scaffolds)', 'N50_(contigs)', 'Mean_scaffold_length',
                'Mean_contig_length', 'Completeness', 'Contamination']

merged_df = merged_df[column_order]

# Write the merged dataframe to a new file
merged_df.to_csv('merged_output.tsv', sep='\t', index=False)
```

### Extract rows from metadata that correspond to dereplicated human-associated genomes that passed QC (see Section 02)

```python
# Read the accession numbers from accession_numbers.txt (list of accession names in dereplicated genomes from directory)
with open('accession_numbers.txt', 'r') as accession_file:
    accession_numbers = set(line.strip() for line in accession_file)

# Create a new file for the matching lines
with open('matching_lines.tsv', 'w') as output_file:
    # Write the header from merged_output.tsv
    with open('merged_output.tsv', 'r') as input_file:
        header = input_file.readline()
        output_file.write(header)

        # Iterate through the lines and write matching lines to the output file
        for line in input_file:
            accession = line.split('\t')[0]
            if accession in accession_numbers:
                output_file.write(line)

print("Matching lines have been extracted and saved to matching_lines.tsv")
```

```python
# Add "missing" to "Strain" entry if genomes is derived from MAG

# Define the input and output file names
input_file = "metadata_bif_genomes_combined.tsv"
output_file = "metadata_bif_genomes_combined_w_strain.tsv"

# Open the input and output files
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    # Read the header line and write it to the output file
    header = infile.readline()
    outfile.write(header)

    # Process the data lines
    for line in infile:
        parts = line.strip().split('\t')  # Split the line into columns
        if len(parts) >= 26:
            strain = parts[4]  # Assuming the strain column is the 5th column (0-based index)

            # Check if the strain is empty and add "missing" if necessary
            if not strain.strip():
                parts[4] = "missing"

            # Write the updated line to the output file
            outfile.write('\t'.join(parts) + '\n')

print("Processing complete. Updated data saved to", output_file)
```

### Extract relevant species and subspecies for metadata table

```python
# Open the input file for reading and the output file for writing
with open("metadata_public_GCA_genomes.tsv", "r") as input_file, open("output_file.tsv", "w") as output_file:
    # Read the header line and write it to the output file
    header = input_file.readline()
    output_file.write(header)
    
    # Process each line in the input file
    for line in input_file:
        # Split the line into columns
        columns = line.strip().split("\t")
        
        # Extract the first two words from the "organismName" column
        organism_name = columns[1]
        genus_species = " ".join(organism_name.split()[:2])
        
        # Update the "genus_species" column
        columns[2] = genus_species
        
        # Write the updated columns to the output file
        output_file.write("\t".join(columns) + "\n")

print("Processing complete.")
```

```python
# Open the input file for reading and the output file for writing
with open("metadata_public_GCA_genomes.tsv", "r") as input_file, open("output_subsp.tsv", "w") as output_file:
    # Read the header line and write it to the output file
    header = input_file.readline()
    output_file.write("subspecies\n")
    
    # Process each line in the input file
    for line in input_file:
        # Split the line into columns
        columns = line.strip().split("\t")
        
        # Extract the organismName and look for "subsp."
        organism_name = columns[1]
        subsp_index = organism_name.find("subsp.")
        
        # If "subsp." is found, extract the word after it
        if subsp_index != -1:
            subsp_name = organism_name[subsp_index + 6:].strip().split()[0]
        else:
            subsp_name = "NA"
        
        # Write the subsp_name to the output file
        output_file.write(subsp_name + "\n")

print("Processing complete.")
```

```python
# Define a dictionary to map organism names to hex codes
color_mapping = {
    "Bifidobacterium adolescentis": "#956cb4",
    "Bifidobacterium angulatum": "#8c613c",
    "Bifidobacterium animalis": "#bade28",
    "Bifidobacterium bifidum": "#440154",
    "Bifidobacterium breve": "#94b0c2",
    "Bifidobacterium catenulatum": "#6c92cc",
    "Bifidobacterium dentium": "#bd8e46",
    "Bifidobacterium gallicum": "#277f8e",
    "Bifidobacterium longum": "#a0da39",
    "Gardnerella vaginalis": "#ffff00",
    "Bifidobacterium pseudocatenulatum": "#55a868",
    "Bifidobacterium pseudolongum": "#c44e52",
    "Bifidobacterium pullorum": "#990000",
    "Bifidobacterium ruminantium": "#ffa500",
    "Bifidobacterium scardovii": "#d6859d",
    "Bifidobacterium thermophilum": "#ff0000",
}

# Input and output file paths
input_file = "add_colours_to_species.txt"
output_file = "species_with_colors.txt"

# Open input and output files
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    # Read and write the header line
    header = infile.readline().strip()
    outfile.write(f"{header}\tcolor\n")
    
    # Process each line in the input file
    for line in infile:
        accession, organism = line.strip().split("\t")
        # Get the hex code for the organism or use a default color if not found
        hex_code = color_mapping.get(organism, "#000000")
        outfile.write(f"{accession}\t{organism}\t{hex_code}\n")

print("Color mapping completed. Output written to", output_file)

```

### Add cleaned entries for isolation source, age group and country of origin

```python
# Isolation source clean

import pandas as pd
import numpy as np

# Read the input file
input_file = "isolation_source_clean.tsv"
df = pd.read_csv(input_file, sep='\t')

# Define a function to process isolation sources
def process_isolation_source(isolation_source):
    if pd.isna(isolation_source) or isolation_source == "missing":
        return "missing"
    
    keywords = ["feces", "stool", "fecal", "gut", "faece", "faeces", "intestinal"]
    if any(keyword in isolation_source for keyword in keywords):
        return "faeces"
    else:
        return isolation_source

# Apply the processing function to the "isolation_source" column and create the "clean" column
df["clean"] = df["isolation_source"].apply(process_isolation_source)

# Write the updated data to the output file
output_file = "isolation_source_clean_updated.tsv"
df.to_csv(output_file, sep='\t', index=False)

print("Processing complete. Updated data saved to", output_file)

```

```python
# Age group

import json

def parse_attributes(json_data):
    try:
        attributes = json_data.get('assemblyInfo', {}).get('biosample', {}).get('attributes', [])
        parsed_attributes = {}

        for attribute in attributes:
            name = attribute.get('name', '')
            value = attribute.get('value', '')
            if name == 'host_age':
                parsed_attributes[name] = value

        return parsed_attributes

    except json.JSONDecodeError:
        print("Invalid JSON format.")
        return None

if __name__ == "__main__":
    file_path = "bif_human_list_for_metadata.txt"
    output_file = "output_host_age.tsv"

    try:
        with open(file_path, 'r') as infile, open(output_file, 'w') as outfile:
            outfile.write("accession\thost_age\n")

            for line in infile:
                try:
                    json_data = json.loads(line)
                    accession = json_data.get('accession', '')
                    attributes_data = parse_attributes(json_data)

                    if attributes_data:
                        host_age = attributes_data.get('host_age', '')
                        outfile.write(f"{accession}\t{host_age}\n")
                except json.JSONDecodeError:
                    print("Invalid JSON format in line:", line)

        print("Output has been written to", output_file)

    except FileNotFoundError:
        print(f"File not found: {file_path}")

#########

## Next script ##

#########

# Combine cleaned age group to isolation source

import pandas as pd
import numpy as np

# Read the input file
input_file = "isolation_source_clean.tsv"
df = pd.read_csv(input_file, sep='\t')

# Define a function to process isolation sources
def process_isolation_source(isolation_source):
    if pd.isna(isolation_source) or isolation_source == "missing":
        return "missing"
    
    keywords_infant = ["infant", "infants", "Infant", "Infants", "baby"]
    keywords_adult = ["adult", "Adult", "adults", "Adults", "maternal"]
    
    if any(keyword in isolation_source for keyword in keywords_infant):
        return "infant"
    elif any(keyword in isolation_source for keyword in keywords_adult):
        return "adult"
    else:
        return isolation_source

# Apply the processing function to the "isolation_source" column and create the "age" column
df["age"] = df["isolation_source"].apply(process_isolation_source)

# Write the updated data to the output file
output_file = "isolation_source_age_updated.tsv"
df.to_csv(output_file, sep='\t', index=False)

print("Processing complete. Updated data saved to", output_file)
```

```python
# Country of origin for isolates
# Unique countries were assigned a latitude and longitude corresponding to their respective capital cities
# Where country was missing, 0,0 was used for the latitude and longitude (see next script), respectively

# Read the original data from the given file
with open("updated_geo_data.txt", "r") as file:
    lines = file.readlines()

# Create a dictionary to store the country-coordinate mappings
country_coordinates = {
    "Sweden": (59.3293, 18.0686),
    "Chile": (-33.4489, -70.6693),
    "China": (39.9042, 116.4074),
    "South Africa": (-25.7461, 28.1881),
    "United Kingdom": (51.5074, -0.1278),
    "Thailand": (13.7563, 100.5018),
    "South Korea": (37.5665, 126.9780),
    "Australia": (-35.2820, 149.1287),
    "Ireland": (53.3498, -6.2603),
    "Mozambique": (-25.8919, 32.6051),
    "India": (28.6139, 77.2090),
    "Canada": (45.4215, -75.6972),
    "Spain": (40.4168, -3.7038),
    "Russia": (55.7558, 37.6173),
    "France": (48.8566, 2.3522),
    "Czech Republic": (50.0755, 14.4378),
    "Netherlands": (52.3676, 4.9041),
    "Belgium": (50.8503, 4.3517),
    "Germany": (52.5200, 13.4050),
    "Kenya": (-1.2921, 36.8219),
    "Mexico": (19.4326, -99.1332),
    "Japan": (35.6895, 139.6917),
    "Korea": (37.5665, 126.9780),  # Duplicate with South Korea, corrected coordinates
    "Italy": (41.9028, 12.4964),
    "Brazil": (-15.7942, -47.8825),
    "Singapore": (1.3521, 103.8198),
    "USA": (38.8951, -77.0364)
}

# Process each line in the data
updated_lines = []
for line in lines:
    fields = line.strip().split("\t")
    country = fields[0]
    if country in country_coordinates and fields[1] != "NA" and fields[2] != "NA":
        new_latitude, new_longitude = country_coordinates[country]
        fields[1] = str(new_latitude)
        fields[2] = str(new_longitude)
    updated_lines.append("\t".join(fields))

# Write the updated data back to the file
with open("updated_geo_data.txt", "w") as file:
    file.write("\n".join(updated_lines))

print("Coordinates updated successfully!")

#########

## Next script ##

#########

import pandas as pd

# Load the TSV file
file_path = './bacteria_metadata.tsv'
df = pd.read_csv(file_path, sep='\t')

# Replace missing (empty) latitude and longitude values with 0
df['latitude'].fillna(0, inplace=True)
df['longitude'].fillna(0, inplace=True)

# Save the updated DataFrame back to a new TSV file
output_path = './updated_06_bacteria_metadata.tsv'
df.to_csv(output_path, sep='\t', index=False)

print(f"Updated file saved to {output_path}")
```

```python
# Assign metadata specific colour for data visualisation

import pandas as pd

# Load bifido_itol_list.txt
bifido_itol = pd.read_csv("bifido_itol_list.txt", sep="\t")

# Load 06_bacteria_metadata.tsv
bacteria_metadata = pd.read_csv("bacteria_metadata.tsv", sep="\t")

# Define palettes
source_palette = {
    "milk": "#D55E00",
    "eye": "#0072B2",
    "dental caries": "#E69F00",
    "missing": "lightgrey",
    "faeces": "#009E73",
    "vagina": "#CC79A7",
    "supplements": "#F0E442",
    "blood": "#56B4E9",
    "urine": "#0072B2"
}

age_palette = {
    "missing": "lightgrey",
    "NA": "lightgrey",
    "child": "#56B4E9",
    "elderly": "#E69F00",
    "infant": "#CC79A7",
    "adult": "#009E73"
}

country_palette = {
    "Italy": "#1f78b4",
    "Netherlands": "#33a02c",
    "South Africa": "#e31a1c",
    "Japan": "#ff7f00",
    "Czech Republic": "#6a3d9a",
    "Mexico": "#a6cee3",
    "USA": "#b2df8a",
    "France": "#fdbf6f",
    "Cambodia": "#cab2d6",
    "Belgium": "#fb9a99",
    "Canada": "#ffff99",
    "Singapore": "#8c510a",
    "Ireland": "#d6859d",
    "Kenya": "#addd8e",
    "South Korea": "#a6a6a6",
    "Chile": "#525252",
    "Brazil": "#ff7f7f",
    "Mozambique": "#a65628",
    "Germany": "#984ea3",
    "Australia": "#999999",
    "Russia": "#2ca25f",
    "Thailand": "#0570b0",
    "Vietnam": "#fdae61",
    "Sweden": "#fee08b",
    "China": "#ffff00",
    "United Kingdom": "#d73027",
    "India": "#4575b4",
    "Spain": "#91bfdb",
    "missing": "lightgrey"
}

# Merge the two dataframes on the Accession column
merged_df = pd.merge(bifido_itol, bacteria_metadata, left_on="Accession_itol", right_on="Accession", how="left")

# Extract relevant columns and add hex color codes
merged_df['isolation_source_clean'] = merged_df['isolation_source_clean'].fillna("missing")
merged_df['isolation_source_color'] = merged_df['isolation_source_clean'].map(source_palette).fillna("lightgrey")

merged_df['age'] = merged_df['age'].fillna("missing")
merged_df['age_color'] = merged_df['age'].map(age_palette).fillna("lightgrey")

merged_df['country'] = merged_df['country'].fillna("missing")
merged_df['country_color'] = merged_df['country'].map(country_palette).fillna("lightgrey")

# Select the relevant columns for the output
output_df = merged_df[['Accession_itol', 'isolation_source_clean', 'isolation_source_color', 'age', 'age_color', 'country', 'country_color']]

# Save the result to a new file
output_df.to_csv("bifido_itol_list_with_metadata.txt", sep="\t", index=False)

print("Script execution complete. Output file saved as 'bifido_itol_list_with_metadata.txt'.")
```

### Add cleaned entries for isolation source, age group and country of origin

### Add Seqfu assembly statistics to metadata (see Section 02)

# 02 | Quality checks and dereplication of bacterial genomes

## # Quality check and dereplicate all genomes at strain level (99.9% ANI)

### 01 | QC - CheckM (v1.2.0) (only genomes ≥90% completeness and **≤5% contamination would be used for final genome database)**

```bash
# CheckM v1.2.0 was run using NBI HPC services
checkm lineage_wf -t <specify no. of threads> -x fna /path/to/genome/directory/ /path/to/desired/output/folder/directory/should/not/exist
```

### 01 | QC - Extracting completeness and contamination scores from CheckM ouput

```python
import csv

input_file = 'bin_stats_ext.tsv'
output_file = 'output.tsv'

with open(input_file, 'r') as f_input, open(output_file, 'w', newline='') as f_output:
    reader = csv.reader(f_input, delimiter='\t')
    writer = csv.writer(f_output, delimiter='\t')
    
    for row in reader:
        accession_number = row[0]
        genome_info = eval(row[1])  # Convert the string to a dictionary
        completeness = genome_info.get('Completeness', '')
        contamination = genome_info.get('Contamination', '')
        
        writer.writerow([accession_number, completeness, contamination])

# Output is tab-delimited and looks like this:
GCA_028203155.1_ASM2820315v1_genomic    92.88888888888889       0.0
GCA_028203125.1_ASM2820312v1_genomic    99.55555555555556       0.0
```

### 01 | QC - Extract the accession numbers / strain codes of strains that have ≥90% completeness and **≤5%**

```python
# Open the input file
with open('output.tsv', 'r') as input_file:
    # Read the lines of the file
    lines = input_file.readlines()

# Create an empty list to store the filtered accession numbers
filtered_accessions = []

# Process each line in the file
for line in lines:
    # Split the line into columns
    columns = line.strip().split('\t')

    # Extract the accession number, completeness, and contamination
    accession = columns[0]
    completeness = float(columns[1])
    contamination = float(columns[2])

    # Check if the completeness is over 90% and contamination is under 5%
    if completeness > 90 and contamination < 5:
        # Add the accession number to the filtered list
        filtered_accessions.append(accession)

# Open the output file
with open('filtered_accessions.txt', 'w') as output_file:
    # Write the filtered accession numbers to the file
    for accession in filtered_accessions:
        output_file.write(accession + '\n')

print("Filtered accessions have been written to 'filtered_accessions.txt'")
```

### 01 | QC - Script for moving files that passed filtering thresholds to relevant directory

```bash
for X in `cat filtered_accessions.txt`; do mv $X.fna path/to/desired/directory/${SX} ; done
```

### 02 | Dereplication of genomes - Command for using dRep (v3.4.3) (99.9% ANI cutoff)

```bash
dRep dereplicate /home/ubuntu/checkm_filtered_all/drep/ -g /home/ubuntu/checkm_filtered_all/*.fna -sa 0.999 -p 16

# -sa option is the cutoff for dereplication (we used 99.9% ANI)
```

## # Step 2: Remove fragmented assemblies using Seqfu

### 03 | Standardise genome and contig names

```python
# Genome and contig names were standardised to include only accession

import os

# Function to standardize genome file names
def standardize_genome_name(filename):
    if filename.startswith("GCA_"):
        parts = filename.split("_")
        if len(parts) > 2:
            return "_".join(parts[:2]) + ".fna"
    
    prefixes = [
        "adolescentis_", "bifidum_", "breve_", "dentium_", 
        "faecale_", "longum_", "longum_sub_infantis_", 
        "pseudocatenulatum_"
    ]
    
    for prefix in prefixes:
        if filename.startswith(prefix):
            return filename[len(prefix):]
    
    return filename

# Function to rename files in the specified directory
def rename_files_in_directory(directory):
    for filename in os.listdir(directory):
        new_name = standardize_genome_name(filename)
        if new_name != filename:
            old_path = os.path.join(directory, filename)
            new_path = os.path.join(directory, new_name)
            os.rename(old_path, new_path)
            print(f"Renamed: {filename} -> {new_name}")

# Use the current directory
directory_path = os.getcwd()
rename_files_in_directory(directory_path)
```

```python
import os

# Function to read the mapping file
def read_mapping(file_path):
    mapping = {}
    with open(file_path, 'r') as f:
        next(f)  # Skip the header line
        for line in f:
            original, new = line.strip().split('\t')
            mapping[original] = new
    return mapping

# Function to rename files based on the mapping
def rename_files(mapping, directory):
    for original, new in mapping.items():
        original_file = os.path.join(directory, f"{original}.fna")
        new_file = os.path.join(directory, f"{new}.fna")
        if os.path.exists(original_file):
            os.rename(original_file, new_file)
            print(f"Renamed {original_file} to {new_file}")
        else:
            print(f"File {original_file} not found")

# Path to the mapping file
mapping_file_path = 'PEARL_strain_name_change.txt'
# Path to the directory containing the .fna files
fna_directory = './'  # Change this to your directory

# Read the mapping and rename the files
mapping = read_mapping(mapping_file_path)
rename_files(mapping, fna_directory)
```

```python
import os
from Bio import SeqIO

def standardize_contig_names(file_path):
    filename = os.path.basename(file_path)
    name_base, _ = os.path.splitext(filename)
    name_base = name_base.replace('.', '_')
    
    # Read the sequences from the file
    sequences = list(SeqIO.parse(file_path, "fasta"))
    
    # Rename the contigs
    for i, record in enumerate(sequences, start=1):
        record.id = f"{name_base}_{i}"
        record.description = ""
    
    # Write the renamed sequences back to the file
    with open(file_path, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
    
    print(f"Processed {filename}: contigs renamed to {name_base}_1, {name_base}_2, etc.")

# Process each .fna file in the current directory
current_directory = os.getcwd()
for file in os.listdir(current_directory):
    if file.endswith(".fna"):
        file_path = os.path.join(current_directory, file)
        standardize_contig_names(file_path)

```

### 03 | Command to use Seqfu

```bash
#!/bin/bash

# Define the output file
output_file="combined_output.txt"

# Create or clear the output file
> $output_file

# Add the header to the output file
echo -e "File\t#Seq\tTotal bp\tAvg\tN50\tN75\tN90\tauN\tMin\tMax" > $output_file

# Loop through each genome file
for genome_file in *.fna; do
  # Run seqfu stats on the genome file and append to the output file
  seqfu stats "$genome_file" | tail -n +2 >> $output_file
done

```

### 03 | Remove genomes with N50 below 50,000 and more than 100 contigs

```python
def process_genomes(input_file, passed_file, failed_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    header = lines[0]  # Save the header line

    passed_genomes = [header]
    failed_genomes = [header]

    for line in lines[1:]:
        columns = line.strip().split('\t')
        filename = columns[0]
        num_contigs = int(columns[1])  # #Seq column
        n50_value = int(columns[4])  # N50 column
        
        # Check the conditions
        if n50_value > 50000 and num_contigs < 100:
            passed_genomes.append(line)
        else:
            failed_genomes.append(line)

    # Write the results to the respective files
    with open(passed_file, 'w') as outfile:
        outfile.writelines(passed_genomes)

    with open(failed_file, 'w') as outfile:
        outfile.writelines(failed_genomes)

# Define the input and output files
input_file = '04_combined_output.txt'
passed_file = 'passed_genomes.txt'
failed_file = 'failed_genomes.txt'

# Process the genomes
process_genomes(input_file, passed_file, failed_file)

```

```python
# Transfer genomes that passed quality filtering

import os
import shutil

# Define the input file and the target directory
input_file = '05_passed_genomes.txt'
target_directory = 'passed_seqfu_genomes'

# Ensure the target directory exists
if not os.path.exists(target_directory):
    print(f"Error: The directory '{target_directory}' does not exist.")
    exit(1)

# Read the input file and extract the filenames from the "File" column
with open(input_file, 'r') as file:
    lines = file.readlines()

# Skip the header line
genome_files = [line.split('\t')[0] for line in lines[1:]]

# Process each genome file
for genome_file in genome_files:
    if os.path.exists(genome_file):
        # Move the file to the target directory
        shutil.move(genome_file, os.path.join(target_directory, genome_file))
        print(f"Moved: {genome_file}")
    else:
        print(f"File not found: {genome_file}")
```

# 03 | Annotation of genomes for standards-compliant output files using Prokka (**v1.14.6)**

## # Annotate all *Bifidobacterium* genomes using Prokka (**v1.14.6)** to generate files such as .faa, .gff and .gbk

```bash
# Following dereplication of human associated Bifidobacterium genomes, all files were deposited into a single directory
# Prokka was run on a loop with the following command

for k in *.fna; do prokka $k --centre --compliant --outdir "$k".prokka.output --prefix PROKKA_$k; echo $k; done

# The .fna extension may be carried over into the filenames. If this is the case, use the following Perl regular expression with the rename command

rename 's/(\.fna)//e' *.faa
```

# 04 | **Phylogenetic profiling of Bifidobacterium genomes using PhyloPhlAn (v3.0.67) and GTDB-Tk (v2.1.0)**

## # Use **PhyloPhlAn (v3.0.67) and GTDB-Tk (v2.1.0) for phylogenetic profiling / QC (to check strains have been appropriately named).**

### 01 | **PhyloPhlAn (v3.0.67) -** Phylogenetic and taxonomic classification of *Bifidobacterium* genomes

```bash
# PhyloPhlAn (v3.0.67) was run as a singularity image in a HPC environment
# .faa outputs from Prokka were used for this analysis

# PhyloPhlAn was run using the following command

IMAGE=/path/to/singulatiy/image/phylophlan-0.1.simg
COMMAND=phylophlan

singularity exec $IMAGE $COMMAND -i /path/to/faa/files/ -d phylophlan \
--databases_folder /path/to/databases/folder/phylophlan_databases/ --diversity low \
-f /path/to/configuration/file/supermatrix_aa.cfg --maas /path/to/mapping/file/phylophlan.tsv \
--fast --output_folder /path/to/output/folder/james_phylophlan_19July2023 --verbose --nproc 64
```

### Please note - the phylogenetic tree was later visualised and edited using iTol (v5)

### 02 | **GTDB-Tk (v2.1.0) - Assigning objective taxonomic classifications to bacterial genomes**

```bash
# GTDB-Tk (v2.1.0) was run as a singularity image in a HPC environment

gtdbtk classify_wf --cpus 20 -x fna --pplacer_cpus 1 --genome_dir /path/to/genome/files --out_dir /path/to/output/directory
```

### Changes made to genome names based on taxonomic classification

```
'Accession'	'Original Species Name'	'Rectified Species Name'
GCA_019129655.1	Bifidobacterium sp.	Bifidobacterium adolescentis
GCA_020686885.1	Bifidobacterium sp.	Bifidobacterium adolescentis
GCA_008669245.1	Bifidobacterium sp.	Bifidobacterium adolescentis
GCA_008662405.1	Bifidobacterium sp.	Bifidobacterium adolescentis
GCA_026016425.1	Bifidobacterium sp.	Bifidobacterium bifidum
GCA_023656805.1	Bifidobacterium sp.	Bifidobacterium breve
GCA_022723235.1	Bifidobacterium sp.	Bifidobacterium catenulatum
GCA_002742425.1	Bifidobacterium sp.	Bifidobacterium catenulatum
GCA_002742445.1	Bifidobacterium sp.	Bifidobacterium catenulatum
GCA_023656765.1	Bifidobacterium sp.	Bifidobacterium catenulatum
GCA_018367525.1	Bifidobacterium catenulatum	Bifidobacterium catenulatum subsp. kashiwanohense
GCA_026016385.1	Bifidobacterium sp.	Bifidobacterium longum
GCA_018368035.1	Bifidobacterium sp.	Bifidobacterium pullorum
GCA_030224545.1	Bifidobacterium sp.	Gardnerella vaginalis
LHP009	Bifidobacterium faecale	Bifidobacterium adolescentis
SRR14193525	Bifidobacterium pseudocatenulatum	Bifidobacterium catenulatum
SRR14193528	Bifidobacterium pseudocatenulatum	Bifidobacterium catenulatum
```

# 05 | Identification of anti-viral defence systems
using DefenseFinder and PADLOC

## # Using DefenseFinder and PADLOC to identify anti-viral defence systems in bifidobacteria. Scripts for data processing are also included.

## 01 | Running Defensefinder and PADLOC

### DefenseFinder: transfer required .faa files to relevant directory

```bash
# Following annotation using Prokka, all .faa files were copied to a separate folder
cp ./*/*.faa /path/to/destination/folder

# Create directory using .faa filenames and move the files to its respective folder
for file in *.faa; do mkdir -- "${file%.faa}"; mv -- "$file" "${file%.faa}"; done
```

### Submit DefenseFinder job on a loop (Bash script)

```bash
#!/bin/bash

# Loop through each directory
for dir in */; do
  # Remove trailing slash from directory name
  dir=${dir%/}
  
  # Construct the .faa file path
  faa_file="${dir}/${dir}.faa"
  
  # Run the defense-finder command if the .faa file exists
  if [[ -f "$faa_file" ]]; then
    echo "Running defense-finder for $faa_file"
    defense-finder run "$faa_file" --out-dir "defensefinder_$dir"
  else
    echo "File $faa_file not found. Skipping."
  fi
done

```

### PADLOC: transfer required .faa and .gff files to relevant directory

```bash
#!/bin/bash

# Source directory where your folders are located
SOURCE_DIR="/path/to/source/directory"

# Destination directory where you want to copy the files
DEST_DIR="/path/to/destination/directory"

# Loop through each folder in the source directory
for folder in "$SOURCE_DIR"/*; do
    if [ -d "$folder" ]; then
        # Extract the folder name (without path)
        folder_name=$(basename "$folder")

        # Construct the corresponding destination folder name
        dest_folder_name="${folder_name}_prokka_output"

        # Check if the destination folder exists; if not, create it
        mkdir -p "$DEST_DIR/$dest_folder_name"

        # Find and copy .gff and .faa files from the source folder to the destination folder
        for file in "$folder"/*.gff "$folder"/*.faa; do
            if [ -f "$file" ]; then
                cp "$file" "$DEST_DIR/$dest_folder_name/"
            else
                echo "Warning: No .gff or .faa files found in $folder"
            fi
        done
    fi
done

echo "Files copied successfully!"
```

### Submit PADLOC job on a loop (Bash script)

```bash
#!/bin/bash

# Loop through each subdirectory
for dir in */ ; do
    # Remove trailing slash from directory name
    dir=${dir%/}

    # Define the necessary file paths
    faa_file="${dir}/${dir}.faa"
    gff_file="${dir}/${dir}.gff"
    crispr_file="${dir}/${dir}_crispr.gff"

    # Check if the required files exist
    if [[ -f "$faa_file" && -f "$gff_file" && -f "$crispr_file" ]]; then
        # Run the padloc command with output directed to the respective directory
        padloc --faa "$faa_file" --gff "$gff_file" --crispr "$crispr_file" --cpu 12 --outdir "$dir"
        echo "Processed padloc for $dir, outputs saved in $dir"
    else
        echo "Required files for $dir not found. Skipping."
    fi
done
```

### Process PADLOC output to match the output from DefenseFinder

```python
import os
import pandas as pd
import glob

# Define the directory containing the padloc folders
directory = './'  # Replace this with the path to the main directory containing the folders

# Get a list of all folders starting with 'padloc_'
folders = [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f)) and f.startswith('padloc_')]

# Iterate through each padloc folder
for folder in folders:
    folder_path = os.path.join(directory, folder)

    # Find the file that ends with '_padloc.tsv'
    tsv_files = glob.glob(os.path.join(folder_path, '*_padloc.tsv'))

    # Check if there is exactly one matching file
    if len(tsv_files) == 1:
        tsv_file = tsv_files[0]

        # Extract the base name of the file (without extension)
        base_name = os.path.splitext(os.path.basename(tsv_file))[0]

        # Read the TSV file
        df = pd.read_csv(tsv_file, sep='\t')

        # Group by system.number
        grouped = df.groupby('system.number')

        # Create a new dataframe for the results
        result = pd.DataFrame(columns=['sys_id', 'type', 'subtype', 'sys_beg', 'sys_end', 'protein_in_syst', 'genes_count', 'name_of_profiles_in_sys', 'seqid'])

        # Process each system
        for sys_num, group in grouped:
            if 'CRISPR_array' not in group['system'].values:
                # Exclude CRISPR arrays
                new_row = pd.DataFrame({
                    'sys_id': [f'{base_name}_{sys_num}'],  # Append the base_name to the system.number
                    'type': [''],  # Initialize the "type" column with empty strings
                    'subtype': [group['system'].iloc[0]],
                    'sys_beg': [group['target.name'].iloc[0]],
                    'sys_end': [group['target.name'].iloc[-1]],
                    'protein_in_syst': [','.join(group['target.name'])],
                    'genes_count': [len(group)],
                    'name_of_profiles_in_sys': [','.join(group['hmm.name'].dropna().unique())],
                    'seqid': [group['seqid'].iloc[0]]  # Add seqid as the last column
                })
                result = pd.concat([result, new_row], ignore_index=True)

        # Sort the result dataframe by sys_id
        result = result.sort_values('sys_id')

        # Reset the index
        result = result.reset_index(drop=True)

        # Define the output file name
        output_file = os.path.join(folder_path, base_name + '_parsed.tsv')

        # Save the result to a new TSV file
        result.to_csv(output_file, sep='\t', index=False)

        print(f"Parsing complete for {tsv_file}. Results saved to '{output_file}'.")

    else:
        print(f"Warning: No or multiple '_padloc.tsv' files found in {folder_path}.")
```

## 02 | Data Processing

### Split RM HNH into a single defence system

```python
import pandas as pd

# Load the data
df = pd.read_csv('01_updated_with_tool_original.tsv', sep='\t')

# Function to check if all proteins in the DefenseFinder entry are present in the padloc entries
def matching_proteins(padloc_proteins, defensefinder_proteins):
    padloc_proteins_set = set(padloc_proteins.split(','))
    defensefinder_proteins_set = set(defensefinder_proteins.split(','))
    return defensefinder_proteins_set.issubset(padloc_proteins_set)

# Identify genomes with padloc RM_type_HNH or RM_type_II
padloc_rm = df[(df['tool_original'] == 'padloc') & (df['subtype'].isin(['RM_type_HNH', 'RM_type_II']))]

# Filter out DefenseFinder RM_Type_II rows with matching genomes and matching protein_in_syst
defense_rm = df[(df['tool_original'] == 'DefenseFinder') & (df['subtype'] == 'RM_Type_II')]

# Loop through DefenseFinder RM_Type_II entries and check for matching proteins in the corresponding padloc entries
rows_to_remove = []
for index, row in defense_rm.iterrows():
    # Get the corresponding padloc entries for the same genome
    padloc_entries = padloc_rm[padloc_rm['Genome'] == row['Genome']]
    
    # Combine the proteins from padloc RM_type_HNH and RM_type_II
    combined_padloc_proteins = ','.join(padloc_entries['protein_in_syst'].unique())
    
    # Check if the DefenseFinder proteins are a subset of the combined padloc proteins
    if matching_proteins(combined_padloc_proteins, row['protein_in_syst']):
        rows_to_remove.append(index)

# Drop the DefenseFinder rows where proteins matched with padloc
df_filtered = df.drop(rows_to_remove)

# Save the filtered data
df_filtered.to_csv('02_deduplicated_HNH_filtered.tsv', sep='\t', index=False)

```

### Deduplicate defence systems

```python
import pandas as pd

# Load the data into a DataFrame
df = pd.read_csv('02_deduplicated_HNH_split.tsv', sep='\t')

# Helper function to get tool_final value based on specific rules
def get_tool_final(subtype, tools, original_tool):
    if subtype == "mza":
        return "padloc"
    elif subtype == "Cas_Type_II-C" or subtype == "Abi2":
        return original_tool
    elif len(tools) == 2:
        return 'both'
    return tools.pop()

# To store non-deduplicated entries
non_deduplicated = []

# To store the deduplicated and concatenated systems
deduplicated = []

# To store successfully deduplicated entries
successfully_deduplicated = []

# Define specific subtype rules
subtype_rules = {
    ("Lamassu-Mrr", "Lamassu-Family"): "Lamassu-Mrr",
    ("Paris_fused", "PARIS_II_merge"): "PARIS_II_merge",
    ("RM_Type_IIG", "BREX_I"): "BREX_I",
    ("Cas_Type_II-C", "Cas_Cluster"): "Cas_Type_II-C",
    ("mza", "RM_Type_II"): "mza",
    ("Septu", "Septu_Type_I"): "Septu_Type_I",
    ("Thoeris_I", "Thoeris_II"): "Thoeris_II",
    ("AbiD", "Abi2"): "Abi2",
    ("Hachiman_Type_I", "Hachiman"): "Hachiman_Type_I",
    ("Cas_Type_II-C", "Cas_Type_II-A"): "Cas_Type_II-C"
}

# Function to handle subtype merging based on rules
def resolve_subtype(subtype1, subtype2):
    if (subtype1, subtype2) in subtype_rules:
        return subtype_rules[(subtype1, subtype2)]
    elif (subtype2, subtype1) in subtype_rules:
        return subtype_rules[(subtype2, subtype1)]
    else:
        return None

# Group by Genome
grouped = df.groupby('Genome')

# Loop through each group to identify overlaps and deduplicate
for genome, group in grouped:
    for index, row in group.iterrows():
        proteins = set(row['protein_in_syst'].split(','))

        # Check for overlaps in previously stored proteins
        overlap_found = False
        preferred_entry = None
        for dedup_index, dedup_row in enumerate(deduplicated):
            if dedup_row['Genome'] == genome:
                existing_proteins = set(dedup_row['protein_in_syst'].split(','))

                # If overlap found, check if subtype matches or needs resolution
                if proteins & existing_proteins:
                    # Check if subtypes match
                    if dedup_row['subtype'] == row['subtype']:
                        # Subtypes match, so merge the systems
                        dedup_row['protein_in_syst'] = ','.join(sorted(existing_proteins | proteins))
                        dedup_row['tool_final'] = get_tool_final(dedup_row['subtype'], set([dedup_row.get('tool_final', ''), row['tool_original']]), dedup_row.get('tool_final', row['tool_original']))
                        overlap_found = True

                        # Log successful deduplication
                        successfully_deduplicated.append({
                            'Genome': genome,
                            'original_sys_id': row['sys_id'],
                            'merged_into_sys_id': dedup_row['sys_id'],
                            'subtype': dedup_row['subtype'],
                            'protein_in_syst': ','.join(sorted(existing_proteins | proteins)),
                            'tool_final': dedup_row['tool_final']
                        })
                        break
                    else:
                        # Check if there is a rule for resolving the subtype conflict
                        resolved_subtype = resolve_subtype(dedup_row['subtype'], row['subtype'])
                        if resolved_subtype:
                            # Determine which row to keep based on resolved subtype
                            if resolved_subtype == dedup_row['subtype']:
                                # Keep the existing `dedup_row` as it has the preferred subtype
                                dedup_row['protein_in_syst'] = ','.join(sorted(existing_proteins | proteins))
                                dedup_row['tool_final'] = get_tool_final(resolved_subtype, set([dedup_row.get('tool_final', ''), row['tool_original']]), row['tool_original'])
                                preferred_entry = dedup_row
                            else:
                                # Replace `dedup_row` with `row` since `row` has the preferred subtype
                                deduplicated[dedup_index] = row.to_dict()
                                deduplicated[dedup_index]['subtype'] = resolved_subtype
                                deduplicated[dedup_index]['protein_in_syst'] = ','.join(sorted(existing_proteins | proteins))
                                deduplicated[dedup_index]['tool_final'] = get_tool_final(resolved_subtype, set([dedup_row.get('tool_final', ''), row['tool_original']]), row['tool_original'])
                                preferred_entry = deduplicated[dedup_index]
                            overlap_found = True

                            # Log successful deduplication
                            successfully_deduplicated.append({
                                'Genome': genome,
                                'original_sys_id': row['sys_id'] if resolved_subtype == row['subtype'] else dedup_row['sys_id'],
                                'merged_into_sys_id': dedup_row['sys_id'] if resolved_subtype == dedup_row['subtype'] else row['sys_id'],
                                'subtype': resolved_subtype,
                                'protein_in_syst': preferred_entry['protein_in_syst'],
                                'tool_final': preferred_entry['tool_final']
                            })
                            break
        
        if not overlap_found:
            # No overlap or conflicting subtypes without resolution, add as new entry
            row_dict = row.to_dict()
            row_dict['tool_final'] = row['tool_original']  # Ensure 'tool_final' is initialized
            deduplicated.append(row_dict)
        
        # Check for overlaps with different subtypes that are unresolved
        for non_dedup_row in deduplicated:
            if non_dedup_row['Genome'] == genome and non_dedup_row['subtype'] != row['subtype']:
                existing_proteins = set(non_dedup_row['protein_in_syst'].split(','))
                if proteins & existing_proteins:
                    # If there is no rule for resolving, log the conflict
                    if not resolve_subtype(non_dedup_row['subtype'], row['subtype']):
                        non_deduplicated.append({
                            'Genome': genome,
                            'protein_overlap': ','.join(sorted(proteins & existing_proteins)),
                            'subtype_1': row['subtype'],
                            'subtype_2': non_dedup_row['subtype']
                        })

# Convert lists to DataFrames for output
deduplicated_df = pd.DataFrame(deduplicated)
non_deduplicated_df = pd.DataFrame(non_deduplicated)
successfully_deduplicated_df = pd.DataFrame(successfully_deduplicated)

# Save the deduplicated results
deduplicated_df.to_csv('03_final_deduplicated_defense_systems.tsv', sep='\t', index=False)

# Save the non-deduplicated overlaps
non_deduplicated_df.to_csv('03_final_non_deduplicated_overlaps.tsv', sep='\t', index=False)

# Save the successfully deduplicated list
successfully_deduplicated_df.to_csv('03_final_successfully_deduplicated_systems.tsv', sep='\t', index=False)
```

### QC | Check for defence systems where deduplication failed

```python
import pandas as pd
from collections import defaultdict

# Load the data
df = pd.read_csv('03_final_deduplicated_defense_systems.tsv', sep='\t')

# Dictionary to store overlapping proteins for each genome
overlaps = defaultdict(list)

# Group by the Genome column
for genome, group in df.groupby('Genome'):
    # Set to track seen proteins within the genome
    seen_proteins = set()
    # Set to track overlapping proteins within the genome
    overlapping_proteins = set()
    
    # Iterate through the proteins in each row of the group
    for proteins in group['protein_in_syst']:
        # Handle NaN values
        if pd.notna(proteins):
            # Split proteins by comma
            protein_list = proteins.split(',')
            
            # Check for overlaps
            for protein in protein_list:
                if protein in seen_proteins:
                    overlapping_proteins.add(protein)
                else:
                    seen_proteins.add(protein)
    
    # If there are overlapping proteins, store them
    if overlapping_proteins:
        overlaps[genome] = list(overlapping_proteins)

# Write the results to a TSV file
with open('04_overlapping_proteins_results.tsv', 'w') as f:
    if overlaps:
        f.write("Genome\tOverlapping Proteins\n")
        for genome, proteins in overlaps.items():
            f.write(f"{genome}\t{', '.join(proteins)}\n")
    else:
        f.write("No overlapping proteins found across genomes.\n")
```

## 03 | Exploring distribution of defence systems

### Get counts for defence system types

```python
import pandas as pd

# Read the input TSV file
input_file = '03_final_deduplicated_defense_systems.tsv'
df = pd.read_csv(input_file, sep='\t')

# Directly use the 'Genome' column for grouping
# Group by 'Genome' and 'type', then count occurrences
counts = df.groupby(['Genome', 'type']).size().reset_index(name='count')

# Reorder columns
counts = counts[['Genome', 'type', 'count']]

# Write the output to a new TSV file
output_file = '05_DS_type_counts_per_genome.tsv'
counts.to_csv(output_file, sep='\t', index=False)

print(f'Counts have been written to {output_file}')
```

### Create defence system matrix for iTol

```python
import pandas as pd

# Define the input and output file names
input_file = "05_DS_type_counts_per_genome.tsv"
output_file = "06_binary_matrix_defence_systems.tsv"

# Read the DefenseFinder output file into a pandas DataFrame
df = pd.read_csv(input_file, sep="\t")

# Pivot the table to create a presence/absence matrix
binary_matrix = df.pivot_table(index='Genome', columns='Type', values='Count', aggfunc='sum', fill_value=0)

# Convert the counts to presence/absence (1/0) - if count > 0, set to 1
binary_matrix = binary_matrix.applymap(lambda x: 1 if x > 0 else 0)

# Write the binary matrix to a tab-separated text file
binary_matrix.to_csv(output_file, sep="\t")

print(f"Binary matrix has been written to {output_file}")
```

```python
# Reorder matrix to match the order of entries in iTol

# Open the itol_genome_names.tsv file and read the node IDs into a list
with open('06_itol_genomes_list.txt', 'r') as f:
    node_ids = [line.strip() for line in f]

# Create a dictionary to store the matches
matches = {}

# Open the defense_systems_matrix.tsv file and read the lines
with open('06_binary_matrix_defence_systems.tsv', 'r') as f:
    # Read the header line containing defense system names
    header = f.readline().strip().split('\t')[1:]
    # Find the indexes of the defense systems in itol_genome_names.tsv
    defense_system_indexes = [header.index(system) for system in header]

    for line in f:
        # Split the line into the genome name and presence of defense systems
        fields = line.strip().split('\t')
        genome_name = fields[0]
        defense_system_presence = fields[1:]

        # Check if the genome name is in the list of node IDs from the other file
        for i, node in enumerate(node_ids):
            if node in genome_name:
                # If it is, add the match to the dictionary
                matches[node] = defense_system_presence

# Write the matches to a new file
with open('07_reordered_matrix_output.tsv', 'w') as f:
    # Write the header line
    f.write("Genomes\t" + "\t".join(header) + "\n")
    # Write the genome names and defense system presence
    for node in node_ids:
        f.write(node + '\t' + '\t'.join(matches[node]) + '\n')
```

### Get defence system counts for each genome (for calculating average number of defence systems in bifidobacteria and per species)

```python
import csv

# Define the input and output file paths
input_file = '05_DS_type_counts_per_genome.tsv'
output_file = '08_genome_total_counts.tsv'

# Create an empty dictionary to store the total counts per genome
genome_counts = {}

# Open the input TSV file and read it
with open(input_file, 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    
    # Iterate through each row and sum the counts by genome
    for row in reader:
        genome = row['Genome']
        count = int(row['Count'])
        
        if genome in genome_counts:
            genome_counts[genome] += count
        else:
            genome_counts[genome] = count

# Write the results to a new TSV file
with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    
    # Write the header
    writer.writerow(['Genome', 'Count'])
    
    # Write the total counts for each genome
    for genome, total_count in genome_counts.items():
        writer.writerow([genome, total_count])

print(f"Output saved to {output_file}")

```

```python
# Add species to accession list

import csv

# Define file paths
counts_file = '08_genome_total_counts.tsv'
metadata_file = 'updated_06_bacteria_metadata.tsv'
output_file = '09_genome_total_counts_with_species.tsv'

# Create a dictionary to store organism names by accession
accession_to_species = {}

# Open and read the metadata file
with open(metadata_file, 'r') as meta_file:
    reader = csv.DictReader(meta_file, delimiter='\t')
    
    # Iterate over the metadata file and extract accession and organismName
    for row in reader:
        accession = row['Accession']
        organism_name = row['organismName']
        
        # Split the organismName and take the second word
        if len(organism_name.split()) > 1:
            second_word = organism_name.split()[1]
        else:
            second_word = organism_name  # In case the name is only one word
        
        # Store the second word by accession
        accession_to_species[accession] = second_word

# Open the counts file and the output file
with open(counts_file, 'r') as counts_file, open(output_file, 'w', newline='') as output_file:
    counts_reader = csv.DictReader(counts_file, delimiter='\t')
    fieldnames = ['Genome', 'Count']
    
    # Create a writer for the output file
    writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()

    # Iterate through each row in the counts file
    for row in counts_reader:
        genome = row['Genome']
        count = row['Count']
        
        # Check if the genome matches an accession in the metadata file
        if genome in accession_to_species:
            # Append the second word of the organism name with an underscore
            new_genome = f"{accession_to_species[genome]}_{genome}"
        else:
            # If no match is found, keep the genome name as is
            new_genome = genome
        
        # Write the updated genome and count to the output file
        writer.writerow({'Genome': new_genome, 'Count': count})

print(f"Output saved to {output_file}")

```

```r
# R script - Shapiro-Wilk Test

# Load necessary libraries
library(readr)

# Read the data from the TSV file
data <- read_tsv("11_genome_total_counts_with_species.tsv")

# Perform the Shapiro-Wilk test on the Count column
shapiro_test <- shapiro.test(data$Count)

# Print the test results to the console (optional)
print(shapiro_test)

# Save the Shapiro-Wilk test results to a text file
output_file <- "shapiro_wilk_test_results.txt"
sink(output_file)  # Redirect the output to the file
cat("Shapiro-Wilk Test for Normality\n")
cat("-------------------------------\n")
cat("W Statistic:", shapiro_test$statistic, "\n")
cat("p-value:", shapiro_test$p.value, "\n")
cat("Alternative hypothesis: The data is not normally distributed.\n")
sink()  # Stop redirecting output

# Confirmation message
cat("Shapiro-Wilk test results saved to", output_file, "\n")

########

# Script for violin plots and stats tests

########

# Load necessary libraries
library(ggplot2)
library(viridis)
library(dplyr)
library(FSA)
library(ggsignif)

# Read the TSV file
data <- read.delim("11_genome_total_counts_with_species.tsv", header = TRUE, sep = "\t")

# Extract the species names from the Genome column
data$Species <- gsub("_.*", "", data$Genome)

# Filter out species with less than 5 entries
species_counts <- table(data$Species)
filtered_species <- names(species_counts[species_counts >= 5])
data <- data[data$Species %in% filtered_species, ]

# Define the species order and labels (filtered)
species_order <- c(
  "adolescentis", "bifidum", "breve", "catenulatum",
  "dentium", "longum", "pseudocatenulatum",
  "pseudolongum"
)
species_labels <- c(
  "B. adolescentis", "B. bifidum", "B. breve", "B. catenulatum",
  "B. dentium", "B. longum", "B. pseudocatenulatum",
  "B. pseudolongum"
)
species_order <- species_order[species_order %in% filtered_species]

# Generate a colorblind-friendly color palette
color_palette <- viridis(length(species_order), option = "D", begin = 0.2, end = 0.8)

# Numbered x-axis
data$SpeciesNumber <- match(data$Species, species_order)

# Create the combined violin and box plot with numbered x-axis
combined_plot <- ggplot(data, aes(x = factor(SpeciesNumber), y = Count, fill = Species)) +
  geom_violin() +
  geom_boxplot(width = 0.1, position = position_dodge(0.75), alpha = 0.5) +
  geom_text(aes(x = SpeciesNumber, y = 0, label = paste("n =", table(data$Species)[Species])),
            vjust = 1.2, hjust = 0.5, size = 3) +  # Add n = X labels
  labs(x = "Species Number", y = "No. of Defense Systems", title = "Defense Systems Count by Species") +
  scale_x_discrete(labels = species_labels) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(15, 30, 15, 15, "pt"),  # Adjust margins
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add border around the plot
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))  # Attach y-axis label to the box
  ) +
  ylim(0, 30) +  # Increased Y-limit to provide more space for significance annotations
  geom_signif(
    comparisons = list(
      c(1, 5),  # adolescentis - dentium
      c(2, 5),  # bifidum - dentium
      c(3, 5),  # breve - dentium
      c(4, 5),  # catenulatum - dentium
      c(5, 6),  # dentium - longum
      c(1, 7),  # adolescentis - pseudocatenulatum
      c(2, 7),  # bifidum - pseudocatenulatum
      c(6, 7)   # longum - pseudocatenulatum
    ),
    map_signif_level = TRUE,
    annotations = c("****", "***", "**", "*", "***", "****", "*", "****"),  # Based on adjusted p-values
    y_position = c(25, 23, 26, 28, 27, 29, 24, 26)  # Adjusted positions to avoid overlapping with violin plots
  )

# Save the combined plot as a high-quality PDF
ggsave("defense_systems_plot_with_significance_updated.pdf", plot = combined_plot, width = 15, height = 10, units = "in", dpi = 300)

# Print a message indicating that the plot has been saved
cat("Plot saved as defense_systems_plot_with_significance_updated.pdf\n")

# Calculate and print median and mean for each species
summary_data <- data %>%
  group_by(Species) %>%
  summarize(median = median(Count), mean = mean(Count))

print("Median and Mean for each species:")
print(summary_data)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(Count ~ Species, data = data)
print("Kruskal-Wallis Test:")
print(kruskal_result)

# Perform Dunn's test if Kruskal-Wallis test is significant
if (kruskal_result$p.value < 0.05) {
  dunn_result <- dunnTest(Count ~ Species, data = data, method = "bh")  # Adjust for multiple comparisons using Benjamini-Hochberg method
  print("Dunn's Test Results:")
  print(dunn_result)
  
  # Save results to a text file
  sink("kruskal_dunn_test_results.txt")
  cat("Kruskal-Wallis Test:\n")
  print(kruskal_result)
  cat("\n\nDunn's Test Results:\n")
  print(dunn_result)
  sink()  # Close the file
  cat("Kruskal-Wallis and Dunn's test results saved to kruskal_dunn_test_results.txt\n")
}
```

### Get percentage occurrence of defence system types

```python
import pandas as pd

# Load the data
file_path = '06_binary_matrix_defence_systems.tsv'  # Ensure the path is correct
data = pd.read_csv(file_path, sep='\t')

# Drop the "Genome" column, as we only want to analyze the defense system columns
defense_systems = data.drop(columns=['Genome'])

# Calculate the percentage occurrence of each defense system
percentage_occurrence = (defense_systems.sum() / len(defense_systems)) * 100

# Sort defense systems by percentage in descending order
sorted_percentage = percentage_occurrence.sort_values(ascending=False)

# Save the results to a txt file
output_file = 'defense_systems_percentage.txt'
with open(output_file, 'w') as file:
    for system, percentage in sorted_percentage.items():
        file.write(f"{system}: {percentage:.2f}%\n")

print(f"Results saved to {output_file}")
```

### Chi-square test for significance of defence system types between individual genomes and species (R script)

```r
### Individual genomes ###

# Load necessary library
library(tidyr)

# Load the data
file_path <- "DS_type_counts_per_genome.tsv"  # Ensure this path is correct
data <- read.table(file_path, sep = "\t", header = TRUE)

# Create the contingency table
contingency_table <- pivot_wider(data, names_from = Type, values_from = Count, values_fill = 0)

# Remove the Genome column to get only the counts for chi-square test
counts_matrix <- as.matrix(contingency_table[, -1])  # Exclude the first column (Genome)

# Perform the Chi-square test with Monte Carlo simulation
chi_square_sim_result <- chisq.test(counts_matrix, simulate.p.value = TRUE, B = 10000)

# Display the results
cat("Chi-square Statistic:", chi_square_sim_result$statistic, "\n")
cat("p-value from Monte Carlo simulation:", chi_square_sim_result$p.value, "\n")
cat("Degrees of Freedom (approximate):", chi_square_sim_result$parameter, "\n")
cat("Expected Frequencies:\n")
print(chi_square_sim_result$expected)
```

```r
### Between species ###

# Load necessary libraries
library(tidyr)
library(dplyr)

# Load the data
file_path <- "binary_matrix_defence_systems_with_species.tsv"  # Ensure the path is correct
data <- read.table(file_path, sep = "\t", header = TRUE)

# Extract species from the "Genome" column
data <- data %>%
  mutate(Species = sub("_.*", "", Genome))  # Extract the species name before "_"

# Group data by Species and calculate the sum of each defense system within each species
species_defense_systems <- data %>%
  select(-Genome) %>%         # Exclude the Genome column
  group_by(Species) %>%
  summarise(across(everything(), sum))  # Sum defense systems within each species

# Convert to matrix format for the Chi-square test
species_defense_matrix <- as.matrix(species_defense_systems[,-1])  # Remove the "Species" column

# Perform the Chi-square test with Monte Carlo simulation
chi_square_sim_result <- chisq.test(species_defense_matrix, simulate.p.value = TRUE, B = 10000)

# Display results
cat("Chi-square Statistic:", chi_square_sim_result$statistic, "\n")
cat("p-value from Monte Carlo simulation:", chi_square_sim_result$p.value, "\n")
cat("Degrees of Freedom (approximate):", chi_square_sim_result$parameter, "\n")
cat("Expected Frequencies:\n")
print(chi_square_sim_result$expected)
```

## 04 | Co-occurrence analysis of defence systems

### Before running CoinFinder, the subtypes corresponding to each genomes were extracted

```python
import pandas as pd

# Read the data from the TSV file
data = pd.read_csv("deduplicated_defense_systems.tsv", sep='\t')

# Select only the 'subtype' and 'Genome' columns
selected_data = data[['subtype', 'Genome']]

# Write the selected data to a new TSV file
selected_data.to_csv("subtype_genome_output_for_coinfinder.tsv", index=False, sep='\t', header=True)

# Optional: print the first few rows to verify the result
print(selected_data.head())
```

### Commands to run Coinfinder (Associate and dissociate modules)

```bash
### Association ###

coinfinder -i subtype_genome_output_for_coinfinder.tsv -p bifidobacterium_tree_coinfinder.txt -o <output prefix> --associate

### Dissociation ###

coinfinder -i subtype_genome_output_for_coinfinder.tsv -p bifidobacterium_tree_coinfinder.txt -o <output prefix> --dissociate

### Association network built in Cytoscape ###
```

### Localisation analysis

```python
# Step 1 - Extract co-occurring pairs of defence systems:

import pandas as pd

# File paths
dedup_file_path = 'deduplicated_defense_systems.tsv'
coinfinder_file_path = 'coinfinder_associate_output_pairs.tsv'

# Read the TSV files
dedup_df = pd.read_csv(dedup_file_path, sep='\t')
coinfinder_df = pd.read_csv(coinfinder_file_path, sep='\t')

# List to store the filtered pairs with sys_beg and sys_end information
result_pairs = []

# Iterate through each unique genome in the dedup_df
for genome in dedup_df['Genome'].unique():
    # Filter rows for the current genome
    genome_data = dedup_df[dedup_df['Genome'] == genome]
    
    # Find matching source-target pairs where both are present in this genome
    for _, row in coinfinder_df.iterrows():
        source_subtype = row['Source']
        target_subtype = row['Target']
        
        # Check if both source and target subtypes exist in the genome data
        source_data = genome_data[genome_data['subtype'] == source_subtype]
        target_data = genome_data[genome_data['subtype'] == target_subtype]
        
        if not source_data.empty and not target_data.empty:
            # Extract sys_beg and sys_end for the source and target
            source_sys_beg = source_data.iloc[0]['sys_beg']
            source_sys_end = source_data.iloc[0]['sys_end']
            target_sys_beg = target_data.iloc[0]['sys_beg']
            target_sys_end = target_data.iloc[0]['sys_end']
            
            # Add the entry to the result list
            result_pairs.append({
                'Genome': genome,
                'Source': source_subtype,
                'Source_sys_beg': source_sys_beg,
                'Source_sys_end': source_sys_end,
                'Target': target_subtype,
                'Target_sys_beg': target_sys_beg,
                'Target_sys_end': target_sys_end
            })

# Convert the result to a DataFrame
result_df = pd.DataFrame(result_pairs)

result_df.to_csv('filtered_coinfinder_pairs_with_sys_info.tsv', sep='\t', index=False)
```

```python
# Step 2 - Standardise CDS codes from Prokka to calculate genomic distance:

import pandas as pd

# File path
input_file_path = '02_filtered_coinfinder_pairs_with_sys_info.tsv'
output_file_path = '03_cleaned_coinfinder_pairs.tsv'

# Read the TSV file
df = pd.read_csv(input_file_path, sep='\t')

# Function to clean the sys_beg and sys_end columns
def clean_sys_column(value):
    # Remove all characters before the last underscore
    return value.split('_')[-1]

# Apply the function to each of the sys_beg and sys_end columns
df['Source_sys_beg'] = df['Source_sys_beg'].apply(clean_sys_column)
df['Source_sys_end'] = df['Source_sys_end'].apply(clean_sys_column)
df['Target_sys_beg'] = df['Target_sys_beg'].apply(clean_sys_column)
df['Target_sys_end'] = df['Target_sys_end'].apply(clean_sys_column)

# Display the cleaned DataFrame
print(df)

# Optionally, save the cleaned DataFrame to a new file
df.to_csv(output_file_path, sep='\t', index=False)
```

```python
# Step 3 - Calculate genomic distances between co-occurring pairs:

import pandas as pd

# File path
input_file_path = '03_cleaned_coinfinder_pairs.tsv'
output_file_path = '04_distance_calculated_coinfinder_pairs.tsv'

# Read the TSV file
df = pd.read_csv(input_file_path, sep='\t')

# Convert relevant columns to integers
df['Source_sys_beg'] = df['Source_sys_beg'].astype(int)
df['Source_sys_end'] = df['Source_sys_end'].astype(int)
df['Target_sys_beg'] = df['Target_sys_beg'].astype(int)
df['Target_sys_end'] = df['Target_sys_end'].astype(int)

# Function to calculate the distance between systems
def calculate_distance(row):
    # Determine which system starts later
    if row['Source_sys_beg'] > row['Target_sys_beg']:
        later_start = row['Source_sys_beg']
        earlier_end = row['Target_sys_end']
    else:
        later_start = row['Target_sys_beg']
        earlier_end = row['Source_sys_end']
    
    # Calculate the distance
    distance = later_start - earlier_end
    return distance

# Apply the function to each row in the DataFrame to calculate the distance
df['Distance'] = df.apply(calculate_distance, axis=1)

# Display the DataFrame with the new Distance column
print(df)

# Optionally, save the modified DataFrame to a new file
df.to_csv(output_file_path, sep='\t', index=False)
```

```python
# Step 4 - Calculate percentage of co-occurring pairs that are within 20 genes:

import pandas as pd

# File path
input_file_path = '04_distance_calculated_coinfinder_pairs.tsv'

# Read the TSV file
df = pd.read_csv(input_file_path, sep='\t')

# Create a list to store the results
results = []

# Group the DataFrame by the unique Source-Target pairs
grouped = df.groupby(['Source', 'Target'])

# Iterate through each group
for (source, target), group in grouped:
    # Total occurrences for this Source-Target pair
    total_occurrences = len(group)
    
    # Count occurrences where Distance is 20 or less (including negative values)
    within_20 = group[group['Distance'] <= 20].shape[0]
    
    # Calculate the percentage
    percentage_within_20 = (within_20 / total_occurrences) * 100
    
    # Add the result to the list
    results.append({
        'Source': source,
        'Target': target,
        'Occurrences': f"{within_20}/{total_occurrences}",
        'Percentage': f"{percentage_within_20:.2f}"
    })

# Convert the results to a DataFrame
results_df = pd.DataFrame(results)

# Display the DataFrame
print(results_df)

# Save the results to a new file
output_file_path = '05_distance_summary.tsv'
results_df.to_csv(output_file_path, sep='\t', index=False)
```

# 06 | Identification of prophages in bifidobacteria

## # Use VIBRANT (v1.2.1), geNomad (v1.6.1) and PhageBoost (v0.1.7) to recover prophages sequences.

Please note - the initial prophage identification was performed before Seqfu filtering and standardisation of genome / contig names.

### Step 1 | Running each tool

```bash
# For each tool, create a list of all genome files in the directory
ls *.fna -1 > {filename}.txt

# Remove ".fna" from each line
sed -i 's/\.fna$//' {filename}.txt

##############

# Script for running VIBRANT
# Run the job on a loop
for X in `cat {filename}.txt` ; do srun VIBRANT_run.py -i ${X}.fna ; done

##############

# Script for running geNomad
for X in `cat {filename}.txt` ; do genomad end-to-end --cleanup ${X}.fna /output/directory/genomad_${X} /genomad_db/ ; done

##############

# Script for running PhageBoost 
for X in `cat all_genomes_list.txt` ; do PhageBoost -f ${X}.fna -o phageboost_${X} ; done
```

### Step 2 | Extracting prophage coordinates from each tool

1) VIBRANT

```python
# Script 1 - 01_extract_phage_info_tsv.py

import os

# Open the output file for writing
with open("output.tsv", "w") as outfile:
    # Write the header row
    outfile.write("Directory\tContig\n")

    # Loop over all directories in the current directory
    for dirname in os.listdir("."):
        # Check if the directory starts with "VIBRANT"
        if dirname.startswith("VIBRANT"):
            # Loop over all subdirectories in the current directory
            for subdirname in os.listdir(dirname):
                # Check if the subdirectory starts with "VIBRANT_phages"
                if subdirname.startswith("VIBRANT_phages"):
                    # Construct the path to the combined.fna file
                    fasta_filename = subdirname.replace("VIBRANT_phages_", "") + ".phages_combined.fna"
                    fasta_path = os.path.join(dirname, subdirname, fasta_filename)
                    # Open the file for reading
                    with open(fasta_path, "r") as infile:
                        # Loop over all lines in the file
                        for line in infile:
                            # Check if the line starts with ">"
                            if line.startswith(">"):
                                # Remove the ">" character from the Contig name
                                contig_name = line.strip()[1:]
                                # Write the directory name and contig name to the output file
                                outfile.write("{}\t{}\n".format(dirname, contig_name))
```

```bash
# Script 2 - Concatenate all integrated coordinates files (VIBRANT_integrated_prophage_coordinates_*.tsv)

cat ./VIBRANT_*/VIBRANT_results*/VIBRANT_integrated_prophage_coordinates_*.tsv/ > concatenated_prophage_coordinates_vibrant.tsv
```

```python
# Script 3 - Remove all sub-headings except the top row (02_remove_words_from_coordinates_file.py)

import csv

# Open the input and output files
with open('concatenated_prophage_coordinates_vibrant.tsv', 'r') as infile, open('concatenated_prophage_coordinates_vibrant_names_removed.tsv', 'w', newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')
    
    # Iterate over each row in the input file
    for row in reader:
        # Remove the columns by value
        row = [col for col in row if col not in ["scaffold", "fragment", "protein start", "protein stop", "protein length", "nucleotide start", "nucleotide stop", "nucleotide length"]]
        
        # Write the modified row to the output file
        writer.writerow(row)
```

```python
# Script 4 - 03_add_coordinates_to_prophage_fragments.py

import csv

# Read the output.tsv file
with open('output.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    output_data = [row for row in reader]

# Read the concatenated_prophage_coordinates_vibrant_names_removed.tsv file
with open('concatenated_prophage_coordinates_vibrant_names_removed.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    coord_data = [row for row in reader]

# Create a dictionary of fragment -> nucleotide_start/nucleotide_stop mappings
coord_dict = {}
for row in coord_data:
    coord_dict[row['fragment']] = (row['nucleotide_start'], row['nucleotide_stop'])

# Add nucleotide_start and nucleotide_stop columns to output_data if there is a 100% match
for row in output_data:
    contig = row['Contig']
    if contig in coord_dict:
        row['nucleotide_start'], row['nucleotide_stop'] = coord_dict[contig]

# Write the updated output_data to a new file
with open('output_with_coords.tsv', 'w') as f:
    fieldnames = list(output_data[0].keys()) + ['nucleotide_start', 'nucleotide_stop']
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(output_data)
```

```python
# Script 5 - 04_extract_contig_name_from_vibrant.py

import pandas as pd

# Read the TSV file into a pandas DataFrame
df = pd.read_csv("output_with_coords.tsv", sep="\t")

# Extract the contig name from the "Contig" column
df["contig_name"] = df["Contig"].str.extract(r"(.+)_fragment")

# If "_fragment" is missing, add the contents of the cell to the "contig_name" column
df.loc[df["contig_name"].isna(), "contig_name"] = df["Contig"]

# Save the updated DataFrame to a new TSV file
df.to_csv("output_with_contig_name.tsv", sep="\t", index=False)

# Print a message to confirm that the file was saved
print("Output saved to output_with_contig_name.tsv")
```

```python
# Script 6 - 05_remove_VIBRANT_from_lines.py

with open("output_with_contig_name_vibrant.tsv", "r") as f_in, open("output_with_contig_name_final.tsv", "w") as f_out:
    for line in f_in:
        if line.startswith("VIBRANT_"):
            line = line.replace("VIBRANT_", "", 1)
        f_out.write(line)
```

## # ** Please note ** - The coordinates file from VIBRANT would later be concatenated with the results from the other identification tools.

2) geNomad

```python
# Script 1 - extract_phage_info_tsv.py

import os

# Set the root directory to the current directory
root_dir = os.getcwd()

# Define the output file for the results
output_file = "output.tsv"

# Open the output file for writing
with open(output_file, "w") as f:
    # Write the header row to the output file
    f.write("Folder\tContig\n")

    # Loop through each folder in the root directory that starts with "genomad"
    for folder in os.listdir(root_dir):
        if folder.startswith("genomad"):
            # Loop through each subdirectory in the folder
            for sub_dir in os.listdir(os.path.join(root_dir, folder)):
                # Check if the subdirectory ends with "summary"
                if sub_dir.endswith("summary"):
                    # Define the path to the "virus.fna" file
                    virus_file = os.path.join(root_dir, folder, sub_dir, f"{sub_dir.split('_summary')[0]}_virus.fna")

                    # Open the file and loop through each line
                    with open(virus_file, "r") as fasta_file:
                        for line in fasta_file:
                            # If the line starts with ">", write the folder name and contig name to the output file
                            if line.startswith(">"):
                                contig_name = line.strip().lstrip(">")
                                f.write(f"{folder}\t{contig_name}\n")
```

```python

# Script 2 - add_prophage_coordinates_to_output.py

import csv

# Open the input and output files
with open('output.tsv', 'r') as input_file, open('output_with_start_end.tsv', 'w', newline='') as output_file:
    # Create a CSV reader and writer objects
    reader = csv.DictReader(input_file, delimiter='\t')
    writer = csv.DictWriter(output_file, fieldnames=reader.fieldnames + ['Start', 'End'], delimiter='\t')
    writer.writeheader()

    # Loop through each row in the input file
    for row in reader:
        # Search for the string "provirus_" in the "Contig" column
        if "provirus_" in row["Contig"]:
            # Extract the first and second numbers after "provirus_"
            start = int(row["Contig"].split("provirus_")[1].split("_")[0])
            end = int(row["Contig"].split("provirus_")[1].split("_")[1])

            # Add the "Start" and "End" columns to the row
            row["Start"] = start
            row["End"] = end

        # Write the row to the output file
        writer.writerow(row)
```

```python
# Script 3 - 03_get_contig_name_from_output_file.py

import csv

# Open the input and output files
with open('output_with_start_end.tsv', 'r') as input_file, open('output_with_contig_name.tsv', 'w', newline='') as output_file:
    # Create a CSV reader and writer objects
    reader = csv.DictReader(input_file, delimiter='\t')
    writer = csv.DictWriter(output_file, fieldnames=reader.fieldnames + ['contig_name'], delimiter='\t')
    writer.writeheader()

    # Loop through each row in the input file
    for row in reader:
        # Check if the "Contig" column contains a "|" symbol
        if '|' in row['Contig']:
            # If it does, add the string before the "|" symbol to the "contig_name" column
            row['contig_name'] = row['Contig'].split('|')[0]
        else:
            # If it doesn't, add the string up until the first space to the "contig_name" column
            row['contig_name'] = row['Contig'].split(' ')[0]

        # Write the updated row to the output file
        writer.writerow(row)
```

```python
# Script 4 - 04_remove_genomad_from_lines.py

with open("output_with_contig_name.tsv", "r") as f_in, open("output_with_contig_name_genomad_final.tsv", "w") as f_out:
    for line in f_in:
        if line.startswith("genomad_"):
            line = line.replace("genomad_", "", 1)
        f_out.write(line)
```

## # ** Please note ** - The coordinates file from geNomad would later be concatenated with the results from the other identification tools.

3) PhageBoost

```python
# Script 1 - extract_phage_info_tsv.py

import os

# Define the parent directory where the phageboost folders are located
parent_dir = os.getcwd()

# Open the output file for writing
with open("output.tsv", "w") as output_file:
    # Write the header row to the output file
    output_file.write("Parent Directory\tContig Name\n")

    # Loop through each folder in the parent directory
    for folder_name in os.listdir(parent_dir):
        # Check if the folder name starts with "phageboost"
        if folder_name.startswith("phageboost"):
            # Get the full path to the folder
            folder_path = os.path.join(parent_dir, folder_name)

            # Loop through each file in the folder
            for file_name in os.listdir(folder_path):
                # Check if the file ends with ".fasta"
                if file_name.endswith(".fasta"):
                    # Get the full path to the file
                    file_path = os.path.join(folder_path, file_name)

                    # Open the file for reading
                    with open(file_path, "r") as fasta_file:
                        # Read the first line of the file to get the contig name
                        contig_name = fasta_file.readline().strip()[1:]

                        # Write the parent directory and contig name to the output file
                        output_file.write(f"{folder_name}\t{contig_name}\n")
```

```python
# Script 2 - add_prophage_coordinates_to_output_2.py

import csv

with open('output.tsv', 'r') as tsvfile:
    reader = csv.DictReader(tsvfile, delimiter='\t')
    fieldnames = reader.fieldnames + ['Start', 'End']
    with open('output_with_start_end_2.tsv', 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in reader:
            if 'phage' in row['Contig']:
                if 'length' in row['Contig']:
                    start = row['Contig'].split('_')[7]
                    end = row['Contig'].split('_')[8]
                else:
                    start = row['Contig'].split('_')[2]
                    end = row['Contig'].split('_')[3]
                row['Start'] = start
                row['End'] = end
            writer.writerow(row)
```

```python
# Script 3 - extract_contig_name_from_output.py

import pandas as pd

# Read the TSV file into a pandas DataFrame
df = pd.read_csv('output_with_start_end_2.tsv', sep='\t')

# Extract the contig name from the "Contig" column
df['contig_name'] = df['Contig'].str.extract(r'(.+)_phage')

# Save the updated DataFrame to a new TSV file
df.to_csv('output_with_contig_name.tsv', sep='\t', index=False)

# Print a message to confirm that the file was saved
print('Output saved to output_with_contig_name.tsv')
```

```python
# Script 4 - 04_remove_PhageBoost_from_lines.py

with open("output_with_contig_name_phageboost.tsv", "r") as f_in, open("output_with_contig_name_phageboost_final.tsv", "w") as f_out:
    for line in f_in:
        if line.startswith("phageboost_"):
            line = line.replace("phageboost_", "", 1)
        f_out.write(line)
```

## # ** Please note ** - The coordinates file from PhageBoost would later be concatenated with the results from the other identification tools.

### Step 3 | Concatenation of prophage coordinates results

### Main Step 1

```bash
# Script 1 - Concatenate prophage coordinate files from each identification tool
# Repeats of headers were manually removed from output file

cat VIBRANT_coordinates.tsv geNomad_coordinates.tsv PhageBoost_coordinates.tsv > prophage_all_data_combined.tsv
```

### Main Step 2

```python
# Script 2 - Edit the concatenated file to include stop, start coordinates for empty cells.
# Where the prophage is an entire contig, there were no coordinates.
# '1' and 'end' were added as the start and end of the prophage.

with open("prophage_all_data_combined.tsv", "r") as f:
    lines = f.readlines()
    header = lines[0].strip().split("\t")
    nucleotide_start_index = header.index("nucleotide_start")
    nucleotide_stop_index = header.index("nucleotide_stop")
    new_lines = [lines[0]]
    for line in lines[1:]:
        fields = line.strip().split("\t")
        if fields[nucleotide_start_index] == "":
            fields[nucleotide_start_index] = "1"
        if fields[nucleotide_stop_index] == "":
            fields[nucleotide_stop_index] = "end"
        new_lines.append("\t".join(fields) + "\n")

with open("prophage_all_data_combined_modified.tsv", "w") as f:
    f.writelines(new_lines)
```

### Main Step 3

## # As there were variations to the names of the genome files when the original analysis was run (VIBRANT vs. geNomad & PhageBoost), extra steps were added to standardise the genome names for concatenation of prophage coordinates from each tool.

** **Please note** ** - The original files for each analysis were the same but issues with filename length arose, which forced us to change it for the original analysis.

Step 1: Extract the list of all the filenames that each prophage was identified in for VIBRANT (prophage_all_data_combined_modified.tsv file) and place in new file (new file = `vibrant_genome_filenames.txt`).

```bash
# vibrant_genome_filenames.txt

GCA_027686185.1_ASM2768618v1_genomic
GCA_027686185.1_ASM2768618v1_genomic
GCA_003437455.1_ASM343745v1_genomic
GCA_028205675.1_ASM2820567v1_genomic
GCA_014871875.1_ASM1487187v1_genomic
GCA_014871875.1_ASM1487187v1_genomic
GCA_014871875.1_ASM1487187v1_genomic
...
```

Step 2: Create list of all genome files used for geNomad and PhageBoost analysis (filename = `all_genomes_list.txt`)

Step 3: Use python script to change VIBRANT filenames to match those used for geNomad and PhageBoost (`add_genus_species_to_genome_files.py`). The output list (`output_genome_filenames_vibrant.txt`) was used to replace the genome filename for the VIBRANT results in the coordinates file.

```python
# Define the file paths
vibrant_file = "vibrant_genome_filenames.txt"
all_genomes_file = "all_genomes_list.txt"

# Read the content of vibrant_genome_filenames.txt
with open(vibrant_file, 'r') as vibrant_file_content:
    vibrant_lines = vibrant_file_content.readlines()

# Read the content of all_genomes_list.txt
with open(all_genomes_file, 'r') as all_genomes_content:
    all_genomes_lines = all_genomes_content.readlines()

# Create a dictionary to store the matches
matches = {}

# Iterate through each line in vibrant_genome_filenames.txt
for line in vibrant_lines:
    line = line.strip()  # Remove leading/trailing whitespace
    # Iterate through each line in all_genomes_list.txt
    for genome_line in all_genomes_lines:
        genome_line = genome_line.strip()
        # Check if the line from vibrant_genome_filenames.txt is a substring of the line in all_genomes_list.txt
        if line in genome_line:
            matches[line] = genome_line
            break  # Stop searching for this line once a match is found

# Create a new file with the matched lines as a new column
output_file = "output_genome_filenames_vibrant.txt"
with open(output_file, 'w') as output_content:
    for line in vibrant_lines:
        line = line.strip()
        # If a match was found, add it as a new column
        if line in matches:
            output_content.write(f"{line}\t{matches[line]}\n")
        else:
            output_content.write(f"{line}\n")

print(f"Matching lines have been saved to {output_file}")
```

Step 4: Contig names from each output were changed to match the original names from the input genome files. They were extracted and added to a new file called `contig_names.txt`. The string “fragment_…” was then removed from the end of each contig name, as these were added following the analysis (script 1 = `remove_fragment_from_contig_names.py`). 

```python
import re

# Function to process a single line
def process_line(line):
    return re.sub(r'_fragment_\d+', '', line)

# Input and output file names
input_file = "contig_names.txt"
output_file = "contig_names_processed.txt"

# Open input and output files
with open(input_file, "r") as input_file, open(output_file, "w") as output_file:
    for line in input_file:
        processed_line = process_line(line.strip())  # Strip to remove trailing newline
        output_file.write(processed_line + "\n")

print("Processing complete. The result is saved in", output_file)
```

Step 5: A sed was also used to extract the first part of the contig name for its unique identifier. This allowed for standardisation of contig names across the three tools and later excision of prophage from cognate host genome.

```bash
sed 's/ .*//' contig_names.txt > modified_contig_names.txt
```

Step 6: We noticed that there were quotation marks surrounding certain contig names, so these were removed from using a separate output file (`modified_contig_names.txt`). The script file was saved as `remove_quote_marks.py`.

```python
# Open the input file for reading
with open("modified_contig_names.txt", "r") as input_file:
    # Read the lines from the input file
    lines = input_file.readlines()

# Open the same file for writing (this will overwrite the original file)
with open("modified_contig_names.txt", "w") as output_file:
    for line in lines:
        # Remove the quotation marks from each line
        modified_line = line.replace('"', '')
        # Write the modified line back to the file
        output_file.write(modified_line)

print("Quotation marks removed from the file.")
```

Step 7: PhageBoost contig names were missing “.1” from the contig names for NCBI genomes. These were re-added using a Python script (`edit_phageboost_contig_names.py`).

```python
# Open the input file for reading and create the output file for writing
with open('phageboost_contig_names.txt', 'r') as input_file, open('output.txt', 'w') as output_file:
    for line in input_file:
        line = line.strip()  # Remove leading/trailing whitespace

        # Check conditions for appending ".1"
        if not (line.startswith("LH_") or line.startswith("contig") or line.startswith("NODE") or line.startswith(".26835")):
            line += ".1"

        # Write the modified line to the output file
        output_file.write(line + '\n')
```

Step 8: The outputs were concatenated manually and the columns rearranged to produce a tab-delimited file called `prophage_all_data_combined.txt` . It looks like this:

```
Genome_file	Contig_prophage_name	nucleotide_start	nucleotide_stop	contig_name	contig_name_edited	Tool
Bifidobacterium_adolescentis_GCA_027686185.1_ASM2768618v1_genomic	"JAQCTB010000003.1 Bifidobacterium adolescentis strain AF102-65 BF3HA05001.Scaf3, whole genome shotgun sequence_fragment_2"	90609	104081	"JAQCTB010000003.1 Bifidobacterium adolescentis strain AF102-65 BF3HA05001.Scaf3, whole genome shotgun sequence"	JAQCTB010000003.1	VIBRANT
Bifidobacterium_adolescentis_GCA_027686185.1_ASM2768618v1_genomic	"JAQCTB010000001.1 Bifidobacterium adolescentis strain AF102-65 BF3HA05001.Scaf1, whole genome shotgun sequence_fragment_10"	882381	922457	"JAQCTB010000001.1 Bifidobacterium adolescentis strain AF102-65 BF3HA05001.Scaf1, whole genome shotgun sequence"	JAQCTB010000001.1	VIBRANT
Bifidobacterium_longum_GCA_003437455.1_ASM343745v1_genomic	"QSQH01000001.1 Bifidobacterium longum strain TM01-1 TM01-1.Scaf1, whole genome shotgun sequence_fragment_3"	111615	129102	"QSQH01000001.1 Bifidobacterium longum strain TM01-1 TM01-1.Scaf1, whole genome shotgun sequence"	QSQH01000001.1	VIBRANT
Bifidobacterium_adolescentis_GCA_028205675.1_ASM2820567v1_genomic	"JAQKOF010000005.1 Bifidobacterium adolescentis strain D53st1_B10_D53t1_180928 NODE_5_length_92005_cov_31.4441, whole genome shotgun sequence_fragment_1"	17791	31152	"JAQKOF010000005.1 Bifidobacterium adolescentis strain D53st1_B10_D53t1_180928 NODE_5_length_92005_cov_31.4441, whole genome shotgun sequence"	JAQKOF010000005.1	VIBRANT
...
```

### Main Step 4

# Combining prophage coordinates from each tool and then excising the sequences from cognate bacterial genome.

Step 1: The previous `prophage_all_data_combined.txt` file was manually edited by changing the order of the columns. The new file was named `combined_prophage_data.txt`. It looks like this:

```
Genome_file	contig_name_edited	Contig_prophage_name	nucleotide_start	nucleotide_stop	Tool
Bifidobacterium_adolescentis_GCA_027686185.1_ASM2768618v1_genomic	JAQCTB010000003.1	JAQCTB010000003.1 Bifidobacterium adolescentis strain AF102-65 BF3HA05001.Scaf3, whole genome shotgun sequence_fragment_2	90609	104081	VIBRANT
Bifidobacterium_adolescentis_GCA_027686185.1_ASM2768618v1_genomic	JAQCTB010000001.1	JAQCTB010000001.1 Bifidobacterium adolescentis strain AF102-65 BF3HA05001.Scaf1, whole genome shotgun sequence_fragment_10	882381	922457	VIBRANT
Bifidobacterium_longum_GCA_003437455.1_ASM343745v1_genomic	QSQH01000001.1	QSQH01000001.1 Bifidobacterium longum strain TM01-1 TM01-1.Scaf1, whole genome shotgun sequence_fragment_3	111615	129102	VIBRANT
Bifidobacterium_adolescentis_GCA_028205675.1_ASM2820567v1_genomic	JAQKOF010000005.1	JAQKOF010000005.1 Bifidobacterium adolescentis strain D53st1_B10_D53t1_180928 NODE_5_length_92005_cov_31.4441, whole genome shotgun sequence_fragment_1	17791	31152	VIBRANT
...
```

Step 2: In the `combined_prophage_data.txt` file, the “Contig_prophage_name” column was altered to replace spaces with underscores to avoid any potential issues with downstream bioinformatics processes using the following Python script (`replace_spaces_with_underscores.py`).

```
# Define the input and output file names
input_file = "combined_prophage_data.txt"
output_file = "combined_prophage_data_updated.txt"

# Initialize an empty list to store the updated data
updated_data = []

# Open and read the input file
with open(input_file, "r") as infile:
    # Read and process the header line
    header = infile.readline().strip().split("\t")
    updated_data.append(header)

    # Process the rest of the lines
    for line in infile:
        # Split the line into columns based on tab delimiters
        columns = line.strip().split("\t")

        # Replace spaces with underscores in the "Contig_prophage_name" column
        contig_prophage_name = columns[2].replace(" ", "_")

        # Update the line with the modified "Contig_prophage_name"
        columns[2] = contig_prophage_name

        # Append the updated line to the data
        updated_data.append(columns)

# Write the updated data to the output file
with open(output_file, "w") as outfile:
    for line in updated_data:
        outfile.write("\t".join(line) + "\n")

print(f"Data has been updated and saved to {output_file}")
```

Step 3: The “1” and “end” in for the prophage coordinates that were an entire contig were changed to “NA” using Python scripts (`replace_end_with_NA.py` & `replace_1_with_NA.py`)

```python
# Define the input and output file names
input_file = "combined_prophage_data_updated.txt"
output_file = "combined_prophage_data_updated_2.txt"

# Initialize an empty list to store the updated data
updated_data = []

# Open and read the input file
with open(input_file, "r") as infile:
    # Read and process the header line
    header = infile.readline().strip().split("\t")
    updated_data.append(header)

    # Process the rest of the lines
    for line in infile:
        # Split the line into columns based on tab delimiters
        columns = line.strip().split("\t")

        # Replace "end" with "NA" in the "nucleotide_stop" column
        nucleotide_stop = columns[4].replace("end", "NA")

        # Update the line with the modified "nucleotide_stop"
        columns[4] = nucleotide_stop

        # Append the updated line to the data
        updated_data.append(columns)

# Write the updated data to the output file
with open(output_file, "w") as outfile:
    for line in updated_data:
        outfile.write("\t".join(line) + "\n")

print(f"Data has been updated and saved to {output_file}")
```

```python
# Define the input and output file names
input_file = "combined_prophage_data_updated_2.txt"
output_file = "modified_combined_prophage_data.txt"

# Open the input and output files
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Process each line in the input file
    for line in infile:
        # Split the line into columns using tab as the delimiter
        columns = line.strip().split('\t')
        
        # Check if the "nucleotide_stop" column contains "NA"
        if columns[4] == "NA":
            # Replace "1" with "NA" in the "nucleotide_start" column
            columns[3] = "NA"
        
        # Write the modified line to the output file
        outfile.write('\t'.join(columns) + '\n')

print("Script completed. Modified data saved to", output_file)
```

Step 4: Use Python script to merge overlapping coordinates of prophages from the same contig, so that the earliest start and and latest end are extracted (`obtain_highest_lowest_coordinates.py`). Input filename changed to `final_combined_prophage_data.txt`.

```python
import pandas as pd
import csv
from collections import defaultdict

# Load the tab-delimited file into a pandas DataFrame
df = pd.read_csv('final_combined_prophage_data.txt', sep='\t', quoting=csv.QUOTE_MINIMAL)

# Group by 'contig_name_edited'
grouped = df.groupby('contig_name_edited')

# Create empty lists to store the results
contig_data = defaultdict(list)

# Define a function to find the maximum overlap within a group
def find_max_overlap(group):
    contig_name_edited = group['contig_name_edited'].iloc[0]
    accession = group['Genome_file'].iloc[0]
    starts = list(group['nucleotide_start'])
    ends = list(group['nucleotide_stop'])

    # Check for missing values ('NA') in 'start' or 'end'
    if any(pd.isna(start) or pd.isna(end) for start, end in zip(starts, ends)):
        # If any 'NA' value is found, replace coordinates with 'NA'
        contig_data[contig_name_edited].append(('NA', 'NA', accession))
    else:
        # Sort the ranges by start value
        sorted_ranges = sorted(zip(starts, ends))

        # Initialize the current range
        current_start, current_end = sorted_ranges[0]

        # Iterate through the sorted ranges
        for start, end in sorted_ranges[1:]:
            if start <= current_end:
                # Overlapping ranges, update the current_end
                current_end = max(current_end, end)
            else:
                # Non-overlapping range found, append the current range to the results
                contig_data[contig_name_edited].append((current_start, current_end, accession))
                # Set the current range to the new range
                current_start, current_end = start, end

        # Append the final current range to the results
        contig_data[contig_name_edited].append((current_start, current_end, accession))

# Apply the find_max_overlap function to each group
grouped.apply(find_max_overlap)

# Create a new DataFrame with the results
result_data = []
for contig_name, data_list in contig_data.items():
    for start, end, accession in data_list:
        result_data.append((contig_name, start, end, accession))
result_df = pd.DataFrame(result_data, columns=['contig_name', 'new_start', 'new_end', 'accession'])

# Save the result DataFrame to a tab-delimited file
result_df.to_csv('final_prophage_coordinates.txt', sep='\t', index=False)

# Print the result DataFrame
print(result_df)
```

Step 5: Prophage sequences were extracted from their respective cognate host genome using a Python script (`extract_prophage_coordinates_from_genome.py` ).

```python
import os
from Bio import SeqIO

# Define the folder location for genome files
genome_folder = "path/to/relvant/folder/"

# Define the output directory for extracted sequences
output_directory = "path/to/relvant/folder/"

# Read the result2.txt file
result_file = "final_prophage_coordinates.txt"

# Create a dictionary to store the start and end coordinates for each contig and accession
coordinates = {}

with open(result_file, "r") as f:
    next(f)  # Skip the header line
    for line in f:
        contig_name, new_start, new_end, accession = line.strip().split("\t")
        # Extract the contig name up to the first space
        contig_name = contig_name.split()[0]
        coordinates[(accession, contig_name, new_start, new_end)] = (new_start, new_end)

# Process each accession and contig
for (accession, contig_name, new_start, new_end), (new_start, new_end) in coordinates.items():
    genome_file = os.path.join(genome_folder, f"{accession}.fna")
    
    if not os.path.exists(genome_file):
        print(f"Genome file not found for accession: {accession}")
        continue

    # Open the genome file using Biopython's SeqIO
    genome_sequences = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    if contig_name not in genome_sequences:
        print(f"Contig {contig_name} not found in genome file {accession}")
        continue

    sequence = genome_sequences[contig_name].seq

    if new_start != "NA" and new_end != "NA":
        start = int(float(new_start))
        end = int(float(new_end))
        sequence = sequence[start - 1:end]

    # Remove spaces and brackets from contig name
    contig_name = contig_name.replace(" ", "_").replace("(", "").replace(")", "")
    
    # Remove ".0" from coordinates in the contig name
    new_start = new_start.replace(".0", "")
    new_end = new_end.replace(".0", "")

    # Specify the full path for the output file with underscores between parts and no ".0" in coordinates
    output_filename = os.path.join(output_directory, f"{accession}_{contig_name}_{new_start}-{new_end}_sequence.fna")
    with open(output_filename, "w") as output_f:
        output_f.write(f">{accession}_{contig_name}_{new_start}-{new_end}\n{sequence}\n")

    print(f"Extracted sequence for {accession}_{contig_name} and saved to {output_filename}")
```

Step 6: Extract prophages from host genomes that were included following Seqfu filtering

```python
# Script 1:

import pandas as pd

# Read the tab-delimited file into a DataFrame
df = pd.read_csv('bacterial_metadata_file.tsv', sep='\t')

# Function to create the Accession_original column
def modify_accession(accession):
    if accession.startswith('GCA_'):
        parts = accession.rsplit('_', 1)
        return parts[0] + '.' + parts[1]
    return accession

# Apply the function to create the new column
df['Accession_original'] = df['Accession_2'].apply(modify_accession)

# Save the updated DataFrame to a new tab-delimited file
df.to_csv('metadata_updated.tsv', sep='\t', index=False)
```

```python
# Script 2:

import pandas as pd

# Read the updated metadata file
metadata_df = pd.read_csv('02_prophage_metadata_file_updated.tsv', sep='\t')

# Read the seqfu genomes list
with open('02_seqfu_genomes_list.txt', 'r') as file:
    seqfu_genomes_list = file.read().splitlines()

# Find matches and create a list of corresponding entries from Prophage_name_matched
matching_prophage_names = metadata_df[metadata_df['Accession_original'].isin(seqfu_genomes_list)]['Prophage_name_matched'].tolist()

# Print the list of matching Prophage_name_matched entries
print(matching_prophage_names)

# Optionally, save the matching Prophage_name_matched entries to a new file
with open('matching_prophage_names.txt', 'w') as output_file:
    for name in matching_prophage_names:
        output_file.write(name + '\n')
```

```python
# Script 3:

import os
import shutil

# Read the list of matching prophage names
with open('03_matching_prophage_names.txt', 'r') as file:
    matching_prophage_names = file.read().splitlines()

# Directory containing the .fna files
source_directory = '.'

# Destination directory
destination_directory = 'seqfu_phages'

# Ensure the destination directory exists
if not os.path.exists(destination_directory):
    os.makedirs(destination_directory)

# Iterate over files in the source directory
for file_name in os.listdir(source_directory):
    # Check if the file has a .fna extension
    if file_name.endswith('.fna'):
        # Check if the base name (without extension) matches any prophage name in the list
        base_name = os.path.splitext(file_name)[0]
        if base_name in matching_prophage_names:
            # Move the matching file to the destination directory
            shutil.move(os.path.join(source_directory, file_name), os.path.join(destination_directory, file_name))

print("Matching genomes have been moved to the seqfu_phages directory.")
```

# 07 | Dereplication of prophages and creation of database

## # Dereplication of prophages at species level (95% ANI) and characterisation of prophages.

### Step 1 | Concatenate prophage sequences into single fasta file

```bash
cat *.fna > bif_prophages_unfiltered_combined.fna
```

### Step 2 | Assess quality of viral sequences and trim contigs where applicable using CheckV (v1.0.1).

```bash
checkv end_to_end <input_file.fna> /path/to/output_directory -t 12
```

### Step 3 | Remove low quality fragments from genome files.

Script 1 - Filter quality summary file to include only “Medium-quality”, “High-quality” and “Complete” prophages (`quality_summary.tsv` file is from CheckV output).

```python
# Open the input file
with open('quality_summary.tsv', 'r') as input_file:
    lines = input_file.readlines()

# Open the output file for writing
with open('filtered_quality_data.tsv', 'w') as output_file:
    # Write the header to the output file
    output_file.write(lines[0])

    # Iterate over each line in the input file, starting from the second line
    for line in lines[1:]:
        # Split the line into columns
        columns = line.strip().split('\t')

        # Extract relevant columns
        contig_length = int(columns[1])
        checkv_quality = columns[9]

        # Check conditions for exporting the data
        if contig_length >= 5000 and checkv_quality in ['Medium-quality', 'High-quality', 'Complete']:
            # Write the line to the output file
            output_file.write(line)

print("Filtered data exported successfully.")

```

Scripts 2 & 3 - Created filtered files for phages labelled as “provirus: yes or no” (from `filtered_quality_data.tsv` file).

```python
# Open the input file
with open('filtered_quality_data.tsv', 'r') as input_file:
    lines = input_file.readlines()

# Open the output file for writing
with open('filtered_quality_data_no_provirus.tsv', 'w') as output_file:
    # Write the header to the output file
    output_file.write(lines[0])

    # Iterate over each line in the input file, starting from the second line
    for line in lines[1:]:
        # Split the line into columns
        columns = line.strip().split('\t')

        # Extract the value of the "provirus" column
        provirus_value = columns[2]

        # Check if the value is "No"
        if provirus_value == "No":
            # Write the line to the output file
            output_file.write(line)

print("Filtered data exported successfully.")

```

```python
# Open the input file
with open('filtered_quality_data.tsv', 'r') as input_file:
    lines = input_file.readlines()

# Open the output file for writing
with open('filtered_quality_data_just_provirus.tsv', 'w') as output_file:
    # Write the header to the output file
    output_file.write(lines[0])

    # Iterate over each line in the input file, starting from the second line
    for line in lines[1:]:
        # Split the line into columns
        columns = line.strip().split('\t')

        # Extract the value of the "provirus" column
        provirus_value = columns[2]

        # Check if the value is "Yes"
        if provirus_value == "Yes":
            # Write the line to the output file
            output_file.write(line)

print("Filtered data exported successfully.")

```

Scripts 4 & 5 - Using matched contig names, extract all prophages that are categorised as “Medium-Quality” or higher in the output `viruses.fna` and the `proviruses.fna` file.

```python
# Read contig names from filtered_quality_data_no_provirus.tsv
tsv_file = "filtered_quality_data_no_provirus.tsv"
contig_names = set()
with open(tsv_file, 'r') as tsv:
    next(tsv)  # Skip header
    for line in tsv:
        contig_id = line.split("\t")[0]
        contig_names.add(contig_id)

# Read contigs from viruses.fna and create matched_viruses.fna
fasta_file = "viruses.fna"
matched_fasta_file = "matched_viruses.fna"
matching_contigs = []
with open(fasta_file, 'r') as fasta, open(matched_fasta_file, 'w') as matched_fasta:
    contig_name = ""
    for line in fasta:
        if line.startswith(">"):
            contig_name = line.strip()[1:]
            if contig_name in contig_names:
                matching_contigs.append(contig_name)
                matched_fasta.write(line)
        elif contig_name in matching_contigs:
            matched_fasta.write(line)

print("Matching contigs extracted to", matched_fasta_file)

```

```python
# Function to parse contig names from contig_id column
def parse_contig_name(contig_id):
    return contig_id.split("_1 ")[0]

# Read filtered_quality_data_just_provirus.tsv file
data_file = "filtered_quality_data_just_provirus.tsv"
contig_names = set()  # Set to store unique contig names
with open(data_file, 'r') as f:
    next(f)  # Skip header
    for line in f:
        contig_id = line.split("\t")[0]
        contig_name = parse_contig_name(contig_id)
        contig_names.add(contig_name)

# Read proviruses.fna and create new fna file with matching contigs
input_fasta = "proviruses.fna"
output_fasta = "matched_proviruses.fna"
with open(input_fasta, 'r') as f_in, open(output_fasta, 'w') as f_out:
    current_contig = ""
    for line in f_in:
        if line.startswith(">"):
            current_contig = line.strip()[1:].split("_1 ")[0]  # Extract contig name
            if current_contig in contig_names:
                f_out.write(line)
                match = True
            else:
                match = False
        else:
            if match:
                f_out.write(line)

print("Matching contigs extracted to", output_fasta)

```

Script 6 - Concatenate each matched virus and provirus file

```bash
cat *.fna > checkv_filtered_prophages_combined.fna
```

Script 7 - List all sequences that are below 5000bp and then remove them from `checkv_filtered_prophages_combined.fna` file.

```python
import pandas as pd

# Read the TSV file into a DataFrame
df = pd.read_csv('quality_summary.tsv', delimiter='\t')

# Iterate through rows and print contig_id and numbers below 5000
for index, row in df.iterrows():
    if row['checkv_quality'] not in ['Low-quality', 'Not-determined']:
        if row['provirus'] == 'No' and row['contig_length'] < 5000:
            print(f"Strain: {row['contig_id']}, Contig Length below 5000: {row['contig_length']}")
        elif row['provirus'] == 'Yes' and row['proviral_length'] < 5000:
            print(f"Strain: {row['contig_id']}, Proviral Length below 5000: {row['proviral_length']}")

```

### Step 4 | Dereplicate phages at MIUViG recommended parameters (95% ANI + 85% AF) using supporting code from CheckV ([https://bitbucket.org/berkeleylab/checkv/src/master/](https://bitbucket.org/berkeleylab/checkv/src/master/)).

Script 1 - Create a blast+ database.

```bash
makeblastdb -in <my_seqs.fna> -dbtype nucl -out <my_db>
```

Script 2 - Use megablast from blast+ package to perform all-vs-all blastn of sequences.

```bash
blastn -query <my_seqs.fna> -db <my_db> -outfmt '6 std qlen slen' -max_target_seqs 10000 -o <my_blast.tsv> -num_threads 16
```

Script 3 - Calculate pairwise ANI by combining local alignments between sequence pairs using `anicalc.py` file.

```bash
anicalc.py -i <my_blast.tsv> -o <my_ani.tsv>
```

**Please note: The output file was used to calculate intergenomic similarities of all prophages prior to dereplication (see Section 08)**

Script 4 - Perform UCLUST-like clustering using the MIUViG recommended-parameters (95% ANI + 85% AF) with `aniclust.py` file.

```bash
aniclust.py --fna <my_seqs.fna> --ani <my_ani.tsv> --out <my_clusters.tsv> --min_ani 95 --min_tcov 85 --min_qcov 0

# The file <my_clusters.tsv> contains the clustering results. The first column is the cluster representative, and the second column contains cluster members. The file <my_ani.tsv> contains the pairwise ANI results. The columns are:

# query_id: query identifier
# target_id: checkv reference genome identifier
# alignment_count: number of blastn alignments
# ani: average nucleotide identity
# query_coverage: percent of query genome covered by alignments
# target_coverage: percent of target genome covered by alignments
```

### Step 5 | Extract cluster representative prophage genomes from clustering results file (`filt_bif_clusters.tsv`).

Script 1 - Extract representative sequences from original concatenated prophages fasta file following CheckV analysis.

```python
# Define the input and output file names
clusters_file = "filt_bif_clusters.tsv"
fasta_file = "checkv_filtered_prophages_combined_over_5000bp.fna"
output_file = "matched_contigs.fna"

# Create a set to store the contigs from the clusters file
contigs_set = set()

# Read contigs from clusters file and add to the set
with open(clusters_file, 'r') as clusters:
    for line in clusters:
        contig_name = line.strip()
        contigs_set.add(contig_name)

# Open output file for writing
with open(output_file, 'w') as output:
    # Read fasta file and write matching contigs to output file
    with open(fasta_file, 'r') as fasta:
        current_contig = ""
        for line in fasta:
            if line.startswith('>'):
                current_contig = line.strip()[1:]
                if current_contig in contigs_set:
                    output.write(line)
            else:
                if current_contig in contigs_set:
                    output.write(line)

print(f"Matching contigs written to {output_file}")
```

Script 2 - Standardise phage names to “Phage_1”, “Phage_2”, etc. (Applied to all viruses, including dereplicated dataset)

```python
import re

# Function to replace non-alphanumeric characters with underscores
def clean_contig_name(contig_name):
    return re.sub(r'[^a-zA-Z0-9]', '_', contig_name)

# Read concatenated_viruses_checkv.fna and create a new fasta file with cleaned contig names
input_fasta = "concatenated_viruses_checkv.fna"
output_fasta = "cleaned_concatenated_viruses_checkv.fna"

with open(input_fasta, 'r') as f_in, open(output_fasta, 'w') as f_out:
    for line in f_in:
        if line.startswith(">"):
            # Extract contig name and clean it
            contig_name = line.strip()[1:]
            cleaned_contig_name = clean_contig_name(contig_name)
            f_out.write(">" + cleaned_contig_name + "\n")
        else:
            f_out.write(line)

print("Contig names cleaned and saved to", output_fasta)
```

Script 3 - Create individual fasta files for each prophage

```python
# Read the cleaned concatenated fasta file
input_fasta = "cleaned_concatenated_viruses_checkv.fna"

with open(input_fasta, 'r') as f:
    current_contig_file = None
    for line in f:
        if line.startswith(">"):
            # Close the previous contig file if exists
            if current_contig_file:
                current_contig_file.close()
            # Extract contig name and open a new file
            contig_name = line.strip()[1:]
            current_contig_file = open(contig_name + ".fna", 'w')
            current_contig_file.write(line)  # Write contig name to new contig file
        else:
            # Write sequence to current contig file
            current_contig_file.write(line)

# Close the last contig file
if current_contig_file:
    current_contig_file.close()

print("Individual .fna files created for each contig.")
```

### The prophage database is now complete and ready for downstream analysis!

# 08 | Characterisation of prophages

### Plot intergenomic similarity using `my_ani.tsv` file from CheckV supporting code output with metadata (R script)

```r
library("reshape2")
library("ComplexHeatmap")
library("gplots")
library("heatmap3")
library("RColorBrewer")

organism_palette <- c(
  "Bifidobacterium scardovii" = "#d6859d",
  "Bifidobacterium longum" = "#a0da39",
  "Bifidobacterium adolescentis" = "#956cb4",
  "Bifidobacterium ruminantium" = "#ffa500",
  "Bifidobacterium piotii" = "#ffff00",
  "Bifidobacterium pseudolongum" = "#c44e52",
  "Bifidobacterium catenulatum" = "#6c92cc",
  "Bifidobacterium thermophilum" = "#ff0000",
  "Bifidobacterium gallicum" = "#277f8e",
  "Bifidobacterium pseudocatenulatum" = "#55a868",
  "Bifidobacterium breve" = "#94b0c2",
  "Bifidobacterium pullorum" = "#990000",
  "Bifidobacterium bifidum" = "#440154",
  "Bifidobacterium animalis" = "#bade28",
  "Bifidobacterium dentium" = "#bd8e46",
  "Bifidobacterium angulatum" = "#8c613c",
  "missing" = "lightgrey"
)

source_palette <- c(
  "milk" = "#D55E00",      # Color-blind friendly
  "eye" = "#0072B2",       # Color-blind friendly
  "dental caries" = "#E69F00", # Color-blind friendly
  "missing" = "lightgrey",
  "faeces" = "#009E73",    # Color-blind friendly
  "vagina" = "#CC79A7",    # Color-blind friendly
  "supplements" = "#F0E442", # Color-blind friendly
  "blood" = "#56B4E9",     # Color-blind friendly
  "urine" = "#0072B2"      # Color-blind friendly
)

age_palette <- c(
  "missing" = "lightgrey",
  "NA" = "lightgrey",
  "child" = "#56B4E9",     # Color-blind friendly
  "elderly" = "#E69F00",   # Color-blind friendly
  "infant" = "#CC79A7",    # Color-blind friendly
  "adult" = "#009E73"      # Color-blind friendly
)

country_palette <- c(
  "Italy" = "#1f78b4",
  "Netherlands" = "#33a02c",
  "South Africa" = "#e31a1c",
  "Japan" = "#ff7f00",
  "Czech Republic" = "#6a3d9a",
  "Mexico" = "#a6cee3",
  "USA" = "#b2df8a",
  "France" = "#fdbf6f",
  "Cambodia" = "#cab2d6",
  "Belgium" = "#fb9a99",
  "Canada" = "#ffff99",
  "Singapore" = "#8c510a",
  "Ireland" = "#d6859d",
  "Kenya" = "#addd8e",
  "South Korea" = "#a6a6a6",
  "Chile" = "#525252",
  "Brazil" = "#ff7f7f",
  "Mozambique" = "#a65628",
  "Germany" = "#984ea3",
  "Australia" = "#999999",
  "Russia" = "#2ca25f",
  "Thailand" = "#0570b0",
  "Vietnam" = "#fdae61",
  "Sweden" = "#fee08b",
  "China" = "#ffff00",
  "United Kingdom" = "#d73027",
  "India" = "#4575b4",
  "Spain" = "#91bfdb",
  "missing" = "lightgrey"
)

# Get data, convert to matrix
data <- read.csv("07_completed_ani_matrix.csv", header = TRUE, row.names = 1)
data_matrix <- as.matrix(data)

# Define a single gradient from 0 to 100
breaks <- seq(min(data_matrix), 100, length.out = 100)
gradient <- colorpanel(length(breaks) - 1, "white", "darkblue")

hm.colors <- gradient

# Define metadata
metadata <- read.table("05_updated_prophage_metadata_734.tsv", header = TRUE, sep = "\t", row.names = 1)

# Extract metadata corresponding to rows in data_matrix
metadata_rows <- metadata[rownames(data_matrix), c("organismName", "isolation_source_clean", "age", "country")]

# Create factors with specific levels
metadata_rows$organismName <- factor(metadata_rows$organismName, levels = names(organism_palette))
metadata_rows$isolation_source_clean <- factor(metadata_rows$isolation_source_clean, levels = names(source_palette))
metadata_rows$age <- factor(metadata_rows$age, levels = names(age_palette))
metadata_rows$country <- factor(metadata_rows$country, levels = names(country_palette))

# Create numeric vectors for color annotations
organism_numeric <- as.numeric(metadata_rows$organismName)
source_numeric <- as.numeric(metadata_rows$isolation_source_clean)
age_numeric <- as.numeric(metadata_rows$age)
country_numeric <- as.numeric(metadata_rows$country)

# Create a color matrix for RowSideColors
row_colors <- cbind(
  organismName = organism_palette[organism_numeric],
  isolation_source_clean = source_palette[source_numeric],
  age = age_palette[age_numeric],
  country = country_palette[country_numeric]
)

# Function to create legend
create_legend <- function(title, labels, colors) {
  legend <- Legend(at = labels, legend_gp = gpar(fill = colors), title = title)
  draw(legend, x = unit(1, "npc") - unit(2, "mm"), just = "right")
}

# Add "B. " to each entry in organism_palette and sort alphabetically
organism_palette <- setNames(organism_palette, paste0("B. ", names(organism_palette)))
organism_palette <- organism_palette[order(names(organism_palette))]

# Sort age_palette and source_palette alphabetically
age_palette <- age_palette[order(names(age_palette))]
source_palette <- source_palette[order(names(source_palette))]

# Create heatmap and legends in a single high-resolution PDF file
pdf("prophage_heatmap_legends_4Sep2024.pdf", width = 15, height = 10)  # Adjust size for compact layout

# Adjust layout to provide more space for legends
layout(matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(3, 1))

# Draw heatmap
heatmap3(data_matrix, scale = "none", trace = "none", col = hm.colors, cexRow = 0.50, cexCol = 0.50,
         key = TRUE, key.title = NA, key.xlab = "Color", key.ylab = NA, dendrogram = "both",
         Rowv = TRUE, Colv = TRUE, symm = TRUE, density.info = "none", RowSideColors = row_colors)

# Draw legends
pushViewport(viewport(layout = grid.layout(3, 1)))

pushViewport(viewport(layout.pos.row = 1))
create_legend("Organism Palette", names(organism_palette), organism_palette)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
create_legend("Age Palette", names(age_palette), age_palette)
popViewport()

pushViewport(viewport(layout.pos.row = 3))
create_legend("Source Palette", names(source_palette), source_palette)
popViewport()

popViewport()

dev.off()
```

### Command to run Pharokka

```bash
for X in `cat phage_list.txt`; do
    input_file="${X}"
    base_name=$(basename "${input_file}" .fna)
    pharokka.py -i "${input_file}" -o "pharokka_${base_name}" -d ~/pharokka_databases -t 12 -p "${base_name}"
done

```

### Command to run ViPTree

```bash
# Run on HPC

#!/bin/bash
#BATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -p XXX-XXX
#SBATCH --mem XG
#SBATCH -t 1-00:00

source package <enter code>

cd /path/to/relevant/directory/

ViPTreeGen ./all_bif_refseq_phages.fna ./viptree_refseq_phages --ncpus 64
```

### Commands to run AcaFinder, dbAPIS and DefenseFinder (—antidefensefinder module)

```bash
# AcaFinder (run on HPC)

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -p XXX-XXX
##SBATCH --mem XG
#SBATCH -t 1-00:00

# Load necessary modules or source the environment
source package <enter code>

# Set the directory containing the Phage_ folders
parent_dir="./"

# Loop through each folder starting with "Phage_"
for folder in "$parent_dir"/Phage_*; do
    # Check if it's a directory
    if [ -d "$folder" ]; then
        # Get the basename of the folder to use in file paths
        base_name=$(basename "$folder")
        
        # Get the path to the .fna, .faa, and .gff files
        fna_file="$folder/${base_name}.fna"
        faa_file="$folder/phanotate.faa"
        gff_file="$folder/${base_name}.gff"

        # Check if the files exist
        if [ -f "$fna_file" ] && [ -f "$faa_file" ] && [ -f "$gff_file" ]; then
            # Run the command
            AcaFind_runner.py --Virus -z /path/to/all_pFam_hmm/ --FNA_file "$fna_file" -o "$folder/${base_name}/" --threads 64
        else
            echo "Error: Missing files in $folder"
        fi
    fi
done
```

```bash
# dbAPIS (run on VM)

#!/bin/bash

# Set the path to your dbAPIS installation
DBAPIS_PATH=/path/to/dbAPIS

# Set the path to your input .faa files
INPUT_DIR=/path/to/faa/files

# Set the path to store the results
OUTPUT_DIR=/path/to/output/directory
mkdir -p $OUTPUT_DIR

# Loop through all .faa files in the input directory
for file in $INPUT_DIR/*.faa; do
    # Extract file name without extension
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*}"

    # Step 1: Run hmmscan
    hmmscan --domtblout $OUTPUT_DIR/hmmscan_$filename_no_ext.out --noali $DBAPIS_PATH/dbAPIS.hmm $file

    # Step 2: Run diamond
    diamond blastp --db $DBAPIS_PATH/APIS_db -q $file -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen -o $OUTPUT_DIR/diamond_$filename_no_ext.out --max-target-seqs 10000

    # Step 3: Parse annotation output
    bash $DBAPIS_PATH/parse_annotation_result.sh $OUTPUT_DIR/hmmscan_$filename_no_ext.out $OUTPUT_DIR/diamond_$filename_no_ext.out

    echo "Processed $filename"
done
```

```bash
# DefenseFinder `--antidefensefinder` module (run on VM)

#!/bin/bash

# Loop through each directory
for dir in */; do
  # Remove trailing slash from directory name
  dir=${dir%/}
  
  # Construct the .faa file path
  faa_file="${dir}/${dir}.faa"
  
  # Run the defense-finder command if the .faa file exists
  if [[ -f "$faa_file" ]]; then
    echo "Running defense-finder for $faa_file"
    defense-finder run "$faa_file" --out-dir "defensefinder_$dir" --antidefensefinder
  else
    echo "File $faa_file not found. Skipping."
  fi
done
```

### geNomad classification to determine if prophages were Caudoviricetes

```bash
# Used `annotate` module
# Run on VM

genomad annotate metagenome.fna genomad_output /path/to/genomad_db
```

### Attempts to determine genus- and species-level classification of prophages using **taxmyPHAGE**

```python
import os
import subprocess

def run_taxmyphage_in_phage_folders():
    # Get the current working directory
    base_dir = os.getcwd()
    
    # List all directories in the current directory
    directories = [d for d in os.listdir(base_dir) if os.path.isdir(d) and d.startswith("Phage_")]
    
    for dir_name in directories:
        # Change to the directory
        os.chdir(dir_name)
        
        # Construct the command
        command = "taxmyphage run -i *.fna -t 12"
        
        # Run the command
        subprocess.run(command, shell=True)
        
        # Change back to the base directory
        os.chdir(base_dir)

# Run the function
run_taxmyphage_in_phage_folders()
```