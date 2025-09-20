import csv

# Input and output file names
input_file = '02_extracted_entries.tsv'
output_file = '03_modified_entries.tsv'

# Function to extract information based on the first column rules
def extract_accession(entry):
    # Remove the ">" symbol
    entry = entry[1:]
    
    if entry.startswith("GCA_"):
        return "_".join(entry.split("_")[:3])  # Up to the third underscore
    elif entry.startswith("C") or entry.startswith("SRR") or entry.startswith("LHP"):
        return entry.split("_")[0]  # Before the first underscore
    elif entry.startswith("LH_"):
        return "_".join(entry.split("_")[:2])  # Before the second underscore
    else:
        return entry  # Default: keep the full entry if no pattern matches

# Open the input and output files
with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    tsv_reader = csv.reader(infile, delimiter='\t')
    tsv_writer = csv.writer(outfile, delimiter='\t')
    
    for row in tsv_reader:
        if row:  # Check if the row is not empty
            # Extract accession and append as a new column
            accession = extract_accession(row[0])
            tsv_writer.writerow(row + [accession])  # Append the new column

print(f"Transformation complete. Results saved in {output_file}")

