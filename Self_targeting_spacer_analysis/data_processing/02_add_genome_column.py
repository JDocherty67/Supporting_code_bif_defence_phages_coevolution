import pandas as pd
import re

# Load the file
input_file = "02_filtered_100_pident.tsv"
output_file = "03_with_genome_column.tsv"

# Read the TSV file
df = pd.read_csv(input_file, sep='\t')

# Function to extract the genome name from qid
def extract_genome(qid):
    if qid.startswith("GCA"):
        # Match pattern like GCA_000154085_1_1_1 and convert last _ to .
        match = re.match(r"^(GCA_\d+_\d+)_\d+", qid)
        return match.group(1).replace("_", ".", 2) if match else qid
    else:
        return qid.split("_")[0]

# Apply the function to the 'qid' column
df['Genome'] = df['qid'].apply(extract_genome)

# Save to a new TSV file
df.to_csv(output_file, sep='\t', index=False)

print(f"Updated data with 'Genome' column saved to {output_file}")
