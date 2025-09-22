import pandas as pd

# Load the TSV file
input_file = "00_filtered_blast_outputs.tsv"
output_file = "02_filtered_100_pident.tsv"

# Read the TSV file
df = pd.read_csv(input_file, sep='\t')

# Drop rows that are entirely empty
df.dropna(how='all', inplace=True)

# Filter rows where pident is exactly 100.000
df_filtered = df[df['pident'] == 100.000]

# Save to a new TSV file
df_filtered.to_csv(output_file, sep='\t', index=False)

print(f"Filtered data saved to {output_file}")
