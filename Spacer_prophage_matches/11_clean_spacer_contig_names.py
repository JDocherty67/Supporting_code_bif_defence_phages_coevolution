import pandas as pd

# Load the TSV file
file_path = "10_final_deduplicated_self_targeting_analysis.tsv"
df = pd.read_csv(file_path, sep="\t")

# Clean the 'Spacers_cleaned' column by removing ':' and the numbers following it
df['Spacers_cleaned'] = df['Spacers_cleaned'].str.replace(r':\d+$', '', regex=True)

# Save the cleaned data to a new TSV file
output_file = "11_final_deduplicated_self_targeting_analysis_cleaned.tsv"
df.to_csv(output_file, sep="\t", index=False)

print(f"Cleaned data saved to {output_file}")

