import pandas as pd

# Load the analysis file and the new metadata file
analysis_file = "07_self_targeting_analysis_with_accession_amended.tsv"
metadata_file = "10_updated_metadata_species_updated.tsv"

analysis_data = pd.read_csv(analysis_file, sep='\t')
metadata = pd.read_csv(metadata_file, sep='\t')

# Merge data to add "Bif_species" by matching "Accession" columns
merged_data = analysis_data.merge(
    metadata[['Accession', 'organismName']],
    on='Accession',
    how='left'
).rename(columns={'organismName': 'Bif_species'})

# Save the result to a new TSV file
merged_data.to_csv("08_self_targeting_analysis_with_bif_species.tsv", sep='\t', index=False)

print("Bif_species column added. Results saved to '07_self_targeting_analysis_with_bif_species.tsv'.")

