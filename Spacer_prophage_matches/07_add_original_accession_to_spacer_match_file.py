import pandas as pd
import re

# Load the analysis file and metadata file
analysis_file = "06_self_targeting_analysis_with_species.tsv"
metadata_file = "04_updated_with_gc_metadata.tsv"

analysis_data = pd.read_csv(analysis_file, sep='\t')
metadata = pd.read_csv(metadata_file, sep='\t')

# Merge to add "Accession" based on matching "Genome" with "Accession_2"
merged_data = analysis_data.merge(
    metadata[['Accession_2', 'Accession_original']],
    left_on='Genome',
    right_on='Accession_2',
    how='left'
)

# Define a function to handle non-matched entries
def fill_accession(row):
    if pd.notna(row['Accession_original']):
        return row['Accession_original']
    else:
        # Modify "Genome" entry if it starts with "GCA_"
        if row['Genome'].startswith("GCA_"):
            return re.sub(r"_(\d+)$", r".\1", row['Genome'])
        else:
            return row['Genome']

# Apply the function to populate the "Accession" column
merged_data['Accession'] = merged_data.apply(fill_accession, axis=1)

# Drop unnecessary columns
merged_data = merged_data.drop(columns=['Accession_2', 'Accession_original'])

# Save the result to a new TSV file
merged_data.to_csv("07_self_targeting_analysis_with_accession_amended.tsv", sep='\t', index=False)

print("Updated 'Accession' column added. Results saved to '06_self_targeting_analysis_with_accession_amended.tsv'.")

