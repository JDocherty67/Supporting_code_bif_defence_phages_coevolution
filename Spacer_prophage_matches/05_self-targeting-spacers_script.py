import pandas as pd

# Load the main file and metadata file
cleaned_entries_file = "04_cleaned_entries.tsv"
metadata_file = "04_updated_with_gc_metadata.tsv"

cleaned_data = pd.read_csv(cleaned_entries_file, sep='\t')
metadata = pd.read_csv(metadata_file, sep='\t')

# Merge data on 'Phage' and 'Genome' with 'Prophage_name_matched' and 'Accession_2' to identify self-targeting spacers
merged_data = cleaned_data.merge(
    metadata,
    left_on=['Phage', 'Genome'],
    right_on=['Prophage_name_matched', 'Accession_2'],
    how='left',
    indicator='self_targeting'
)

# Mark as self-targeting if there is a match in the metadata (merge result is 'both')
merged_data['self_targeting'] = merged_data['self_targeting'].apply(lambda x: 'Yes' if x == 'both' else 'No')

# Select relevant columns, including the new 'self_targeting' column
output_data = merged_data[[
    'Genome', 'Spacer', 'Phage', 'e-value', 'spacer_start', 'spacer_end', 
    'phage_start', 'phage_end', 'PAM', 'PAM_reverse', 'self_targeting'
]]

# Save the result to a new TSV file
output_data.to_csv("05_self_targeting_analysis.tsv", sep='\t', index=False)

print("Self-targeting analysis complete. Results saved to '04_self_targeting_analysis.tsv'.")

