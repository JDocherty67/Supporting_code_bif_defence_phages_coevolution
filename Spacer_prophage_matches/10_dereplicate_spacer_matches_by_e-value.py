import pandas as pd

# Load the analysis file
analysis_file = "08_self_targeting_analysis_updated_same.tsv"
analysis_data = pd.read_csv(analysis_file, sep='\t')

# Convert e-value to numeric for comparison
analysis_data['e-value'] = analysis_data['e-value'].astype(float)

# Sort by 'Spacer', 'Phage', and then by 'e-value' in ascending order
analysis_data = analysis_data.sort_values(by=['Spacer', 'Phage', 'e-value'])

# Drop duplicates by keeping the lowest e-value entry for each 'Spacer'-'Phage' combination
deduplicated_data = analysis_data.drop_duplicates(subset=['Spacer', 'Phage'], keep='first')

# Save the deduplicated data to a new TSV file
deduplicated_data.to_csv("10_final_deduplicated_self_targeting_analysis.tsv", sep='\t', index=False)

print("Dereplication complete. Results saved to '08_final_deduplicated_self_targeting_analysis.tsv'.")
