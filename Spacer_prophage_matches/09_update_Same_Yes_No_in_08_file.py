import pandas as pd

# Load the analysis file
analysis_file = "08_self_targeting_analysis_with_bif_species.tsv"
analysis_data = pd.read_csv(analysis_file, sep='\t')

# Update the 'Same' column based on comparison of 'Phage_host_species' and 'Bif_species'
analysis_data['Same'] = analysis_data.apply(
    lambda row: 'Yes' if row['Phage_host_species'] == row['Bif_species'] else 'No',
    axis=1
)

# Save the updated data to a new TSV file
analysis_data.to_csv("08_self_targeting_analysis_updated_same.tsv", sep='\t', index=False)

print("Updated 'Same' column based on species comparison. Results saved to '08_self_targeting_analysis_updated_same.tsv'.")

