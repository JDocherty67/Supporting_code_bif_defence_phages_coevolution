import pandas as pd

# Load the files
file1 = "11_final_deduplicated_self_targeting_analysis_cleaned.tsv"
file2 = "12_concatenated_crisprs_all_with_genome.tsv"

df1 = pd.read_csv(file1, sep='\t')
df2 = pd.read_csv(file2, sep='\t')

# Merge data based on matches between 'Spacers_cleaned' and 'CRISPR'
merged_df = df1.merge(
    df2[['CRISPR', 'Prediction', 'Subtype', 'Subtype_probability']],
    left_on='Spacers_cleaned', 
    right_on='CRISPR', 
    how='left'
)

# Drop the extra 'CRISPR' column added by merge
merged_df = merged_df.drop(columns=['CRISPR'])

# Save the merged data to a new file
output_file = "12_final_deduplicated_self_targeting_analysis_cleaned_with_crispr_info.tsv"
merged_df.to_csv(output_file, sep='\t', index=False)

print(f"File saved as {output_file}")

