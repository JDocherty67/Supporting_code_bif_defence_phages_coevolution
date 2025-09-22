import pandas as pd

# Load the file
input_file = "05_with_defense_info.tsv"
output_file = "06_with_STSS_final.tsv"

# Read the data
df = pd.read_csv(input_file, sep='\t')

# Create STSS_final column
df['STSS_final'] = df.apply(
    lambda row: row['defense_profile_name'] if pd.notna(row['defense_profile_name']) and str(row['defense_profile_name']).strip() != '' else row['product'],
    axis=1
)

# Save the updated file
df.to_csv(output_file, sep='\t', index=False)
print(f"Updated file with STSS_final column saved as {output_file}")
