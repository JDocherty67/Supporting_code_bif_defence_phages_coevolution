import pandas as pd

# Input files
faa_file = "04_with_faa_product.tsv"
defense_file = "04_final_deduplicated_defense_systems.tsv"
output_file = "05_with_defense_info.tsv"

# Load both files
faa_df = pd.read_csv(faa_file, sep="\t")
defense_df = pd.read_csv(defense_file, sep="\t")

# Add new columns to faa_df
faa_df['defense_type'] = None
faa_df['defense_subtype'] = None
faa_df['defense_profile_name'] = None

# Iterate through each row in the faa_df
for idx, row in faa_df.iterrows():
    genome = row['Genome']
    faa_id = row['faa_ID']

    # Filter defense systems for the same genome
    matches = defense_df[defense_df['Genome'] == genome]

    for _, defense_row in matches.iterrows():
        protein_list = str(defense_row['protein_in_syst']).split(",")
        
        if faa_id in protein_list:
            index = protein_list.index(faa_id)
            profile_list = str(defense_row['name_of_profiles_in_sys']).split(",")

            faa_df.at[idx, 'defense_type'] = defense_row['type']
            faa_df.at[idx, 'defense_subtype'] = defense_row['subtype']
            # Defensive check in case profile count doesn't match
            if index < len(profile_list):
                faa_df.at[idx, 'defense_profile_name'] = profile_list[index]
            else:
                faa_df.at[idx, 'defense_profile_name'] = profile_list[0]  # Fallback to first
            break  # Stop after first match

# Save to new file
faa_df.to_csv(output_file, sep='\t', index=False)
print(f"Output written to {output_file}")

