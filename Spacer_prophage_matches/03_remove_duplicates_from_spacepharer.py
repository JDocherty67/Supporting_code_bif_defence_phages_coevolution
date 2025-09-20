import pandas as pd

# Load the TSV file
file_path = "03_modified_entries.tsv"
data = pd.read_csv(file_path, sep='\t')

# Sort by e-value (ascending) to prioritize rows with lower e-values
data['e-value'] = data['e-value'].astype(float)
data_sorted = data.sort_values(by='e-value', ascending=True)

# Drop duplicates based on "Spacer" and "Phage" while keeping the entry with the lowest e-value
cleaned_data = data_sorted.drop_duplicates(subset=['Spacer', 'Phage'], keep='first')

# Save the cleaned data to a new TSV file
cleaned_data.to_csv("04_cleaned_entries.tsv", sep='\t', index=False)

print("Duplicates removed based on e-value. Cleaned data saved to '04_cleaned_entries.tsv'.")

