# Open and read the file
input_file = '01_spacepharer_output.tsv'
output_file = '02_extracted_entries.tsv'

# Open the input file and output file
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        # Check if the line starts with '>'
        if line.startswith('>'):
            outfile.write(line)  # Write the line to the output file
            
print(f"Extraction complete. Results saved in {output_file}")
