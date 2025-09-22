import pandas as pd
import os
import re

# Input and output files
input_file = "03_with_genome_column.tsv"
gff_dir = "prokka_gffs"
output_file = "04_with_faa_product.tsv"

# Read the main table
df = pd.read_csv(input_file, sep='\t')

# Add empty columns for new data
df['faa_ID'] = None
df['product'] = None

# Process each genome group
for genome, group_df in df.groupby('Genome'):
    gff_file = os.path.join(gff_dir, f"{genome}.gff")
    if not os.path.isfile(gff_file):
        print(f"WARNING: GFF file not found for {genome}")
        continue

    # Parse the GFF file
    gff_entries = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            fields = line.strip().split('\t')
            contig = fields[0]
            attributes = fields[8]

            # Extract ID and product from attributes
            id_match = re.search(r'ID=([^;]+)', attributes)
            product_match = re.search(r'product=([^;]+)', attributes)

            if id_match:
                gff_entries.append({
                    'contig': contig,
                    'faa_ID': id_match.group(1),
                    'product': product_match.group(1) if product_match else None
                })

    # Organize GFF data by contig
    contig_dict = {}
    for entry in gff_entries:
        contig_dict.setdefault(entry['contig'], []).append(entry)

    # Match each row in the dataframe
    for idx, row in df[df['Genome'] == genome].iterrows():
        sseqid = row['sseqid']
        suffix_match = re.search(r'_(\d+)$', str(sseqid))
        if not suffix_match:
            continue
        contig_suffix = suffix_match.group(1)

        # Find GFF contig with matching suffix
        matched_contigs = [c for c in contig_dict if c.endswith(f"_{contig_suffix}")]
        if not matched_contigs:
            continue

        # Use the first matching CDS for the contig
        contig = matched_contigs[0]
        entry = contig_dict[contig][0]
        df.at[idx, 'faa_ID'] = entry['faa_ID']
        df.at[idx, 'product'] = entry['product']

# Save the updated file
df.to_csv(output_file, sep='\t', index=False)
print(f"Output written to {output_file}")
