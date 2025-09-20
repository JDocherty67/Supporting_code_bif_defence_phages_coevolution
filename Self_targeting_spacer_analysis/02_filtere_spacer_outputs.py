import pandas as pd
import os
import glob

# === Load CRISPR metadata ===
crispr_file = "06_concatenated_crisprs_all_with_genome.tsv"
crispr_df = pd.read_csv(crispr_file, sep="\t", dtype=str)

# Ensure numeric types for coordinates
crispr_df["Start"] = crispr_df["Start"].astype(int)
crispr_df["End"] = crispr_df["End"].astype(int)

# Build lookup table: {(genome, contig): [(start, end), ...]}
region_dict = {}
for _, row in crispr_df.iterrows():
    key = (row["Genome"], row["Contig"])
    region_dict.setdefault(key, []).append((row["Start"], row["End"]))

# === Process each BLAST output ===
blast_files = glob.glob("blast_outputs/*_blast.tsv")

for blast_file in blast_files:
    genome_file_id = os.path.basename(blast_file).replace("_blast.tsv", "")  # e.g. GCA_028201875.1

    filtered_lines = []

    with open(blast_file) as f:
        for line in f:
            qid, sseqid, sstart, send, sstrand, pident, length, evalue = line.strip().split("\t")
            sstart, send = int(sstart), int(send)
            sstart, send = min(sstart, send), max(sstart, send)

            # === Map genome and contig to CRISPR metadata format ===
            if genome_file_id.startswith("GCA_") and "." in genome_file_id:
                # GCA_028201875.1 → GCA_028201875_1
                genome_lookup = genome_file_id.replace(".", "_")
                # Contig: first 4 underscore-separated parts
                contig_lookup = "_".join(sseqid.split("_")[:4])
            else:
                genome_lookup = genome_file_id
                # Contig: first 2 parts for LHP, C, SRR, etc.
                contig_lookup = "_".join(sseqid.split("_")[:2])

            key = (genome_lookup, contig_lookup)
            regions = region_dict.get(key, [])

            # === Check for overlap with any CRISPR region ===
            overlap = any(sstart <= end and send >= start for start, end in regions)

            if not overlap:
                filtered_lines.append(line.strip())

    # Write filtered result
    output_file = blast_file.replace("_blast.tsv", "_blast_filtered.tsv")
    with open(output_file, "w") as out:
        out.write("\n".join(filtered_lines) + "\n")

    print(f"✅ {os.path.basename(blast_file)} → {os.path.basename(output_file)} ({len(filtered_lines)} hits retained)")
