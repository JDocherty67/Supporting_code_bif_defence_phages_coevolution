#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -p X-node
#SBATCH --mem=XG
#SBATCH -t 1-00:00
#SBATCH -o logs/blast_only_%j.out
#SBATCH -e logs/blast_only_%j.err

module load blast+

mkdir -p tmp_spacers blast_outputs logs

echo "Splitting spacers by genome..."

awk '
BEGIN {FS=":"}
/^>/ {
    header = $0;
    seq = "";
    getline seq_line;

    if (match(header, /^>GCA_[0-9]+_[0-9]+_[0-9]+/)) {
        split(substr(header, 2), a, "_");
        genome_id = a[1]"_"a[2]"."a[3];  # Replace third underscore with dot
    } else if (match(header, /^>LHP[0-9]+/)) {
        split(substr(header, 2), a, "_");
        genome_id = a[1];
    } else if (match(header, /^>C[0-9]+/)) {
        split(substr(header, 2), a, "_");
        genome_id = a[1];
    } else if (match(header, /^>SRR[0-9]+/)) {
        split(substr(header, 2), a, "_");
        genome_id = a[1];
    } else {
        next;
    }

    out = "tmp_spacers/"genome_id"_spacers.fa";
    print header >> out;
    print seq_line >> out;
}
' 05_filtered_spacers.fa

echo "Starting BLAST runs..."

for spacer_file in tmp_spacers/*_spacers.fa; do
    genome_id=$(basename "$spacer_file" | sed 's/_spacers.fa//')
    genome_file="${genome_id}.fna"
    output_file="blast_outputs/${genome_id}_blast.tsv"

    if [[ ! -f "$genome_file" ]]; then
        echo "Genome file $genome_file not found. Skipping $genome_id."
        continue
    fi

    echo "Running BLAST for $genome_id..."

    makeblastdb -in "$genome_file" -dbtype nucl

    blastn -query "$spacer_file" \
           -db "$genome_file" \
           -evalue 1e-5 \
           -task blastn-short \
           -outfmt "6 qseqid sseqid sstart send sstrand pident length evalue" \
           -max_target_seqs 100 \
           -num_threads 16 \
           -out "$output_file"

    echo "Finished $genome_id"
done

echo "All BLAST jobs complete. Output in blast_outputs/"
