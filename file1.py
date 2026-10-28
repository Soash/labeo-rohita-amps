from Bio import SeqIO

input_fasta = "Labeo_rohita_immune.fasta"
output_fasta = "Labeo_rohita_AMPs.fasta"

valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
min_length = 15
max_length = 50

unique_fragments = set()

with open(output_fasta, "w") as out_file:
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequence = str(record.seq).upper()

        # Sliding window fragmentation — every possible fragment length 15 to 50
        for i in range(len(sequence)):
            for j in range(i + min_length, min(i + max_length + 1, len(sequence) + 1)):
                fragment = sequence[i:j]

                if set(fragment).issubset(valid_amino_acids) and fragment not in unique_fragments:
                    unique_fragments.add(fragment)
                    out_file.write(f">rohu_fragment_{len(unique_fragments)}\n{fragment}\n")

print(f"Generated {len(unique_fragments)} unique fragments (lengths 15-50)!")



