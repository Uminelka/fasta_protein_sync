from Bio import SeqIO

fasta_file = 'organism_9606.fasta'

sizes_sequences = [100, 1000, 10000]

sequences = []

with open(fasta_file) as f:
    for sequence in SeqIO.parse(f, 'fasta'):
        sequences.append(sequence)

        if len(sequences) >= 10000:
            break

for size in sizes_sequences:
    fasta_files = f"example_{size}.fasta"
    
    with open(fasta_files, "w") as f:
        SeqIO.write(sequences[:size], f, "fasta")
        
    print(f, size)