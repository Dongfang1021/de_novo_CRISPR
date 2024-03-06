from datasketch import MinHash, MinHashLSH

# Parameters
num_perm = 128  # Number of permutations for MinHash
threshold = 0.8  # Similarity threshold

# Create an LSH index
lsh = MinHashLSH(threshold=threshold, num_perm=num_perm)

# Function to convert sequences to MinHash objects
def sequence_to_minhash(sequence, num_perm=128):
    m = MinHash(num_perm=num_perm)
    for d in sequence:
        m.update(d.encode('utf8'))
    return m

# Add sequences to LSH
for idx, seq in enumerate(sequences):
    lsh.insert(f"seq{idx}", sequence_to_minhash(seq, num_perm=num_perm))

# Query for similar pairs (example for one sequence)
query_seq = sequence_to_minhash(sequences[0], num_perm=num_perm)
similar_sequences = lsh.query(query_seq)

print("Similar sequences:", similar_sequences)
