import random

def random_dna(length):
    return ''.join(random.choice("ACGT") for _ in range(length))

def reverse_complement(seq):
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]


def make_transposon():
    IR_left = random_dna(random.randint(6, 10))
    IR_right = reverse_complement(IR_left)

    core = random_dna(random.randint(20, 40))

    DR = random_dna(random.randint(3, 6))

    transposon_sequence = DR + IR_left + core + IR_right + DR
    return transposon_sequence, IR_left, IR_right, DR


def insert_transposons(base_seq, num_TEs=3):
    positions = []
    seq = base_seq

    for _ in range(num_TEs):
        transposon, IR_left, IR_right, DR = make_transposon()
        pos = random.randint(0, len(seq))
        seq = seq[:pos] + transposon + seq[pos:]
        positions.append((pos, pos + len(transposon), IR_left, IR_right, DR))

    return seq, positions

def detect_transposons(sequence, min_IR=6, max_IR=10):
    found = []

    n = len(sequence)

    for i in range(n - min_IR):
        # try every possible IR length
        for L in range(min_IR, max_IR + 1):

            if i + L > n:
                continue
            
            IR_left = sequence[i:i+L]
            IR_right = reverse_complement(IR_left)

            for j in range(i + L, n - L):
                if sequence[j:j+L] == IR_right:

                    DR_left = sequence[i - 3:i] if i >= 3 else None
                    DR_right = sequence[j+L:j+L+3] if j+L+3 <= n else None

                    # If DRs match at least partially, accept
                    if DR_left is not None and DR_right is not None and DR_left[:2] == DR_right[:2]:
                        found.append((i-3, j+L+3))
    
    return found


base_length = random.randint(200, 400)
base_seq = random_dna(base_length)

full_seq, inserted_positions = insert_transposons(base_seq, num_TEs=random.randint(3,4))

print("\n=== INSERTED TRANSPOSONS (TRUE POSITIONS) ===")
for pos in inserted_positions:
    print(f"Start={pos[0]}, End={pos[1]}")

detected = detect_transposons(full_seq)

print("\n=== DETECTED TRANSPOSONS ===")
for d in detected:
    print(f"Detected start={d[0]}, end={d[1]}")

print(f"\nFinal DNA length: {len(full_seq)} bp")
