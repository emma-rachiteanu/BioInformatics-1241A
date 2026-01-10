# Positive log-likelihood scores indicate similarity to the exon–intron boundary model.
# Clear local maxima above background suggest the presence of a splice site signal.

import math

motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT"
]

alphabet = ['A', 'C', 'G', 'T']
num_sequences = len(motifs)
motif_length = len(motifs[0])

# 2
count_matrix = {}
for base in alphabet:
    count_matrix[base] = [0] * motif_length

for seq in motifs:
    for i in range(motif_length):
        count_matrix[seq[i]][i] += 1

# 3
frequency_matrix = {}
for base in alphabet:
    frequency_matrix[base] = []
    for i in range(motif_length):
        frequency_matrix[base].append(count_matrix[base][i] / num_sequences)

# 4
null_prob = 0.25
log_likelihood_matrix = {}

for base in alphabet:
    log_likelihood_matrix[base] = []
    for i in range(motif_length):
        p = frequency_matrix[base][i]
        if p == 0:
            log_likelihood_matrix[base].append(0)
        else:
            log_likelihood_matrix[base].append(math.log(p / null_prob))

# 5
S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

def score_window(window, ll_matrix):
    score = 0.0
    for i in range(len(window)):
        score += ll_matrix[window[i]][i]
    return score

window_scores = []
for i in range(len(S) - motif_length + 1):
    window = S[i:i + motif_length]
    score = score_window(window, log_likelihood_matrix)
    window_scores.append((i, window, score))

# Interpretation:
# Windows with the highest scores represent candidate exon–intron borders.
# Scores significantly higher than surrounding windows indicate a true signal.

best_window = max(window_scores, key=lambda x: x[2])

print("Sliding window scores:")
for pos, window, score in window_scores:
    print(f"Position {pos}: {window} -> {score:.2f}")

print("\nBest scoring window:")
print(f"Position {best_window[0]}: {best_window[1]} -> {best_window[2]:.2f}")
