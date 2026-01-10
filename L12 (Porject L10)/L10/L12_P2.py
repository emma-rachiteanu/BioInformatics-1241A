# Peaks in the score profile indicate likely functional motif locations.
# X-axis (L) = genome position, Y-axis (P) = log-likelihood score.

import math
import matplotlib.pyplot as plt
import urllib.request

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
null_prob = 0.25

count_matrix = {b: [0]*motif_length for b in alphabet}
for seq in motifs:
    for i in range(motif_length):
        count_matrix[seq[i]][i] += 1

freq_matrix = {b: [] for b in alphabet}
for b in alphabet:
    for i in range(motif_length):
        freq_matrix[b].append(count_matrix[b][i] / num_sequences)

log_ll_matrix = {b: [] for b in alphabet}
for b in alphabet:
    for i in range(motif_length):
        p = freq_matrix[b][i]
        if p == 0:
            log_ll_matrix[b].append(0)
        else:
            log_ll_matrix[b].append(math.log(p / null_prob))

influenza_ids = [
    "NC_002023", "NC_007373", "NC_007366", "NC_007371", "NC_007367",
    "NC_007372", "NC_007370", "NC_007368", "NC_007374", "NC_007369"
]

def download_fasta(genome_id):
    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={genome_id}&db=nuccore&report=fasta"
    response = urllib.request.urlopen(url)
    fasta = response.read().decode("utf-8")
    lines = fasta.splitlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    return sequence

genomes = []
for gid in influenza_ids:
    seq = download_fasta(gid)
    genomes.append((gid, seq))

def score_window(window, ll_matrix):
    score = 0.0
    for i in range(len(window)):
        base = window[i]
        if base in ll_matrix:
            score += ll_matrix[base][i]
    return score

screenshot_index = 2

for genome_id, sequence in genomes:
    positions = []
    scores = []

    for i in range(len(sequence) - motif_length + 1):
        window = sequence[i:i + motif_length]
        score = score_window(window, log_ll_matrix)
        positions.append(i)
        scores.append(score)

    plt.figure()
    plt.plot(positions, scores)
    plt.xlabel("L (genomic location)")
    plt.ylabel("P (log-likelihood score)")
    plt.title(f"Motif signal profile – {genome_id}")

    plt.savefig(f"Screenshot_{screenshot_index}.jpg", dpi=300, bbox_inches="tight")
    screenshot_index += 1

    plt.show()
