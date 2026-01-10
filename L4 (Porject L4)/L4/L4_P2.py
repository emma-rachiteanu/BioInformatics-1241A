import matplotlib.pyplot as plt
from collections import Counter

GENETIC_CODE = {
    "UUU":"Phe","UUC":"Phe","UUA":"Leu","UUG":"Leu",
    "UCU":"Ser","UCC":"Ser","UCA":"Ser","UCG":"Ser",
    "UAU":"Tyr","UAC":"Tyr","UAA":"Stop","UAG":"Stop",
    "UGU":"Cys","UGC":"Cys","UGA":"Stop","UGG":"Trp",
    "CUU":"Leu","CUC":"Leu","CUA":"Leu","CUG":"Leu",
    "CCU":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "CAU":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
    "CGU":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "AUU":"Ile","AUC":"Ile","AUA":"Ile","AUG":"Met",
    "ACU":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "AAU":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
    "AGU":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
    "GUU":"Val","GUC":"Val","GUA":"Val","GUG":"Val",
    "GCU":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "GAU":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "GGU":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly"
}


def read_fasta(path):
    seq = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    return seq.replace("T", "U")


def count_codons(seq):
    codons = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3 and "N" not in codon:
            codons.append(codon)
    return Counter(codons)


def save_barplot(counter, title, screenshot_id):
    top = counter.most_common(10)
    labels = [c for c, _ in top]
    values = [v for _, v in top]

    plt.figure()
    plt.bar(labels, values)
    plt.title(title)
    plt.xlabel("Codon")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(f"Screenshot_{screenshot_id}.jpg", dpi=300)
    plt.close()


if __name__ == "__main__":
    covid_seq = read_fasta("covid.fasta")
    flu_seq   = read_fasta("influenza.fasta")

    covid_counts = count_codons(covid_seq)
    flu_counts   = count_codons(flu_seq)

    save_barplot(covid_counts, "Top 10 Codons – COVID", 1)
    save_barplot(flu_counts, "Top 10 Codons – Influenza", 2)

    print("Top 3 amino acids (COVID):")
    aa_covid = Counter(GENETIC_CODE[c] for c in covid_counts if c in GENETIC_CODE and GENETIC_CODE[c] != "Stop")
    print(aa_covid.most_common(3))

    print("\nTop 3 amino acids (Influenza):")
    aa_flu = Counter(GENETIC_CODE[c] for c in flu_counts if c in GENETIC_CODE and GENETIC_CODE[c] != "Stop")
    print(aa_flu.most_common(3))
