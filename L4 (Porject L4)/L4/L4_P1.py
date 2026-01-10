GENETIC_CODE = {
    # Phenylalanine / Leucine
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",

    # Serine / Proline
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",

    # Tyrosine / Histidine / Glutamine
    "UAU": "Tyr", "UAC": "Tyr",
    "CAU": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",

    # Cysteine / Tryptophan
    "UGU": "Cys", "UGC": "Cys", "UGG": "Trp",

    # Isoleucine / Methionine (start)
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile",
    "AUG": "Met",

    # Threonine / Asparagine / Lysine
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",

    # Valine / Alanine
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",

    # Aspartic acid / Glutamic acid
    "GAU": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",

    # Glycine
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",

    # Arginine / Serine
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser",
    "AGA": "Arg", "AGG": "Arg",

    # Stop codons
    "UAA": "Stop", "UAG": "Stop", "UGA": "Stop"
}


def translate_rna(rna):
    rna = rna.upper().replace(" ", "")
    protein = []

    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        amino = GENETIC_CODE.get(codon)

        if amino is None:
            protein.append("???")
            continue

        if amino == "Stop":
            break

        protein.append(amino)

    return "-".join(protein)


if __name__ == "__main__":
    seq = input("Enter RNA coding sequence:\n")
    print("Protein sequence:")
    print(translate_rna(seq))
