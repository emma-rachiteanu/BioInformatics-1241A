# Laboratory 1 Problem 3
# Use AI to adapt your current algorithm in order to make an app that takes a FASTA file and read the sequence content from it and display the relative percentages
# for the symbols present in the alphabet of the sequence.
# Note: FASTA represents a file format that contains DNA, ARN or protons sequence. Thus, it contains the information for your input.

def read_fasta(file_path):
    """Reads a FASTA file and returns the concatenated sequence as a string."""
    sequence = ""
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue  # skip headers and empty lines
            sequence += line.upper()
    return sequence

def calculate_percentages(seq):
    letters = [i for i in seq if i.isalpha()]
    alph = sorted(list(set(letters)))
    tot = len(letters)
    perc = {i: letters.count(i)/tot*100 for i in alph}
    return alph, perc

def main():
    fasta_file = "sample.fasta"  # FASTA file path
    sequence = read_fasta(fasta_file)
    a, p = calculate_percentages(sequence)
    print("Alphabet of the sequence: " + ''.join(a))
    print("Percentages:")
    for i in a:
        print(f"{i}: {p[i]:.2f}%")

if __name__ == "__main__":
    main()
