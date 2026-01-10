"""
Download from Moodle the FASTA file containing promotor sequences and use it as an input for an updated application which is able to
generate and save the objective digital stain (ODS) inside the folder. The promotor file can be found inside the PromKappa package on
Moodle or GitHub inside folder "bin"
"""
import matplotlib.pyplot as plt
import os

WINDOW_SIZE = 30

def load_fasta(path):
    seqs = {}
    name = None
    buff = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buff)
                    buff = []
                name = line[1:].strip()
            else:
                buff.append(line.upper())

    if name is not None:
        seqs[name] = "".join(buff)

    return seqs

def sliding_windows(seq, w):
    return [seq[i:i+w] for i in range(len(seq)-w+1)]

def cg_percent(seq):
    return (seq.count("C") + seq.count("G")) * 100 / len(seq)

def raw_ic(seq):
    n=len(seq)
    if n<2: return 0
    d={}
    for c in seq: d[c]=d.get(c,0)+1
    num=sum([v*(v-1) for v in d.values()])
    return (num*100)/(n*(n-1))

TEST_SEQ="CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
SCALE = 27.53/raw_ic(TEST_SEQ)

def kappa_ic(seq):
    return raw_ic(seq) * SCALE

def pattern(seq):
    wins = sliding_windows(seq, WINDOW_SIZE)
    xs=[]; ys=[]
    for w in wins:
        xs.append(cg_percent(w))
        ys.append(kappa_ic(w))
    cx = sum(xs)/len(xs) if xs else 0
    cy = sum(ys)/len(ys) if ys else 0
    return xs, ys, (cx, cy)

def save_pattern(name, xs, ys, center, outdir):
    plt.figure(figsize=(6,6))
    plt.scatter(xs, ys, s=15)
    plt.scatter(center[0], center[1], c="red", marker="x")
    plt.title(name)
    plt.xlabel("C+G %")
    plt.ylabel("Kappa IC")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{name}_pattern.png"))
    plt.close()

def main():
    fasta_path = "Promotori lista completa.fasta"

    if not os.path.exists(fasta_path):
        print("FASTA file missing!")
        return

    seqs = load_fasta(fasta_path)
    print("Loaded promoters:", list(seqs.keys()))

    if not os.path.exists("chart"):
        os.mkdir("chart")

    centers = {}

    for name, seq in seqs.items():
        xs, ys, center = pattern(seq)
        centers[name] = center
        save_pattern(name, xs, ys, center, "chart")

    plt.figure(figsize=(7,7))
    for name, (x,y) in centers.items():
        plt.scatter(x, y)
        plt.text(x, y, name)
    plt.title("Centers of Weight (Promoters)")
    plt.xlabel("C+G % (center)")
    plt.ylabel("Kappa IC (center)")
    plt.tight_layout()

    plt.show()

if __name__=="__main__":
    main()
