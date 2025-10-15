"""
Lab 3 problem 2

Design an application that uses the sliding window method in order to read the melting temperature tm over the sequence s. Use a sliding window of 8 positions and choose
a FASTA file as input. Plot the results.
"""
import math
import matplotlib.pyplot as plt

def tm1(seq):
    seq = seq.upper()
    A = seq.count("A")
    T = seq.count("T")
    G = seq.count("G")
    C = seq.count("C")
    return 4*(G+C) + 2*(A+T)

def tm2(seq):
    seq = seq.upper()
    length = len(seq)
    G = seq.count("G")
    C = seq.count("C")
    GCpercent = ((G + C) / length) * 100
    Na_conc = 0.01
    val = 81.5 + 16.6 * math.log10(Na_conc) + 0.41 * GCpercent - (600 / length)
    return val

fname = "sample.fasta"

with open(fname) as f:
    stuff = f.read()

stuff = stuff.replace("\n", "")
if ">" in stuff:
    seq = stuff.split(">")[1]
    seq = seq[seq.find("A"):]
else:
    seq = stuff

seq = seq.strip().upper()

if len(seq) < 8:
    print("seq too short")
else:
    print("sliding window of 8 bases:\n")
    tm1_vals = []
    tm2_vals = []
    x = []
    for i in range(len(seq) - 7):
        window = seq[i:i+8]
        t1 = tm1(window)
        t2 = tm2(window)
        tm1_vals.append(t1)
        tm2_vals.append(t2)
        x.append(i)
        print(i, "-", i+8, window, "=> tm1:", round(t1,2), "C ; tm2:", round(t2,2), "C")

    plt.plot(x, tm1_vals, label="TM formula 1", marker="o")
    plt.plot(x, tm2_vals, label="TM formula 2", marker="o")
    plt.xlabel("Window start position")
    plt.ylabel("Melting Temp (°C)")
    plt.title("Sliding window TM over DNA sequence (8 bases)")
    plt.legend()
    plt.show()

print("done")
