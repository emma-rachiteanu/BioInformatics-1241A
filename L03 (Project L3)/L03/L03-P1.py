"""
Lab 3 problem 1

The melting temp (tm) is the temp at which one half of a particular DNA will disociate and become a single strand of DNA primerland and seq are of critical importance in 
designing the parameters of a successful amplification. The tm of a nucleic acid duplex increases both with its length and with increasing GC content.

tm = 4 (G+ C) + 2 (A + T)
tm = 81.5 + 16.6 (log10([Na+])) + .41*(%GC) - 600/length 

Implement an app that computes the tm of a DNA seq, by using one of these formula or both of them. 

Input = a string of DNA
Output = temp in Celcius
"""

import math

dna = input("DNA sequence (A,T,G,C): ")

A = dna.upper().count("A")
T = dna.upper().count("T")
G = dna.upper().count("G")
C = dna.upper().count("C")

tm1 = 4 * (G + C) + 2 * (A + T)

length = len(dna)
if length == 0:
    print("???")
else:
    GCpercent = ((G + C) / length) * 100
    Na_conc = 0.01
    tm2 = 81.5 + 16.6 * math.log10(Na_conc) + 0.41 * GCpercent - (600 / length)

    print("F1:", tm1, "C")
    print("F2:", round(tm2, 2), "C")