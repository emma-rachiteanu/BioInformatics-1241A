#   P2
#   A DNA sequence is given: S = "ACGGGCATATGCGC". Make an application which is able to show the percentage of the components from the alphabet of the sequence S.
#   In other words, the input of the sequence S and the output is the alphabet of the sequence and the percentage of each letter in the alphabet found in the sequence S.

def calculate_percentages(seq):
    letters=[]
    for i in seq:
        if i.isalpha():
            letters+=[i.upper()]
    alph=[]
    for i in letters:
        if i not in alph:
            alph+=[i]
    alph.sort()
    perc={}
    tot=len(letters)
    for i in alph:
        c=0
        for j in letters:
            if j==i:
                c+=1
        perc[i]=c/tot*100
    return alph,perc

def main():
    s="ACGGGCATATGCGC"
    a,p=calculate_percentages(s)
    print("Alphabet of the sequence: "+''.join(a))
    print("Percentages:")
    for i in a:
        print(i,":",p[i],"%")

main()


