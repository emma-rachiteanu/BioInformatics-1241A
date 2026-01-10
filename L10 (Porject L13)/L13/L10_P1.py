import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WINDOW_SIZE = 30
MAX_PROMOTERS = 5

S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"

def load_fasta(path, limit):
    seqs = {}
    name=None
    buff=[]
    count=0
    with open(path,"r") as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if count>=limit: break
                if name is not None:
                    seqs[name]="".join(buff)
                    count+=1
                    if count>=limit: break
                name=line[1:].strip()
                buff=[]
            else:
                if count<limit:
                    buff.append(line.upper())
    if name is not None and count<limit:
        seqs[name]="".join(buff)
    return seqs

def sliding_windows(seq,w):
    return [seq[i:i+w] for i in range(len(seq)-w+1)]

def cg_percent(seq):
    if len(seq)==0: return 0
    return (seq.count("C")+seq.count("G"))*100/len(seq)

def raw_ic(seq):
    if len(seq)<2: return 0
    d={}
    for x in seq: d[x]=d.get(x,0)+1
    n=len(seq)
    num=sum([v*(v-1) for v in d.values()])
    return num*100/(n*(n-1))

_raw_S = raw_ic(S)
IC_SCALE = 27.53/_raw_S if _raw_S!=0 else 1

def kappa_ic(seq):
    return raw_ic(seq)*IC_SCALE

def pattern(seq):
    wins = sliding_windows(seq,WINDOW_SIZE)
    xs=[]; ys=[]
    for w in wins:
        xs.append(cg_percent(w))
        ys.append(kappa_ic(w))
    if xs:
        cx=sum(xs)/len(xs); cy=sum(ys)/len(ys)
    else:
        cx=0; cy=0
    return xs,ys,(cx,cy)

def main():
    print("Global C+G% for S =", cg_percent(S))
    print("Global Kappa IC for S =", kappa_ic(S))

    xs,ys,c = pattern(S)

    plt.figure()
    plt.scatter(xs,ys,s=15)
    plt.scatter([c[0]],[c[1]],marker='x')
    plt.title("Pattern of S")
    plt.xlabel("C+G%"); plt.ylabel("Kappa IC")
    plt.savefig("pattern_S.png")
    plt.close()

    promoters = load_fasta("Promotori lista completa.fasta", MAX_PROMOTERS)

    centers = {}
    for name,seq in promoters.items():
        px,py,pc = pattern(seq)
        centers[name] = pc

        plt.figure()
        plt.scatter(px,py,s=10)
        plt.scatter(pc[0],pc[1],marker='x')
        plt.title(name)
        plt.xlabel("C+G%"); plt.ylabel("Kappa IC")
        plt.savefig(f"{name}_pattern.png")
        plt.close()

    plt.figure()
    for n,(x,y) in centers.items():
        plt.scatter(x,y)
        plt.text(x,y,n)
    plt.title("Centers of Patterns (Promoters)")
    plt.xlabel("C+G% center"); plt.ylabel("Kappa IC center")
    plt.savefig("promoters_centers.png")
    plt.close()

if __name__=="__main__":
    main()
