import requests
import math
import matplotlib.pyplot as plt

# influenza A/PR/8/34(H1N1)
FLU_SEGMENTS = [
    "NC_002023.1",
    "NC_002021.1",
    "NC_002022.1",
    "NC_002017.1",
    "NC_002019.1",
    "NC_002018.1",
    "NC_002016.1",
    "NC_002020.1",
]

COVID_ACCESSION = "NC_045512.2"

img_idx = 10


def save_and_show(fig):
    global img_idx
    fname = f"Screenshot_{img_idx}.jpg"
    fig.savefig(fname, dpi=200, bbox_inches="tight")
    print("Saved:", fname)
    plt.show()
    plt.close(fig)
    img_idx += 1


def download_fasta(accession):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        "?db=nuccore&id=" + accession + "&rettype=fasta&retmode=text"
    )
    r = requests.get(url)
    if r.status_code != 200:
        raise RuntimeError("download failed:" + accession)
    seq = []
    for line in r.text.splitlines():
        if not line or line.startswith(">"):
            continue
        seq.append(line.strip())
    return "".join(seq).upper()


def smith_waterman_score(a, b, match=2, mismatch=-1, gap=-2):
    n = len(a)
    m = len(b)
    H = []
    for i in range(n + 1):
        H.append([0] * (m + 1))
    best = 0
    for i in range(1, n + 1):
        ai = a[i - 1]
        Hi = H[i]
        Him1 = H[i - 1]
        for j in range(1, m + 1):
            s = match if ai == b[j - 1] else mismatch
            diag = Him1[j - 1] + s
            up = Him1[j] + gap
            left = Hi[j - 1] + gap
            v = diag
            if up > v:
                v = up
            if left > v:
                v = left
            if v < 0:
                v = 0
            Hi[j] = v
            if v > best:
                best = v
    return best


def make_windows(seq, size, step):
    wins = []
    i = 0
    while i < len(seq):
        chunk = seq[i:i+size]
        if len(chunk) < size * 0.25:
            break
        wins.append((i, chunk))
        i += step
    return wins


def build_similarity_map(seqA, seqB, win_size, step, match=2, mismatch=-1, gap=-2):
    winsA = make_windows(seqA, win_size, step)
    winsB = make_windows(seqB, win_size, step)
    mat = []

    for startA, wA in winsA:
        row = []
        for startB, wB in winsB:
            score = smith_waterman_score(wA, wB, match, mismatch, gap)
            smax = match * min(len(wA), len(wB))
            if smax <= 0:
                smax = 1
            row.append(score / float(smax))
        mat.append(row)

    return mat, winsA, winsB


def zscore_matrix(mat):
    vals = [v for row in mat for v in row]
    mean = sum(vals) / len(vals)
    var = sum((x - mean)**2 for x in vals) / len(vals)
    std = var**0.5 if var > 0 else 1.0

    out = []
    for row in mat:
        out.append([(v-mean)/std for v in row])
    return out, mean, std

def gc_series(seq, win, step):
    xs = []
    ys = []
    i = 0
    while i + win <= len(seq):
        frag = seq[i:i+win]
        g = frag.count("G")
        c = frag.count("C")
        gc = (g + c) / float(len(frag))
        xs.append(i)
        ys.append(gc)
        i += step
    return xs, ys


def entropy_series(seq, win, step):
    xs = []
    ys = []
    i = 0
    while i + win <= len(seq):
        frag = seq[i:i+win]
        count = {"A":0, "C":0, "G":0, "T":0}
        for ch in frag:
            if ch in count:
                count[ch]+=1
        H=0
        total = float(len(frag))
        for b in "ACGT":
            p = count[b]/total if total>0 else 0
            if p>0:
                H -= p*math.log(p,2)
        xs.append(i)
        ys.append(H)
        i += step
    return xs, ys


def base_composition(seq):
    d = {b:0 for b in "ACGT"}
    for ch in seq:
        if ch in d:
            d[ch]+=1
    t = float(sum(d.values()))
    return [d[b]/t for b in "ACGT"]


def all_kmers(k):
    if k == 1:
        return ["A","C","G","T"]
    prev = all_kmers(k-1)
    res=[]
    for p in prev:
        for b in "ACGT":
            res.append(p+b)
    return res


def kmer_counts(seq, k):
    d = {}
    for i in range(len(seq)-k+1):
        s = seq[i:i+k]
        if any(ch not in "ACGT" for ch in s):
            continue
        d[s] = d.get(s, 0) + 1
    return d


def kmer_vector(counts, k):
    km = all_kmers(k)
    vec=[]
    total = float(sum(counts.values())) if counts else 1.0
    for s in km:
        vec.append(counts.get(s,0)/total)
    return km, vec


def kmer_matrix_4mer(vec):
    # turn 256-mer vector into 16×16 FCGR
    rows=[]
    i=0
    for _r in range(16):
        row=[]
        for _c in range(16):
            row.append(vec[i])
            i+=1
        rows.append(row)
    return rows


def dotplot_points(seq1, seq2, k=10):
    pos_map={}
    for j in range(len(seq2)-k+1):
        s=seq2[j:j+k]
        if any(ch not in "ACGT" for ch in s):
            continue
        pos_map.setdefault(s,[]).append(j)
    xs=[]; ys=[]
    for i in range(len(seq1)-k+1):
        s=seq1[i:i+k]
        if s in pos_map:
            for j in pos_map[s]:
                xs.append(i)
                ys.append(j)
    return xs, ys


def cgr_points(seq, max_points=20000):
    corners={"A":(0,0),"C":(0,1),"G":(1,1),"T":(1,0)}
    x=0.5; y=0.5
    xs=[]; ys=[]
    step = max(1, len(seq)//max_points)
    idx=0
    for ch in seq:
        if ch not in corners:
            continue
        cx,cy = corners[ch]
        x=(x+cx)/2
        y=(y+cy)/2
        if idx % step == 0:
            xs.append(x)
            ys.append(y)
        idx+=1
    return xs,ys


def fcgr_matrix(seq, k=4):
    c = kmer_counts(seq,k)
    km, vec = kmer_vector(c,k)
    return kmer_matrix_4mer(vec)


def pearson(a, b):
    mA = sum(a)/len(a)
    mB = sum(b)/len(b)
    num = sum((a[i]-mA)*(b[i]-mB) for i in range(len(a)))
    dA = math.sqrt(sum((x-mA)**2 for x in a))
    dB = math.sqrt(sum((x-mB)**2 for x in b))
    return num/(dA*dB) if dA>0 and dB>0 else 0


def main():
    print("Downloading influenza...")
    flu_parts=[]
    for acc in FLU_SEGMENTS:
        print(" ",acc)
        flu_parts.append(download_fasta(acc))
    flu_seq="".join(flu_parts)

    print("Downloading SARS-CoV-2...")
    covid_seq = download_fasta(COVID_ACCESSION)

    print("Flu length:",len(flu_seq))
    print("Covid length:",len(covid_seq))

    sim1500, wins_flu_1500, wins_cov_1500 = build_similarity_map(
        flu_seq, covid_seq, 1500, 1200
    )
    z1500, _, _ = zscore_matrix(sim1500)

    # NLA
    flat_raw = [v for row in sim1500 for v in row]
    NLA = sum(flat_raw)/len(flat_raw) if flat_raw else 0
    print("\nNLA (Normalized Local Alignment):", NLA)

    # KPSS
    counts_f = kmer_counts(flu_seq,4)
    counts_c = kmer_counts(covid_seq,4)
    km,vec_f = kmer_vector(counts_f,4)
    _, vec_c = kmer_vector(counts_c,4)

    dot = sum(vec_f[i]*vec_c[i] for i in range(len(vec_f)))
    nf = math.sqrt(sum(v*v for v in vec_f))
    nc = math.sqrt(sum(v*v for v in vec_c))
    KPSS = dot/(nf*nc) if nf>0 and nc>0 else 0

    print("KPSS (k-mer profile similarity):", KPSS)

    # GCEPC
    xs_f_gc, ys_f_gc = gc_series(flu_seq,1000,500)
    xs_c_gc, ys_c_gc = gc_series(covid_seq,1000,500)
    xs_f_en, ys_f_en = entropy_series(flu_seq,1000,500)
    xs_c_en, ys_c_en = entropy_series(covid_seq,1000,500)

    n = min(len(ys_f_gc), len(ys_c_gc))
    m = min(len(ys_f_en), len(ys_c_en))

    rho_gc = pearson(ys_f_gc[:n], ys_c_gc[:n])
    rho_en = pearson(ys_f_en[:m], ys_c_en[:m])
    GCEPC = 0.5*(rho_gc + rho_en)

    print("GCEPC (GC+Entropy correlation):", GCEPC)


    fig, ax = plt.subplots(figsize=(12,6))
    img = ax.imshow(z1500, origin="lower", aspect="auto")
    ax.set_title("Local Alignment Similarity (z-score, 1500bp windows)")
    plt.colorbar(img, ax=ax)
    save_and_show(fig)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(xs_f_gc, ys_f_gc, label="Influenza")
    ax.plot(xs_c_gc, ys_c_gc, label="SARS-CoV-2")
    ax.legend()
    ax.set_title("GC content")
    save_and_show(fig)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(xs_f_en, ys_f_en, label="Influenza")
    ax.plot(xs_c_en, ys_c_en, label="SARS-CoV-2")
    ax.legend()
    ax.set_title("Entropy")
    save_and_show(fig)

    comp_f = base_composition(flu_seq)
    comp_c = base_composition(covid_seq)
    fig, ax = plt.subplots(figsize=(6,4))
    bases = ["A","C","G","T"]
    x = range(4)
    w = 0.35
    ax.bar([i-w/2 for i in x], comp_f, w, label="Influenza")
    ax.bar([i+w/2 for i in x], comp_c, w, label="SARS-CoV-2")
    ax.set_xticks(x)
    ax.set_xticklabels(bases)
    ax.legend()
    ax.set_title("Base Composition")
    save_and_show(fig)

    fig, (a1,a2) = plt.subplots(1,2,figsize=(10,4))
    fcgr_f = fcgr_matrix(flu_seq,4)
    fcgr_c = fcgr_matrix(covid_seq,4)
    a1.imshow(fcgr_f,origin="lower")
    a2.imshow(fcgr_c,origin="lower")
    a1.set_title("Influenza FCGR (4-mer)")
    a2.set_title("SARS-CoV-2 FCGR (4-mer)")
    save_and_show(fig)

    xs, ys = dotplot_points(flu_seq, covid_seq, k=10)
    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(xs, ys, s=1)
    ax.set_title("Dotplot 10-mer matches")
    save_and_show(fig)

    xf, yf = cgr_points(flu_seq)
    xc, yc = cgr_points(covid_seq)
    fig, (axf, axc) = plt.subplots(1,2,figsize=(8,4))
    axf.scatter(xf,yf,s=1)
    axc.scatter(xc,yc,s=1)
    axf.set_title("Influenza CGR")
    axc.set_title("SARS-CoV-2 CGR")
    save_and_show(fig)

    print("Finished. All screenshots saved.")


if __name__ == "__main__":
    main()
