import random
import urllib.request
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# 1
NCBI_ID = "NC_000913.3"
SEQ_LEN = random.randint(1000, 3000)
CIRCULAR_DNA = True

url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={NCBI_ID}&db=nuccore&report=fasta"
fasta = urllib.request.urlopen(url).read().decode("utf-8")
seq_full = "".join(line.strip().upper() for line in fasta.splitlines() if not line.startswith(">"))

# random slice
max_start = max(0, len(seq_full) - SEQ_LEN)
start = random.randint(0, max_start) if max_start > 0 else 0
dna = seq_full[start:start + SEQ_LEN]
n = len(dna)

print(f"NCBI ID: {NCBI_ID}")
print(f"DNA length used: {n} bp (slice start index: {start})")

# 2
enzymes = {
    "EcoRI":   ("GAATTC", 1),  # G^AATTC
    "BamHI":   ("GGATCC", 1),  # G^GATCC
    "HindIII": ("AAGCTT", 1),  # A^AGCTT
    "TaqI":    ("TCGA",   1),  # T^CGA
    "HaeIII":  ("GGCC",   2),  # GG^CC
}

# 3
def find_cut_positions(seq: str, recog: str, cut_offset: int):
    cuts = []
    i = 0
    while True:
        j = seq.find(recog, i)
        if j == -1:
            break
        p = j + cut_offset
        if 0 < p < len(seq):
            cuts.append(p)
        i = j + 1
    return sorted(set(cuts))

def fragments_from_cuts(seq_len: int, cuts, circular: bool):
    cuts = sorted(set(c for c in cuts if 0 < c < seq_len))

    if not cuts:
        return [seq_len]

    if circular:
        pts = cuts + [cuts[0] + seq_len]
        frags = []
        for a, b in zip(pts, pts[1:]):
            frags.append(b - a)
        return frags
    else:
        pts = [0] + cuts + [seq_len]
        frags = []
        for a, b in zip(pts, pts[1:]):
            frags.append(b - a)
        return frags

digest = {}
print("\n--- Digest output ---")
for name, (site, offset) in enzymes.items():
    cuts = find_cut_positions(dna, site, offset)
    frags = fragments_from_cuts(n, cuts, CIRCULAR_DNA)
    frags_sorted = sorted(frags, reverse=True)
    digest[name] = {"site": site, "cuts": cuts, "fragments": frags_sorted}

    print(f"\n{name} (site={site})")
    print(f"Number of cleavages: {len(cuts)}")
    if cuts:
        print("Cut positions (1-based):", [c + 1 for c in cuts])
    else:
        print("Cut positions: none")
    print("Fragment lengths:", frags_sorted)

# 4
bp_min, bp_max = 50, max(3000, n)

def y_from_bp(bp):
    bp = max(bp_min, min(bp_max, bp))
    a = math.log10(bp_min)
    b = math.log10(bp_max)
    t = (math.log10(bp) - a) / (b - a)
    return 0.08 + (1.0 - t) * 0.84

def draw_band(ax, xc, y, w, h, strength):
    ax.add_patch(Rectangle((xc - w/2, y - (h*2.6)/2), w, h*2.6,
                           facecolor="white", edgecolor="none", alpha=0.05 * strength))
    ax.add_patch(Rectangle((xc - w/2, y - h/2), w, h,
                           facecolor="white", edgecolor="none", alpha=0.55 + 0.35 * strength))

ladder_bp = (
    list(range(100, 1001, 100)) +
    [1500, 2000, 2500, 3000]
)
ladder_bp = [bp for bp in ladder_bp if bp_min <= bp <= bp_max]

lane_labels = ["M"] + list(enzymes.keys())
lane_count = len(lane_labels)

fig, ax = plt.subplots(figsize=(10.8, 4.8), dpi=300)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
fig.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.24)

gel_x0, gel_y0 = 0.10, 0.18
gel_w, gel_h = 0.88, 0.64

ax.add_patch(Rectangle((gel_x0, gel_y0), gel_w, gel_h,
                       facecolor="#0a0a0a", edgecolor="black", linewidth=2.0))

lane_w = gel_w / lane_count
lane_centers = [gel_x0 + lane_w * (i + 0.5) for i in range(lane_count)]

well_h = 0.06 * gel_h
well_y = gel_y0 + gel_h - well_h
for i in range(lane_count):
    wx = gel_x0 + lane_w * i + lane_w * 0.18
    ww = lane_w * 0.64
    ax.add_patch(Rectangle((wx, well_y), ww, well_h, facecolor="#151515", edgecolor="none"))

for i in range(1, lane_count):
    x = gel_x0 + lane_w * i
    ax.plot([x, x], [gel_y0, gel_y0 + gel_h], linewidth=0.6, alpha=0.10)

band_h = gel_h * 0.012
band_w_ladder = lane_w * 0.55
band_w_sample = lane_w * 0.62

for bp in ladder_bp:
    y = gel_y0 + gel_h * y_from_bp(bp)
    strength = 1.0 if bp in (500, 1500) else 0.75
    draw_band(ax, lane_centers[0], y, band_w_ladder, band_h, strength)

for lane_i, enz in enumerate(enzymes.keys(), start=1):
    frags = digest[enz]["fragments"]

    max_bands = 18
    frags = frags[:max_bands]

    for L in frags:
        y = gel_y0 + gel_h * y_from_bp(L)
        strength = 0.75 + 0.25 * random.random()
        draw_band(ax, lane_centers[lane_i], y, band_w_sample, band_h, strength)

for i, label in enumerate(lane_labels):
    ax.text(lane_centers[i], gel_y0 - 0.055, label, ha="center", va="top")

scale_x = gel_x0 - 0.035
ax.plot([scale_x, scale_x], [gel_y0, gel_y0 + gel_h], linewidth=2.2, color="black")
for bp in (3000, 1500, 500, 100):
    y = gel_y0 + gel_h * y_from_bp(bp)
    ax.plot([scale_x - 0.010, scale_x + 0.010], [y, y], linewidth=2.2, color="black")
    ax.text(scale_x - 0.015, y, f"{bp}", ha="right", va="center")

ax.text(gel_x0, gel_y0 - 0.095, "Top (wells)", ha="left", va="top")
ax.text(gel_x0, gel_y0 - 0.15, "↓ smaller fragments run further", ha="left", va="top")

ax.set_title("Simulated gel electrophoresis (restriction digests)", fontsize=22, pad=12)

plt.savefig("Screenshot_2.jpg", dpi=300, bbox_inches="tight", pad_inches=0.10)
plt.close()
print("\nSaved: Screenshot_2.jpg")
