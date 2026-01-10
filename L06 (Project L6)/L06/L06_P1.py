import random
import urllib.request
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# 1
ncbi_id = "NC_000913.3"
url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={ncbi_id}&db=nuccore&report=fasta"
fasta = urllib.request.urlopen(url).read().decode("utf-8")
seq_full = "".join(line.strip().upper() for line in fasta.splitlines() if not line.startswith(">"))

seq = seq_full[:random.randint(1000, 3000)]
n = len(seq)
print(f"Loaded {n} bp from {ncbi_id}")

# 2
def log_uniform_int(lo, hi):
    a = math.log10(lo)
    b = math.log10(hi)
    return int(round(10 ** random.uniform(a, b)))

samples = []
lengths = []
for _ in range(10):
    L = log_uniform_int(100, min(3000, n))
    s = random.randint(0, n - L)
    samples.append(seq[s:s+L])
    lengths.append(L)

# 3
print("Sample fragment lengths (bp):", lengths)

# 4
bp_min, bp_max = 100, 3000

def y_norm(bp):
    bp = max(bp_min, min(bp_max, bp))
    a = math.log10(bp_min)
    b = math.log10(bp_max)
    t = (math.log10(bp) - a) / (b - a)
    return 0.08 + (1.0 - t) * 0.84

ladder_bp = list(range(100, 1100, 100)) + [1500, 2000, 3000]

fig = plt.figure(figsize=(8.2, 4.8), dpi=300)
ax = plt.gca()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

gel_x0, gel_y0 = 0.12, 0.18
gel_w, gel_h = 0.84, 0.70

ax.add_patch(Rectangle((gel_x0, gel_y0), gel_w, gel_h, facecolor="#0b0b0b", edgecolor="black", linewidth=2.0))

lanes = 11  # M + 10 samples
lane_w = gel_w / lanes
lane_centers = [gel_x0 + lane_w * (i + 0.5) for i in range(lanes)]

# wells
well_y = gel_y0 + gel_h - 0.06 * gel_h
well_h = 0.045 * gel_h
for i in range(lanes):
    wx = gel_x0 + lane_w * i + lane_w * 0.18
    ww = lane_w * 0.64
    ax.add_patch(Rectangle((wx, well_y), ww, well_h, facecolor="#151515", edgecolor="none"))

for i in range(1, lanes):
    x = gel_x0 + lane_w * i
    ax.plot([x, x], [gel_y0, gel_y0 + gel_h], linewidth=0.5, alpha=0.12)

def to_plot_y(yN):
    return gel_y0 + gel_h * yN

def draw_smear_band(xc, y, w, h, strength):
    ax.add_patch(Rectangle((xc - w/2, y - (h*2.8)/2), w, h*2.8, facecolor="white", edgecolor="none", alpha=0.08 * strength))
    ax.add_patch(Rectangle((xc - w/2, y - h/2), w, h, facecolor="white", edgecolor="none", alpha=0.55 + 0.35 * strength))
    ax.add_patch(Rectangle((xc - w/2, y + h*0.15), w, h*0.18, facecolor="white", edgecolor="none", alpha=0.25))

band_h = gel_h * 0.010

# ladder lane (M)
for bp in ladder_bp:
    y = to_plot_y(y_norm(bp))
    w = lane_w * 0.55
    strength = 0.9 if bp in (500, 1500) else 0.6
    draw_smear_band(lane_centers[0], y, w, band_h, strength)

# sample lanes
for idx, L in enumerate(lengths, start=1):
    y = to_plot_y(y_norm(L))
    w = lane_w * 0.62
    strength = 0.65 + random.random() * 0.30
    draw_smear_band(lane_centers[idx], y, w, band_h * 1.1, strength)

# lane labels
ax.text(lane_centers[0], gel_y0 - 0.035, "M", ha="center", va="top", fontsize=11, color="black")
for i in range(1, lanes):
    ax.text(lane_centers[i], gel_y0 - 0.035, f"S{i}", ha="center", va="top", fontsize=10, color="black")

scale_x = gel_x0 - 0.03
ax.plot([scale_x, scale_x], [gel_y0, gel_y0 + gel_h], linewidth=2, color="black")
for bp in (3000, 1500, 500, 100):
    y = to_plot_y(y_norm(bp))
    ax.plot([scale_x - 0.01, scale_x + 0.01], [y, y], linewidth=2, color="black")
    ax.text(scale_x - 0.015, y, f"{bp}", ha="right", va="center", fontsize=10, color="black")

ax.text(gel_x0, gel_y0 + gel_h + 0.045, "Simulated gel electrophoresis", fontsize=20, ha="left", va="bottom", color="black")
ax.text(gel_x0, gel_y0 - 0.070, "Top (wells) → fragments start here", fontsize=12, ha="left", va="top", color="black")
ax.text(gel_x0, gel_y0 - 0.100, "↓ smaller fragments migrate further", fontsize=12, ha="left", va="top", color="black")

plt.savefig("Screenshot_1.jpg", dpi=300, bbox_inches="tight", pad_inches=0.12)
plt.close()
print("Saved Screenshot_1.jpg")
