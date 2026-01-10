import random
import time
import heapq

# 1
with open("seq.fasta") as f:
    seq = "".join(line.strip().upper() for line in f if not line.startswith(">"))

seq = seq[:2000]
n = len(seq)

print("[1] Loaded DNA sequence")
print(f"   length used: {n} bp")
print("=" * 40)

# 2
num_reads = 2000
min_len, max_len = 100, 150
reads = []

for _ in range(num_reads):
    L = random.randint(min_len, max_len)
    s = random.randint(0, n - L)
    reads.append(seq[s:s+L])

print("[2–3] Random reads created")
print(f"   reads stored: {len(reads)}")
print("=" * 40)

# 3
def max_overlap(a: str, b: str, min_k: int) -> int:
    max_k = min(len(a), len(b))
    for k in range(max_k, min_k - 1, -1):
        if a.endswith(b[:k]):
            return k
    return 0

# 4
min_overlap = 30
max_merges = 800
use_subset = 600

active = reads[:use_subset]

alive = [True] * len(active)
version = [0] * len(active)

prefix_index = {}
def add_to_index(idx: int):
    if not alive[idx]:
        return
    p = active[idx][:min_overlap]
    prefix_index.setdefault(p, []).append(idx)

def rebuild_index():
    prefix_index.clear()
    for i in range(len(active)):
        if alive[i]:
            add_to_index(i)

rebuild_index()

heap = []

def push_candidates(i: int):
    if not alive[i]:
        return
    suf = active[i][-min_overlap:]
    candidates = prefix_index.get(suf, [])
    for j in candidates:
        if i == j or not alive[j]:
            continue
        k = max_overlap(active[i], active[j], min_overlap)
        if k > 0:
            heapq.heappush(heap, (-k, i, j, version[i], version[j]))

for i in range(len(active)):
    push_candidates(i)

print("[4] Greedy assembly (indexed + heap)")
print(f"   working reads: {sum(alive)} | min_overlap: {min_overlap}")
print("=" * 40)

t0 = time.time()
merge_count = 0

while heap and merge_count < max_merges:
    negk, i, j, vi, vj = heapq.heappop(heap)
    k = -negk

    if not (0 <= i < len(active) and 0 <= j < len(active)):
        continue
    if not (alive[i] and alive[j]):
        continue
    if version[i] != vi or version[j] != vj:
        continue

    merged = active[i] + active[j][k:]

    alive[i] = False
    alive[j] = False
    version[i] += 1
    version[j] += 1

    active.append(merged)
    alive.append(True)
    version.append(0)
    new_idx = len(active) - 1

    add_to_index(new_idx)
    push_candidates(new_idx)

    merge_count += 1
    if merge_count <= 25 or merge_count % 50 == 0:
        print(f"   merge {merge_count:3d}: overlap={k:3d} -> new contig length={len(merged)} | alive={sum(alive)}")

t1 = time.time()

longest = 0
for i in range(len(active)):
    if alive[i]:
        if len(active[i]) > longest:
            longest = len(active[i])

print("=" * 40)
print("Done")
print(f"merges performed: {merge_count}")
print(f"remaining fragments: {sum(alive)}")
print(f"longest contig: {longest} bp")
print(f"runtime: {t1 - t0:.2f} sec")


"""
The main problem with this approach is that the reconstruction is ambiguous.
Because real DNA contains many repeated or very similar regions, the short fragments can overlap in several equally valid ways. As a result, the algorithm cannot determine the correct order of the fragments and may either
assemble an incorrect sequence or stop early when no unambiguous overlap can be found.

"""