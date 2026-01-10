import matplotlib.pyplot as plt
import numpy as np

s1 = "ACCGTGAAGCCAATAC"
s2 = "AGCGTGCAGCCAATAC"

gap = -1
match = 1
mismatch = -1

def nw_align(a, b):
    n = len(a)
    m = len(b)
    score = []
    for i in range(n+1):
        row = []
        for j in range(m+1):
            row.append(0)
        score.append(row)

    for i in range(1, n+1):
        score[i][0] = score[i-1][0] + gap
    for j in range(1, m+1):
        score[0][j] = score[0][j-1] + gap

    for i in range(1, n+1):
        for j in range(1, m+1):
            if a[i-1] == b[j-1]:
                diag = score[i-1][j-1] + match
            else:
                diag = score[i-1][j-1] + mismatch
            up = score[i-1][j] + gap
            left = score[i][j-1] + gap
            v = diag
            if up > v: v = up
            if left > v: v = left
            score[i][j] = v

    i = n
    j = m
    al1 = ""
    al2 = ""
    path = [(i, j)]

    while i > 0 or j > 0:
        cur = score[i][j]
        moved = False

        if i > 0 and j > 0:
            d = match if a[i-1] == b[j-1] else mismatch
            if cur == score[i-1][j-1] + d:
                al1 = a[i-1] + al1
                al2 = b[j-1] + al2
                i -= 1; j -= 1
                path.append((i, j))
                continue

        if i > 0 and cur == score[i-1][j] + gap:
            al1 = a[i-1] + al1
            al2 = "-" + al2
            i -= 1
            path.append((i, j))
            continue

        if j > 0 and cur == score[i][j-1] + gap:
            al1 = "-" + al1
            al2 = b[j-1] + al2
            j -= 1
            path.append((i, j))
            continue

        break

    return al1, al2, score, path

al1, al2, score, path = nw_align(s1, s2)

print("Alignment:")
print(al1)
print("".join("|" if al1[i]==al2[i] else " " for i in range(len(al1))))
print(al2)

matches = sum(1 for i in range(len(al1)) if al1[i] == al2[i])
print("\nMatches:", matches)
print("Length :", len(al1))
print("Similarity:", round(matches/len(al1)*100, 2), "%")

# ---- Visualization -----

mat = np.array(score)
p = np.zeros_like(mat)
for (i,j) in path:
    p[i][j] = 1

fig, axes = plt.subplots(1, 2, figsize=(12,6))

ax = axes[0]
im = ax.imshow(mat, cmap="inferno")
ax.set_title("Score Matrix (Needleman–Wunsch)")
ax.set_xlabel("Seq2")
ax.set_ylabel("Seq1")
plt.colorbar(im, ax=ax)

ax2 = axes[1]
ax2.imshow(np.zeros_like(mat)+1, cmap="Wistia")  
for (i,j) in path:
    ax2.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, color="red"))
ax2.set_title("Traceback Path")
ax2.set_xlabel("Seq2")
ax2.set_ylabel("Seq1")

plt.tight_layout()
plt.show()
