import math
import re
import json
import tkinter as tk
from tkinter import scrolledtext
from collections import defaultdict, Counter

AUTHOR_A_TEXT = r"""
Somnoroase păsărele
Pe la cuiburi se adună,
Se ascund în rămurele —
Noapte bună!
Luna plină-n zarea lin
Varsă razele-n păduri,
Iar vântul șoapte-n codru spune.
"""

AUTHOR_B_TEXT = r"""
Leoaică tânără, iubirea
mi-a sărit în față.
Mă pândise-n încordare
mai demult.
Mi-a cerut o inimă
și o trecere prin foc,
iar eu nu aveam
decât cuvinte și tăcere.
"""

MY_TEXT = r"""
Somnoroase păsărele în codru se adună,
iubirea mi-a sărit pe frunze, fără noapte bună,
luna plină varsă raze în cuvinte și tăcere,
iar vântul spune-n șoaptă o trecere prin durere.
"""

WORD_RE = re.compile(r"[a-zăâîșțţ]+", flags=re.IGNORECASE)

def tokenize(text):
    text = text.lower().strip()
    tokens = WORD_RE.findall(text)
    out = []
    i = 0
    while i < len(tokens):
        if tokens[i]:
            out.append(tokens[i])
        i += 1
    return out

def bigrams(tokens):
    pairs = []
    i = 0
    while i < len(tokens) - 1:
        pairs.append((tokens[i], tokens[i + 1]))
        i += 1
    return pairs

def transition_counts(tokens):
    trans = defaultdict(Counter)
    pairs = bigrams(tokens)
    i = 0
    while i < len(pairs):
        a, b = pairs[i]
        trans[a][b] += 1
        i += 1
    return trans

def build_vocab(a, b):
    vocab = set()
    i = 0
    while i < len(a):
        vocab.add(a[i])
        i += 1
    i = 0
    while i < len(b):
        vocab.add(b[i])
        i += 1
    vocab.add("<UNK>")
    return sorted(vocab)

def smoothed_prob(counts, prev, nxt, V, alpha):
    out = counts.get(prev)
    if out is None:
        return alpha / (alpha * V)

    total = 0
    for k in out:
        total += out[k]

    numer = out.get(nxt, 0) + alpha
    denom = total + alpha * V

    if denom == 0:
        return alpha / (alpha * V)

    return numer / denom

def build_llr(tokens_a, tokens_b, alpha):
    counts_a = transition_counts(tokens_a)
    counts_b = transition_counts(tokens_b)
    vocab = build_vocab(tokens_a, tokens_b)
    vocab_set = set(vocab)
    V = len(vocab)

    def llr(prev, nxt):
        p = prev if prev in vocab_set else "<UNK>"
        q = nxt if nxt in vocab_set else "<UNK>"
        pa = smoothed_prob(counts_a, p, q, V, alpha)
        pb = smoothed_prob(counts_b, p, q, V, alpha)
        pa = max(pa, 1e-12)
        pb = max(pb, 1e-12)
        return math.log(pa) - math.log(pb)

    meta = {
        "alpha": alpha,
        "vocab_size": V
    }
    return llr, meta

def avg_llr(tokens, llr_func):
    pairs = bigrams(tokens)
    if not pairs:
        return 0.0

    s = 0.0
    i = 0
    while i < len(pairs):
        a, b = pairs[i]
        s += llr_func(a, b)
        i += 1
    return s / len(pairs)

def windows(tokens, w, step):
    n = len(tokens)
    if n <= w:
        yield 0, tokens
        return

    i = 0
    while i <= n - w:
        yield i, tokens[i:i + w]
        i += step

def label(score, eps):
    if score > eps:
        return "A"
    if score < -eps:
        return "B"
    return "N"

def make_report():
    tokens_a = tokenize(AUTHOR_A_TEXT)
    tokens_b = tokenize(AUTHOR_B_TEXT)
    tokens_m = tokenize(MY_TEXT)

    llr_func, meta = build_llr(tokens_a, tokens_b, alpha=0.5)

    n = len(tokens_m)
    if n < 20:
        W, STEP = 5, 1
    else:
        W, STEP = 7, 2

    EPS = 0.25

    out = []
    out.append("STYLOMETRY ANALYSIS (WORD TRANSITIONS)")
    out.append("=" * 90)
    out.append(f"Author A tokens : {len(tokens_a)}")
    out.append(f"Author B tokens : {len(tokens_b)}")
    out.append(f"Analysed tokens : {len(tokens_m)}")
    out.append("")
    out.append(f"alpha = {meta['alpha']}, vocab size = {meta['vocab_size']}")
    out.append(f"window = {W}, step = {STEP}, eps = {EPS}")
    out.append("")
    out.append("WINDOW RESULTS")
    out.append("-" * 90)
    out.append(f"{'start':>6} {'end':>6} {'avgLLR':>10} {'label':>6}  preview")

    votes = [0.0] * len(tokens_m)
    counts = [0] * len(tokens_m)

    for start, chunk in windows(tokens_m, W, STEP):
        sc = avg_llr(chunk, llr_func)
        lab = label(sc, EPS)
        end = start + len(chunk) - 1
        preview = " ".join(chunk[:6])
        out.append(f"{start:6d} {end:6d} {sc:10.4f} {lab:>6}  {preview}")

        j = start
        while j <= end and j < len(tokens_m):
            votes[j] += sc
            counts[j] += 1
            j += 1

    out.append("")
    out.append("TOKEN TIMELINE")
    out.append("-" * 90)
    out.append(f"{'idx':>4} {'token':<16} {'vote':>10} {'lab':>4}")

    labels = []
    i = 0
    while i < len(tokens_m):
        v = votes[i] / counts[i] if counts[i] > 0 else 0.0
        lab = label(v, EPS)
        labels.append(lab)
        out.append(f"{i:4d} {tokens_m[i]:<16} {v:10.4f} {lab:>4}")
        i += 1

    out.append("")
    out.append("SUMMARY")
    out.append("-" * 90)
    out.append(f"A-like: {labels.count('A') * 100.0 / len(labels):.2f}%")
    out.append(f"B-like: {labels.count('B') * 100.0 / len(labels):.2f}%")
    out.append(f"Neutral: {labels.count('N') * 100.0 / len(labels):.2f}%")

    return "\n".join(out)

def show_popup(text):
    root = tk.Tk()
    root.title("Exercise 2 — Stylometry Report")
    root.geometry("1100x750")

    box = scrolledtext.ScrolledText(root, wrap=tk.WORD, font=("Consolas", 10))
    box.pack(fill=tk.BOTH, expand=True)
    box.insert(tk.END, text)
    box.config(state=tk.DISABLED)

    root.mainloop()

if __name__ == "__main__":
    report = make_report()
    show_popup(report)
