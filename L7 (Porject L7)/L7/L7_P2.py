"""
P2
Download the influenza genomes. For each genome plot on a chart the met frequent repetitions
"""

import json
import random
from collections import Counter
from urllib import request, parse
from contextlib import closing
import matplotlib.pyplot as plt


def grab_sequence(min_size=1000, max_size=3000, pool=100, api=None, rng_seed=None):
    """Query NCBI for a random genomic record in a given size range."""
    if rng_seed is not None:
        random.seed(rng_seed)

    try:
        query = f"{min_size}:{max_size}[SLEN] AND biomol_genomic[PROP] NOT mitochondrial[Title] NOT chloroplast[Title]"
        params = {"db": "nuccore", "term": query, "retmode": "json", "retmax": str(pool)}
        if api:
            params["api_key"] = api

        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + parse.urlencode(params)
        with closing(request.urlopen(search_url, timeout=10)) as resp:
            listing = json.load(resp)

        id_options = listing.get("esearchresult", {}).get("idlist", [])
        if not id_options:
            raise RuntimeError("Search returned no usable IDs.")

        chosen_id = random.choice(id_options)

        fetch_params = {"db": "nuccore", "id": chosen_id, "rettype": "fasta", "retmode": "text"}
        if api:
            fetch_params["api_key"] = api

        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + parse.urlencode(fetch_params)
        with closing(request.urlopen(fetch_url, timeout=10)) as resp:
            fasta_block = resp.read().decode("utf-8", errors="ignore")

        title, seq = read_fasta(fasta_block)
        if not seq:
            raise RuntimeError("Downloaded sequence is empty.")
        return {"id": title, "seq": seq}

    except Exception as err:
        print(f"[NCBI ERROR] {err}")
        print("→ Using synthetic fallback sequence.")
        fake_len = random.randint(min_size, max_size)
        fake_seq = build_random_dna(fake_len)
        return {"id": "SYNTH_FAKE", "seq": fake_seq}


def read_fasta(raw):
    lines = [ln.strip() for ln in raw.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith(">"):
        return ("UNKNOWN", "")

    header = lines[0][1:]
    letters = "".join(lines[1:]).upper()
    letters = "".join(b for b in letters if b in "ACGTN")
    return (header, letters)


def build_random_dna(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def fetch_influenza_set(n=10):
    """Download several influenza genome segments from NCBI."""
    collected = []
    print(f"Retrieving {n} genomes...")
    for index in range(n):
        item = grab_sequence(min_size=1000, max_size=15000, pool=100)
        collected.append(item["seq"])
        print(f"  {index+1}/{n}: {item['id'][:50]}")
    return collected


def most_repeated_kmer(s, low=3, high=10):
    """Return the most frequent substring in the given DNA sequence."""
    countmap = Counter()

    for size in range(low, high + 1):
        for i in range(len(s) - size + 1):
            frag = s[i:i+size]
            if "N" not in frag:
                countmap[frag] += 1

    for frag, ct in countmap.most_common():
        if ct > 1:
            return frag, ct

    return None, 0


print("Downloading influenza sequences...")
genome_batch = fetch_influenza_set(10)

summary = []
print("\nScanning sequences for common repeats:")

for idx, genome in enumerate(genome_batch, start=1):
    rep, freq = most_repeated_kmer(genome)
    summary.append((idx, rep, freq))
    print(f"Seq {idx}: repeat '{rep}' → {freq} occurrences")


indexes = [item[0] for item in summary]
counts = [item[2] for item in summary]

plt.figure(figsize=(10, 6))
plt.plot(indexes, counts, marker='o', linewidth=2)
plt.scatter(indexes, counts, s=100)

plt.title("Most Common Repeat Frequency per Genome (Influenza)")
plt.xlabel("Genome Index")
plt.ylabel("Repeat Occurrence")
plt.grid(True, linestyle="--", alpha=0.4)
plt.xticks(indexes)

plt.tight_layout()
plt.show()
