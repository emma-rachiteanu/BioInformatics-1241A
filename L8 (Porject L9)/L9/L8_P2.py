"""
P2

Load from NCBI a total of 3 bacterial genomes. Use these genomes as input for your app. Also, modify your app in order to be able to handle the amount of information from these genomes.
Your application must detect transposable elements and the results form the output must show their position and their length.
In this case, the inverted repeats are unknown, unline the previous assignment in which they were known.
Note: the following cases must be taken into consideration
1. Transposon involving (one transposon being inside another transposon)
2. Transposon overalapping (one transposon beginning inside of the other transposon but continuing outside of it)
Note 2: The minimum size of the inverted repeats must be of four bases and the maximum size of the inverted repeats must be of six. So, anything between 4 and 6

In these cases we may use a sliding window. We write a small function which transforms this into an ivnerted repeat and with the result we make another sliding window in which we
want to find that inverted repeat somewhere.
"""

import requests
from bisect import bisect_left
from collections import defaultdict



NCBI_EFETCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
)

GENOMES = {
    "NC_000913.3": "Escherichia coli str. K-12 substr. MG1655",
    "NC_000964.3": "Bacillus subtilis subsp. subtilis str. 168",
    "NC_000962.3": "Mycobacterium tuberculosis H37Rv",
}


def fetch_fasta_sequence(accession: str) -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }

    print(f"Downloading {accession} from NCBI...")
    r = requests.get(NCBI_EFETCH_URL, params=params, timeout=60)
    r.raise_for_status()

    seq_lines = []
    for line in r.text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            # header line
            continue
        seq_lines.append(line.strip())
    sequence = "".join(seq_lines).upper()
    print(f"  -> downloaded {len(sequence):,} bp")
    return sequence



def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def sliding_windows(seq: str, k: int):
    n = len(seq)
    for i in range(n - k + 1):
        yield i, seq[i:i + k]


class TransposonDetector:

    def __init__(
        self,
        min_ir_len: int = 4,
        max_ir_len: int = 6,
        min_core: int = 50,
        max_core: int = 5000,
    ):
        self.min_ir_len = min_ir_len
        self.max_ir_len = max_ir_len
        self.min_core = min_core
        self.max_core = max_core

    def _index_kmers(self, seq: str, k: int):

        index = defaultdict(list)
        for pos, kmer in sliding_windows(seq, k):
            index[kmer].append(pos)
        return index

    def detect(self, seq: str):

        n = len(seq)
        results = []

        for k in range(self.min_ir_len, self.max_ir_len + 1):
            print(f"  Indexing {k}-mers...")
            index = self._index_kmers(seq, k)

            for kmer, positions in index.items():
                rc = reverse_complement(kmer)
                if rc not in index:
                    continue

                rc_positions = index[rc]

                for i in positions:
                    min_j = i + k + self.min_core
                    max_j = i + k + self.max_core

                    j_start_idx = bisect_left(rc_positions, min_j)

                    for j_idx in range(j_start_idx, len(rc_positions)):
                        j = rc_positions[j_idx]
                        if j > max_j:
                            break

                        core_len = j - (i + k)
                        te_start = i
                        te_end = j + k
                        te_len = te_end - te_start

                        results.append(
                            {
                                "start": te_start,
                                "end": te_end,
                                "length": te_len,
                                "ir_len": k,
                                "core_len": core_len,
                            }
                        )

        results.sort(key=lambda r: (r["start"], r["end"]))
        return results


def main():
    detector = TransposonDetector(
        min_ir_len=4,
        max_ir_len=6,
        min_core=50,
        max_core=5000,
    )

    for acc, name in GENOMES.items():
        print("\n" + "=" * 70)
        print(f"Genome: {name} ({acc})")
        print("=" * 70)

        seq = fetch_fasta_sequence(acc)

        print("Detecting candidate transposable elements...")
        tes = detector.detect(seq)

        print(f"\nTotal TEs found: {len(tes):,}")
        print("First 20 (1-based coordinates):")
        print("idx\tstart\tend\tlength\tIR_len\tcore_len")
        for idx, te in enumerate(tes[:20], start=1):
            start1 = te["start"] + 1
            end1 = te["end"]
            print(
                f"{idx}\t{start1}\t{end1}\t{te['length']}\t"
                f"{te['ir_len']}\t{te['core_len']}"
            )

        out_name = f"{acc}_te_candidates.tsv"
        with open(out_name, "w") as out:
            out.write("start_0based\tend_0based\tlength\tIR_len\tcore_len\n")
            for te in tes:
                out.write(
                    f"{te['start']}\t{te['end']}\t"
                    f"{te['length']}\t{te['ir_len']}\t{te['core_len']}\n"
                )
        print(f"\nAll candidates saved to: {out_name}")


if __name__ == "__main__":
    main()