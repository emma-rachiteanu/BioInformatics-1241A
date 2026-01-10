"""
Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).
Implement a software application that detects repetition between 3b and 6b in this DNA sequence
NOTE: Repetitive sequences refer to patterns that repeat N times
Minimum number of repetitions is 2
"""

import random, json
from urllib import request, parse
from contextlib import closing

# fasta

def parse_fasta(fasta_text):
    lines = [l.strip() for l in fasta_text.splitlines() if l.strip()]
    if not lines or not lines[0].startswith('>'):
        return ('UNKNOWN', '')
    header = lines[0][1:]
    seq = ''.join(lines[1:]).upper()
    seq = ''.join(c for c in seq if c in 'ACGTN')
    return (header, seq)


def random_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))


def fetch_ncbi_sequence(min_len=1000, max_len=3000, retmax=100, api_key=None, seed=None):
    if seed is not None:
        random.seed(seed)
    try:
        term = f"{min_len}:{max_len}[SLEN] AND biomol_genomic[PROP] NOT mitochondrial[Title] NOT chloroplast[Title]"
        params = {
            'db': 'nuccore',
            'term': term,
            'retmode': 'json',
            'retmax': str(retmax)
        }
        if api_key:
            params['api_key'] = api_key
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + parse.urlencode(params)
        with closing(request.urlopen(url, timeout=10)) as resp:
            data = json.load(resp)
        ids = data.get('esearchresult', {}).get('idlist', [])
        if not ids:
            raise RuntimeError('no ids :(')
        pick = random.choice(ids)

        params = {
            'db': 'nuccore',
            'id': pick,
            'rettype': 'fasta',
            'retmode': 'text'
        }
        if api_key:
            params['api_key'] = api_key
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + parse.urlencode(params)
        with closing(request.urlopen(url, timeout=10)) as resp:
            text = resp.read().decode('utf-8','ignore')

        head, seq = parse_fasta(text)
        if not seq or not (min_len <= len(seq) <= max_len):
            raise RuntimeError('bad len')
        return {'accession': head, 'sequence': seq}

    except Exception:
        L = random.randint(min_len, max_len)
        return {'accession': 'SYNTHETIC_RANDOM_SEQ', 'sequence': random_dna(L)}


# repeat detector

def find_repeats(seq, a=3, b=6, reps=2):
    out = []
    ln = len(seq)
    for size in range(a, b+1):
        for i in range(ln - size):
            pt = seq[i:i+size]
            c = 1
            j = i + size
            while seq[j:j+size] == pt:
                c += 1
                j += size
            if c >= reps:
                out.append((pt, i, c))
    return out


# run the thing
if __name__ == '__main__':
    data = fetch_ncbi_sequence()
    seq = data['sequence']
    print('id:', data['accession'])
    print('len:', len(seq))
    r = find_repeats(seq)
    for p, st, c in r:
        print("pattern:", p, "| length:", len(p), "| repeats:", c, "| position:", st)