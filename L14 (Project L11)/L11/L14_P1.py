import math

ALPHABET = ["A", "C", "G", "T"]
VALID = set(ALPHABET)

def count_transitions(seq):
    seq = seq.strip().upper()

    counts = {}
    totals = {}
    i = 0
    while i < len(ALPHABET):
        a = ALPHABET[i]
        totals[a] = 0
        counts[a] = {}
        j = 0
        while j < len(ALPHABET):
            b = ALPHABET[j]
            counts[a][b] = 0
            j += 1
        i += 1

    k = 0
    while k < len(seq) - 1:
        a = seq[k]
        b = seq[k + 1]
        if a in VALID and b in VALID:
            counts[a][b] += 1
            totals[a] += 1
        k += 1

    return counts, totals


def transition_probs(counts, totals):
    Tr = {}
    i = 0
    while i < len(ALPHABET):
        a = ALPHABET[i]
        Tr[a] = {}
        denom = totals[a]
        j = 0
        while j < len(ALPHABET):
            b = ALPHABET[j]
            if denom == 0:
                Tr[a][b] = 0.0
            else:
                Tr[a][b] = counts[a][b] / denom
            j += 1
        i += 1
    return Tr


def log_base(x, base):
    return math.log(x) / math.log(base)


def log_likelihood_matrix(Tr_plus, Tr_minus):
    beta = {}
    i = 0
    while i < len(ALPHABET):
        a = ALPHABET[i]
        beta[a] = {}
        j = 0
        while j < len(ALPHABET):
            b = ALPHABET[j]
            p = Tr_plus[a][b]
            q = Tr_minus[a][b]

            if p == 0.0 and q == 0.0:
                beta[a][b] = 0.0
            elif q == 0.0 and p > 0.0:
                beta[a][b] = float("inf")
            elif p == 0.0 and q > 0.0:
                beta[a][b] = float("-inf")
            else:
                beta[a][b] = log_base(p / q, 2.0)

            j += 1
        i += 1
    return beta


def score_sequence(seq, beta):
    seq = seq.strip().upper()
    s = 0.0
    k = 0
    while k < len(seq) - 1:
        a = seq[k]
        b = seq[k + 1]
        if a in VALID and b in VALID:
            s += beta[a][b]
        k += 1
    return s


def fmt_cell(v, digits):
    if v == float("inf"):
        return "+INF".rjust(10)
    if v == float("-inf"):
        return "-INF".rjust(10)
    return f"{v:10.{digits}f}"


def print_table(title, M, digits):
    w = 10
    line = "═" * (6 + (w + 1) * len(ALPHABET))
    print(line)
    print(title)
    print(line)

    header = "     "
    j = 0
    while j < len(ALPHABET):
        header += ALPHABET[j].rjust(w) + " "
        j += 1
    print(header)

    i = 0
    while i < len(ALPHABET):
        a = ALPHABET[i]
        row = a.rjust(3) + " |"
        j = 0
        while j < len(ALPHABET):
            b = ALPHABET[j]
            row += fmt_cell(M[a][b], digits) + " "
            j += 1
        print(row)
        i += 1

    print(line)
    print("")


def print_counts(title, counts, totals):
    line = "─" * 64
    print(line)
    print(title)
    print(line)
    i = 0
    while i < len(ALPHABET):
        a = ALPHABET[i]
        out = []
        j = 0
        while j < len(ALPHABET):
            b = ALPHABET[j]
            out.append(f"{b}:{counts[a][b]}")
            j += 1
        print(f"{a} -> " + ", ".join(out) + f"   (total {totals[a]})")
        i += 1
    print(line)
    print("")


def decide(score):
    if score > 0:
        return "CpG ISLAND (+) more likely"
    if score < 0:
        return "NON-ISLAND (-) more likely"
    return "No preference (score = 0)"


def main():
    S1 = "ATCGATTCGATATCATACACGTAT"
    S2 = "CTCGACTAGTATGAAGTCCACGCTTG"
    S  = "CAGGTTGGAAACGTAA"

    c_plus, t_plus = count_transitions(S1)
    c_minus, t_minus = count_transitions(S2)

    Tr_plus = transition_probs(c_plus, t_plus)
    Tr_minus = transition_probs(c_minus, t_minus)

    beta = log_likelihood_matrix(Tr_plus, Tr_minus)

    print("CpG MODEL COMPARISON (transition LLR)")
    print("")
    print_counts("Raw transition counts: CpG+ model (S1)", c_plus, t_plus)
    print_counts("Raw transition counts: CpG- model (S2)", c_minus, t_minus)

    print_table("Tr+ (probabilities from S1)", Tr_plus, digits=4)
    print_table("Tr- (probabilities from S2)", Tr_minus, digits=4)
    print_table("beta = log2(Tr+/Tr-)", beta, digits=4)

    sc = score_sequence(S, beta)
    print("Test sequence:", S)
    print(f"Total LLR score = {sc:.6f}")
    print("Decision:", decide(sc))


if __name__ == "__main__":
    main()
