import json
import tkinter as tk
from tkinter import ttk, messagebox, filedialog


def init_counts(symbols):
    counts = {}
    i = 0
    while i < len(symbols):
        j = 0
        while j < len(symbols):
            counts[f"{symbols[i]}->{symbols[j]}"] = 0
            j += 1
        i += 1
    return counts


def transition_counts(seq, symbols):
    seq = seq.strip().upper()
    counts = init_counts(symbols)
    i = 0
    while i < len(seq) - 1:
        a = seq[i]
        b = seq[i + 1]
        key = f"{a}->{b}"
        if key in counts:
            counts[key] += 1
        i += 1
    return counts


def transition_probs(counts, symbols):
    totals = {s: 0 for s in symbols}
    for k, v in counts.items():
        left = k.split("->")[0]
        totals[left] += v

    probs = {}
    i = 0
    while i < len(symbols):
        a = symbols[i]
        j = 0
        while j < len(symbols):
            b = symbols[j]
            key = f"{a}->{b}"
            denom = totals[a]
            probs[key] = (counts[key] / denom) if denom > 0 else 0.0
            j += 1
        i += 1
    return probs, totals


class DNAApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Exercise 2 — DNA Transition Matrix")
        self.geometry("1120x760")

        self.symbols = ["A", "C", "G", "T"]
        self.last_payload = None

        self._build_ui()
        self._fill_example()
        self._setup_tables()

    def _build_ui(self):
        top = ttk.Frame(self, padding=10)
        top.pack(fill=tk.X)

        ttk.Label(top, text="DNA sequence (expected length: 50):").pack(side=tk.LEFT)
        ttk.Button(top, text="Compute", command=self._compute).pack(side=tk.LEFT, padx=10)
        ttk.Button(top, text="Save JSON…", command=self._save).pack(side=tk.LEFT)

        ttk.Separator(self).pack(fill=tk.X)

        main = ttk.Frame(self, padding=10)
        main.pack(fill=tk.BOTH, expand=True)

        left = ttk.Labelframe(main, text="Input", padding=10)
        left.pack(side=tk.LEFT, fill=tk.Y)

        right = ttk.Frame(main)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(10, 0))

        self.seq_box = tk.Text(left, width=38, height=6, wrap="word")
        self.seq_box.pack(fill=tk.X)

        self.info = ttk.Label(left, text="", justify=tk.LEFT)
        self.info.pack(anchor="w", pady=(10, 0))

        nb = ttk.Notebook(right)
        nb.pack(fill=tk.BOTH, expand=True)

        tab1 = ttk.Frame(nb, padding=10)
        tab2 = ttk.Frame(nb, padding=10)
        nb.add(tab1, text="Counts")
        nb.add(tab2, text="Probabilities")

        self.count_tree = ttk.Treeview(tab1, show="headings", height=12)
        self.count_tree.pack(fill=tk.BOTH, expand=True)

        self.prob_tree = ttk.Treeview(tab2, show="headings", height=12)
        self.prob_tree.pack(fill=tk.BOTH, expand=True)

        self.bottom = tk.Text(self, height=6, wrap="word")
        self.bottom.pack(fill=tk.X, padx=10, pady=(0, 10))
        self.bottom.config(state=tk.DISABLED)

    def _fill_example(self):
        seq = "ATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
        self.seq_box.delete("1.0", tk.END)
        self.seq_box.insert(tk.END, seq)

    def _setup_tables(self):
        cols = ["from"] + self.symbols
        self.count_tree.configure(columns=cols)
        self.prob_tree.configure(columns=cols)

        for tree in (self.count_tree, self.prob_tree):
            for c in cols:
                tree.heading(c, text=c)
                tree.column(c, width=120 if c != "from" else 80, anchor="center")

    def _clear_tree(self, tree):
        for item in tree.get_children():
            tree.delete(item)

    def _compute(self):
        seq = self.seq_box.get("1.0", tk.END).strip().upper()

        if len(seq) == 0:
            messagebox.showerror("Missing input", "Please enter a DNA sequence.")
            return

        bad = []
        i = 0
        while i < len(seq):
            if seq[i] not in self.symbols:
                bad.append(seq[i])
            i += 1

        if bad:
            messagebox.showerror("Invalid symbols", "Only A, C, G, T are allowed.")
            return

        if len(seq) != 50:
            if not messagebox.askyesno("Length differs", f"Your sequence length is {len(seq)}.\nContinue anyway?"):
                return

        counts = transition_counts(seq, self.symbols)
        probs, totals = transition_probs(counts, self.symbols)

        self._clear_tree(self.count_tree)
        self._clear_tree(self.prob_tree)

        i = 0
        while i < len(self.symbols):
            a = self.symbols[i]
            row_c = [a]
            row_p = [a]
            j = 0
            while j < len(self.symbols):
                b = self.symbols[j]
                key = f"{a}->{b}"
                row_c.append(str(int(counts[key])))
                row_p.append(f"{float(probs[key]):.5f}")
                j += 1
            self.count_tree.insert("", "end", values=row_c)
            self.prob_tree.insert("", "end", values=row_p)
            i += 1

        comp = {s: seq.count(s) for s in self.symbols}

        self.last_payload = {
            "sequence": seq,
            "sequence_length": len(seq),
            "nucleotide_counts": comp,
            "transition_counts": counts,
            "transition_probabilities": probs,
        }

        self.info.config(text=f"Length: {len(seq)}\nA:{comp['A']}  C:{comp['C']}  G:{comp['G']}  T:{comp['T']}")

        self.bottom.config(state=tk.NORMAL)
        self.bottom.delete("1.0", tk.END)
        self.bottom.insert(tk.END, "Computed transition counts and probabilities.\n")
        self.bottom.insert(tk.END, "Edit the sequence and click Compute again whenever you want.\n")
        self.bottom.config(state=tk.DISABLED)

    def _save(self):
        if not self.last_payload:
            messagebox.showinfo("Nothing to save", "Click Compute first.")
            return

        path = filedialog.asksaveasfilename(
            title="Save JSON",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if not path:
            return

        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(self.last_payload, f, indent=4)
            messagebox.showinfo("Saved", f"Saved to:\n{path}")
        except Exception as e:
            messagebox.showerror("Save error", str(e))


if __name__ == "__main__":
    DNAApp().mainloop()
