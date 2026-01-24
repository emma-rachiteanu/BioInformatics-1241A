import json
import re
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from collections import defaultdict


def tokenize(text):
    return re.findall(r"[A-Za-z]+", text.lower())


def build_symbols(tokens):
    w2s = {}
    s2w = {}
    k = 1
    i = 0
    while i < len(tokens):
        w = tokens[i]
        if w not in w2s:
            sym = "W" + str(k)
            w2s[w] = sym
            s2w[sym] = w
            k += 1
        i += 1
    return w2s, s2w


def count_pairs(tokens, w2s):
    cnt = defaultdict(int)
    i = 0
    while i < len(tokens) - 1:
        a = w2s[tokens[i]]
        b = w2s[tokens[i + 1]]
        cnt[(a, b)] += 1
        i += 1
    return cnt


def prob_pairs(cnt):
    totals = defaultdict(int)
    for (a, b), c in cnt.items():
        totals[a] += c

    probs = {}
    for (a, b), c in cnt.items():
        probs[(a, b)] = c / totals[a] if totals[a] > 0 else 0.0
    return probs


class WordApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Exercise 3 — Word Transition Statistics")
        self.geometry("1240x800")

        self.topn_var = tk.IntVar(value=12)
        self.last_payload = None

        self._ui()
        self._fill_example()

    def _ui(self):
        top = ttk.Frame(self, padding=10)
        top.pack(fill=tk.X)

        ttk.Label(top, text="Top transitions:").pack(side=tk.LEFT)
        ttk.Spinbox(top, from_=1, to=200, textvariable=self.topn_var, width=6).pack(side=tk.LEFT, padx=6)

        ttk.Button(top, text="Analyze", command=self._analyze).pack(side=tk.LEFT, padx=10)
        ttk.Button(top, text="Save JSON…", command=self._save).pack(side=tk.LEFT)

        ttk.Separator(self).pack(fill=tk.X)

        main = ttk.Frame(self, padding=10)
        main.pack(fill=tk.BOTH, expand=True)

        left = ttk.Labelframe(main, text="Text input (editable)", padding=10)
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        right = ttk.Frame(main)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(10, 0))

        self.text_box = tk.Text(left, height=14, wrap="word")
        self.text_box.pack(fill=tk.BOTH, expand=True)

        self.stats = ttk.Label(left, text="", justify=tk.LEFT)
        self.stats.pack(anchor="w", pady=(10, 0))

        nb = ttk.Notebook(right)
        nb.pack(fill=tk.BOTH, expand=True)

        tab_map = ttk.Frame(nb, padding=10)
        tab_seq = ttk.Frame(nb, padding=10)
        tab_top = ttk.Frame(nb, padding=10)

        nb.add(tab_map, text="Word ↔ Symbol")
        nb.add(tab_seq, text="Symbol Sequence")
        nb.add(tab_top, text="Top Transitions")

        self.map_tree = ttk.Treeview(tab_map, columns=("sym", "word"), show="headings", height=12)
        self.map_tree.heading("sym", text="symbol")
        self.map_tree.heading("word", text="word")
        self.map_tree.column("sym", width=110, anchor="center")
        self.map_tree.column("word", width=220, anchor="w")
        self.map_tree.pack(fill=tk.BOTH, expand=True)

        self.seq_view = tk.Text(tab_seq, wrap="word")
        self.seq_view.pack(fill=tk.BOTH, expand=True)
        self.seq_view.config(state=tk.DISABLED)

        self.top_tree = ttk.Treeview(tab_top, columns=("from", "to", "count", "p", "pair"), show="headings", height=12)
        for c, w, a in [("from", 160, "w"), ("to", 160, "w"), ("count", 90, "center"), ("p", 90, "center"), ("pair", 120, "center")]:
            self.top_tree.heading(c, text=c)
            self.top_tree.column(c, width=w, anchor=a)
        self.top_tree.pack(fill=tk.BOTH, expand=True)

        self.bottom = tk.Text(self, height=6, wrap="word")
        self.bottom.pack(fill=tk.X, padx=10, pady=(0, 10))
        self.bottom.config(state=tk.DISABLED)

    def _fill_example(self):
        text = (
            "Lorem ipsum dolor sit amet consectetur adipiscing elit sed do eiusmod tempor incididunt "
            "ut labore et dolore magna aliqua. Ut enim ad minim veniam quis nostrud exercitation "
            "ullamco laboris nisi ut aliquip ex ea commodo consequat."
        )
        self.text_box.delete("1.0", tk.END)
        self.text_box.insert(tk.END, text)

    def _clear_tree(self, tree):
        for item in tree.get_children():
            tree.delete(item)

    def _analyze(self):
        text = self.text_box.get("1.0", tk.END).strip()
        if not text:
            messagebox.showerror("Missing input", "Please paste or type some text.")
            return

        tokens = tokenize(text)
        if len(tokens) < 2:
            messagebox.showerror("Not enough words", "Need at least 2 words to form transitions.")
            return

        w2s, s2w = build_symbols(tokens)
        cnt = count_pairs(tokens, w2s)
        probs = prob_pairs(cnt)

        self._clear_tree(self.map_tree)
        self._clear_tree(self.top_tree)

        items = list(w2s.items())
        items.sort(key=lambda x: int(x[1][1:]))

        i = 0
        while i < len(items):
            w, s = items[i]
            self.map_tree.insert("", "end", values=(s, w))
            i += 1

        seq = [w2s[w] for w in tokens]
        self.seq_view.config(state=tk.NORMAL)
        self.seq_view.delete("1.0", tk.END)

        per_line = 14
        i = 0
        while i < len(seq):
            self.seq_view.insert(tk.END, " → ".join(seq[i:i+per_line]) + "\n")
            i += per_line
        self.seq_view.config(state=tk.DISABLED)

        pairs = list(cnt.items())
        pairs.sort(key=lambda x: x[1], reverse=True)

        topn = int(self.topn_var.get())
        shown = pairs[:topn]

        i = 0
        while i < len(shown):
            (a, b), c = shown[i]
            p = probs[(a, b)]
            self.top_tree.insert("", "end", values=(s2w[a], s2w[b], c, f"{p:.4f}", f"{a}->{b}"))
            i += 1

        self.stats.config(text=f"Words: {len(tokens)}   Unique: {len(w2s)}   Unique transitions: {len(cnt)}")

        pretty_counts = {}
        pretty_probs = {}
        for (a, b), c in cnt.items():
            key = f"{s2w[a]} -> {s2w[b]} ({a}->{b})"
            pretty_counts[key] = c
        for (a, b), p in probs.items():
            key = f"{s2w[a]} -> {s2w[b]} ({a}->{b})"
            pretty_probs[key] = p

        self.last_payload = {
            "original_text": text,
            "text_length": len(text),
            "word_count": len(tokens),
            "unique_word_count": len(w2s),
            "word_to_symbol": w2s,
            "symbol_to_word": s2w,
            "symbol_sequence": seq,
            "transition_counts": pretty_counts,
            "transition_probabilities": pretty_probs,
        }

        self.bottom.config(state=tk.NORMAL)
        self.bottom.delete("1.0", tk.END)
        self.bottom.insert(tk.END, "Analysis complete.\n")
        self.bottom.insert(tk.END, "Change the text or Top transitions value, then click Analyze again.\n")
        self.bottom.config(state=tk.DISABLED)

    def _save(self):
        if not self.last_payload:
            messagebox.showinfo("Nothing to save", "Click Analyze first.")
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
    WordApp().mainloop()
