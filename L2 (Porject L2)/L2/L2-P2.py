"""
Lab 2 Problem 2

Design an application by using AI, which contains a GUI that allows the user to select a fasta file. The content
of the fasta file should be analyzed by using a sliding window of 30 positions. The content for each sliding window 
should be used in order to extract/compute the relative frequencies of the nucleoids formed in the aphabet of the 
sequence. The output of the app should be a chart containing 4 signals, one signal for each symbol for the alphabet of the sequence.
In order to plot this chart first, the vector that contains the relative frequencies for each symbol must be ccalculated
And only then plotted
"""

import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

ALPH = ("A","C","G","T")

# Hardcoded FASTA file path
FASTA_FILE = "sample.fasta"

def load_fasta_file(path):
    lines=[]
    with open(path,"r",encoding="utf-8") as f:
        for l in f:
            l=l.strip()
            if not l or l.startswith(">"):
                continue
            lines.append(l)
    return "".join(lines).upper()

def sliding(seq,k):
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]

def rel_freq(win):
    L=len(win)
    if L==0:
        return {b:0.0 for b in ALPH}
    cnt={b:0 for b in ALPH}
    for ch in win:
        if ch in cnt:
            cnt[ch]+=1
    return {b:cnt[b]/L for b in ALPH}

def freq_vectors(seq,k):
    vecs={b:[] for b in ALPH}
    for w in sliding(seq,k):
        f=rel_freq(w)
        for b in ALPH:
            vecs[b].append(f[b])
    return vecs

class FastaApp(ttk.Frame):
    def __init__(self, root):
        super().__init__(root)
        root.title("FASTA window freqs")
        root.geometry("900x600")
        self.pack(fill=tk.BOTH, expand=True)

        self.seq = ""
        self.vecs = None

        self.k_var = tk.StringVar(value="30")
        # FIX: pass self as master and use value=...
        self.status_var = tk.StringVar(master=self, value="Using hardcoded file: "+FASTA_FILE)

        self._make_controls()
        self._make_plot()

    def _make_controls(self):
        bar=ttk.Frame(self)
        bar.pack(fill=tk.X,padx=8,pady=6)

        ttk.Label(bar,text="Window k:").pack(side=tk.LEFT,padx=(0,4))
        ttk.Entry(bar,textvariable=self.k_var,width=6).pack(side=tk.LEFT)

        ttk.Button(bar,text="Analyze & Plot",command=self.analyze).pack(side=tk.LEFT,padx=8)
        ttk.Label(self,textvariable=self.status_var,anchor="w").pack(fill=tk.X,padx=8,pady=(0,6))

    def _make_plot(self):
        self.fig=Figure(figsize=(6.5,4.5),dpi=100)
        self.ax=self.fig.add_subplot(111)
        self.ax.set_xlabel("Window index")
        self.ax.set_ylabel("Rel freq")
        self.ax.set_title("A/C/G/T frequencies")

        self.canvas=FigureCanvasTkAgg(self.fig,master=self)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH,expand=True,padx=8,pady=8)

    def analyze(self):
        try:
            k=int(self.k_var.get())
            if k<=0: raise ValueError
        except:
            messagebox.showerror("Error","k must be positive integer")
            return

        try:
            self.seq=load_fasta_file(FASTA_FILE)
        except Exception as e:
            messagebox.showerror("Error",f"Cannot read FASTA:\n{e}")
            return

        if len(self.seq)<k:
            messagebox.showwarning("Too short",f"Seq len={len(self.seq)} < k={k}")
            self._plot_empty()
            return

        self.status_var.set(f"Computing vectors for k={k}...")
        self.update_idletasks()

        self.vecs=freq_vectors(self.seq,k)
        self._plot_vecs(self.vecs)
        self.status_var.set(f"Done: len={len(self.seq)}, windows={len(self.vecs['A'])}, k={k}")

    def _plot_empty(self):
        self.ax.clear()
        self.ax.set_xlabel("Window index")
        self.ax.set_ylabel("Rel freq")
        self.ax.set_title("A/C/G/T frequencies")
        self.canvas.draw_idle()

    def _plot_vecs(self,vecs):
        self.ax.clear()
        n=len(next(iter(vecs.values())))
        x=list(range(n))
        self.ax.plot(x,vecs["A"],label="A")
        self.ax.plot(x,vecs["C"],label="C")
        self.ax.plot(x,vecs["G"],label="G")
        self.ax.plot(x,vecs["T"],label="T")
        self.ax.set_ylim(0,1)
        self.ax.set_xlabel("Window index")
        self.ax.set_ylabel("Rel freq")
        self.ax.set_title("A/C/G/T frequencies")
        self.ax.grid(True,linestyle=":",linewidth=0.5)
        self.ax.legend(loc="upper right")
        self.canvas.draw_idle()


def main():
    root=tk.Tk()
    app=FastaApp(root)
    root.mainloop()

if __name__=="__main__":
    main()
