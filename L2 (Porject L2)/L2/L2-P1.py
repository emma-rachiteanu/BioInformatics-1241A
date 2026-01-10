"""
Lab 2 problem 1

Make a brute force engine which is able to generate all the dinucleotides and trinucleotides combinations.
Search for each combination inside sequence S and calculate their relative frequencies.
S = "ATTGTCCCAATCTGTTG"
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox

def make_dinucs():
    nucs = ['A','C','G','T']
    d = []
    for x in nucs:
        for y in nucs:
            d.append(x+y)
    return d

def make_trinucs():
    n = ['A','C','G','T']
    t = []
    for a in n:
        for b in n:
            for c in n:
                t.append(a+b+c)
    return t

def count_pat(seq, pat):
    cnt = 0
    for i in range(len(seq)-len(pat)+1):
        if seq[i:i+len(pat)] == pat:
            cnt +=1
    return cnt

def percent_pat(seq, pat):
    pos = len(seq)-len(pat)+1
    if pos<=0:
        return 0
    return count_pat(seq, pat)/pos*100

class NucThing:
    def __init__(self, root):
        self.root = root
        self.root.title("kjqnswjk")
        self.root.geometry("800x600")
        self.seq_def = "ATTGTCCCAATCTGTTG"
        self.make_ui()

    def make_ui(self):
        ftop = ttk.Frame(self.root)
        ftop.pack(fill=tk.X, padx=8, pady=6)
        ttk.Label(ftop, text="DNA seq:").pack(side=tk.LEFT)
        self.e_seq = ttk.Entry(ftop, width=80)
        self.e_seq.pack(side=tk.LEFT,padx=6)
        self.e_seq.insert(0,self.seq_def)
        ttk.Button(ftop,text="Go",command=self.do_analysis).pack(side=tk.LEFT,padx=4)
        ttk.Button(ftop,text="Clear",command=self.clear_all).pack(side=tk.LEFT,padx=4)

        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill=tk.BOTH,expand=True,padx=8,pady=6)

        self.f2 = ttk.Frame(self.nb)
        self.nb.add(self.f2,text="Di (k=2)")
        self.t2 = scrolledtext.ScrolledText(self.f2,wrap=tk.NONE,font=("Consolas",10))
        self.t2.pack(fill=tk.BOTH,expand=True,padx=4,pady=4)

        self.f3 = ttk.Frame(self.nb)
        self.nb.add(self.f3,text="Tri (k=3)")
        self.t3 = scrolledtext.ScrolledText(self.f3,wrap=tk.NONE,font=("Consolas",10))
        self.t3.pack(fill=tk.BOTH,expand=True,padx=4,pady=4)

        self.fs = ttk.Frame(self.nb)
        self.nb.add(self.fs,text="Summary")
        self.ts = scrolledtext.ScrolledText(self.fs,wrap=tk.WORD,font=("Consolas",10))
        self.ts.pack(fill=tk.BOTH,expand=True,padx=4,pady=4)

        self.status = tk.StringVar()
        self.status.set("Ready")
        ttk.Label(self.root,textvariable=self.status,anchor="w").pack(fill=tk.X,padx=8,pady=(0,6))

    def clear_all(self):
        self.t2.delete(1.0,tk.END)
        self.t3.delete(1.0,tk.END)
        self.ts.delete(1.0,tk.END)
        self.status.set("Cleared")

    def do_analysis(self):
        s = self.e_seq.get().strip().upper()
        if not s:
            messagebox.showwarning("Oops","Enter DNA seq plz")
            return
        if not all(c in "ACGT" for c in s):
            messagebox.showerror("Err","Only A,C,G,T allowed")
            return
        self.status.set("Analyzing...")
        self.root.update_idletasks()
        self.clear_all()
        self.analyze_di(s)
        self.analyze_tri(s)
        self.make_summary(s)
        self.status.set(f"Done, len={len(s)}")

    def analyze_di(self, s):
        res=[]
        for d in make_dinucs():
            c = count_pat(s,d)
            p = percent_pat(s,d)
            res.append((d,c,p))
        res.sort(key=lambda x:(-x[2],x[0]))
        lines=[f"Len:{len(s)}","Windows k=2:"+str(max(len(s)-1,0)),"","Di   Count   %","-"*22]
        for d,c,p in res:
            lines.append(f"{d:<6} {c:>5} {p:>6.2f}")
        self.t2.insert("1.0","\n".join(lines))
        self.d_res=res

    def analyze_tri(self, s):
        res=[]
        for t in make_trinucs():
            c=count_pat(s,t)
            p=percent_pat(s,t)
            res.append((t,c,p))
        res.sort(key=lambda x:(-x[2],x[0]))
        lines=[f"Len:{len(s)}","Windows k=3:"+str(max(len(s)-2,0)),"","Tri    Count   %","-"*24]
        for t,c,p in res:
            lines.append(f"{t:<7} {c:>5} {p:>6.2f}")
        self.t3.insert("1.0","\n".join(lines))
        self.t_res=res

    def make_summary(self,s):
        lines=[f"Seq: {s}",f"Len: {len(s)}","","Bases:"]
        for b in 'ACGT':
            c=s.count(b)
            pct=(c/len(s)*100 if len(s)>0 else 0)
            lines.append(f"  {b}: {c} ({pct:.2f}%)")
        lines.append("\nTop 5 Di:")
        for i,(d,c,p) in enumerate(self.d_res[:5],1):
            lines.append(f"  {i}. {d} {c} ({p:.2f}%)")
        lines.append("\nTop 5 Tri:")
        for i,(t,c,p) in enumerate(self.t_res[:5],1):
            lines.append(f"  {i}. {t} {c} ({p:.2f}%)")
        self.ts.insert("1.0","\n".join(lines))

def main():
    r=tk.Tk()
    app=NucThing(r)
    r.mainloop()

if __name__=="__main__":
    main()
