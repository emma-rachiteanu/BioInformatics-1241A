import json
import tkinter as tk
from tkinter import ttk, messagebox, filedialog


def mat_vec_mul(A, x):
    n = len(A)
    y = [0.0] * n
    i = 0
    while i < n:
        s = 0.0
        j = 0
        while j < n:
            s += float(A[i][j]) * float(x[j])
            j += 1
        y[i] = s
        i += 1
    return y


def predict(A, x0, steps):
    seq = []
    seq.append([float(v) for v in x0])
    cur = [float(v) for v in x0]
    k = 0
    while k < steps:
        cur = mat_vec_mul(A, cur)
        seq.append(cur)
        k += 1
    return seq


def is_square(A):
    if not A:
        return False
    n = len(A)
    i = 0
    while i < n:
        if len(A[i]) != n:
            return False
        i += 1
    return True


def parse_grid(entries):
    A = []
    i = 0
    while i < len(entries):
        row = []
        j = 0
        while j < len(entries[i]):
            row.append(float(entries[i][j].get().strip()))
            j += 1
        A.append(row)
        i += 1
    return A


def parse_vec(entries):
    v = []
    i = 0
    while i < len(entries):
        v.append(float(entries[i].get().strip()))
        i += 1
    return v


class MatrixPredictorApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Exercise 1 — Discrete Prediction (x[t+1] = A·x[t])")
        self.geometry("1100x740")

        self.n_var = tk.IntVar(value=3)
        self.steps_var = tk.IntVar(value=5)

        self._build_top()
        self._build_body()
        self._build_result()

        self._make_grid()

    def _build_top(self):
        top = ttk.Frame(self, padding=10)
        top.pack(fill=tk.X)

        ttk.Label(top, text="Matrix size n:").pack(side=tk.LEFT)
        ttk.Spinbox(top, from_=1, to=12, textvariable=self.n_var, width=5, command=self._make_grid).pack(side=tk.LEFT, padx=6)

        ttk.Label(top, text="Steps:").pack(side=tk.LEFT, padx=(18, 0))
        ttk.Spinbox(top, from_=1, to=50, textvariable=self.steps_var, width=5).pack(side=tk.LEFT, padx=6)

        ttk.Button(top, text="Rebuild input grid", command=self._make_grid).pack(side=tk.LEFT, padx=(18, 6))
        ttk.Button(top, text="Run prediction", command=self._run).pack(side=tk.LEFT, padx=6)
        ttk.Button(top, text="Save JSON…", command=self._save_json).pack(side=tk.LEFT, padx=6)

        ttk.Separator(self).pack(fill=tk.X)

    def _build_body(self):
        body = ttk.Frame(self, padding=10)
        body.pack(fill=tk.BOTH, expand=False)

        left = ttk.Labelframe(body, text="Matrix A", padding=10)
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 8))

        right = ttk.Labelframe(body, text="Initial vector x0", padding=10)
        right.pack(side=tk.LEFT, fill=tk.Y, padx=(8, 0))

        self.matrix_frame = ttk.Frame(left)
        self.matrix_frame.pack(fill=tk.BOTH, expand=True)

        self.vec_frame = ttk.Frame(right)
        self.vec_frame.pack(fill=tk.Y)

        hint = ttk.Label(right, text="Tip: edit values and press Run.\nDefaults are a valid example.", justify=tk.LEFT)
        hint.pack(anchor="w", pady=(12, 0))

    def _build_result(self):
        area = ttk.Labelframe(self, text="Results", padding=10)
        area.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 10))

        self.tree = ttk.Treeview(area, columns=(), show="headings", height=10)
        self.tree.pack(fill=tk.BOTH, expand=True)

        self.summary = tk.Text(area, height=6, wrap="word")
        self.summary.pack(fill=tk.X, pady=(10, 0))
        self.summary.config(state=tk.DISABLED)

        self.last_payload = None

    def _clear_matrix_widgets(self):
        for w in self.matrix_frame.winfo_children():
            w.destroy()
        for w in self.vec_frame.winfo_children():
            w.destroy()

    def _make_grid(self):
        try:
            n = int(self.n_var.get())
            if n < 1:
                raise ValueError
        except Exception:
            messagebox.showerror("Invalid n", "Matrix size must be a positive integer.")
            return

        self._clear_matrix_widgets()

        self.A_entries = []
        self.x_entries = []

        hdr = ttk.Frame(self.matrix_frame)
        hdr.grid(row=0, column=0, columnspan=n+1, sticky="w", pady=(0, 6))
        ttk.Label(hdr, text="A[i,j] values (float)").pack(side=tk.LEFT)

        j = 0
        while j < n:
            ttk.Label(self.matrix_frame, text=f"j={j+1}", width=8, anchor="center").grid(row=1, column=j+1, padx=2, pady=2)
            j += 1

        i = 0
        while i < n:
            ttk.Label(self.matrix_frame, text=f"i={i+1}", width=8, anchor="center").grid(row=i+2, column=0, padx=2, pady=2)
            row = []
            j = 0
            while j < n:
                e = ttk.Entry(self.matrix_frame, width=10)
                e.grid(row=i+2, column=j+1, padx=2, pady=2)
                row.append(e)
                j += 1
            self.A_entries.append(row)
            i += 1

        ttk.Label(self.vec_frame, text="x0 entries").grid(row=0, column=0, sticky="w", pady=(0, 6))
        i = 0
        while i < n:
            ttk.Label(self.vec_frame, text=f"x0[{i+1}]:").grid(row=i+1, column=0, sticky="e", padx=(0, 6), pady=2)
            e = ttk.Entry(self.vec_frame, width=14)
            e.grid(row=i+1, column=1, pady=2)
            self.x_entries.append(e)
            i += 1

        self._fill_example()

        self._setup_result_columns(n)

    def _fill_example(self):
        n = len(self.A_entries)
        example_A = [
            [0.70, 0.20, 0.10],
            [0.30, 0.40, 0.30],
            [0.20, 0.30, 0.50],
        ]
        example_x0 = [1.0, 0.0, 0.0]

        i = 0
        while i < n:
            j = 0
            while j < n:
                val = 0.0
                if n == 3:
                    val = example_A[i][j]
                else:
                    val = 1.0 if i == j else 0.0
                self.A_entries[i][j].delete(0, tk.END)
                self.A_entries[i][j].insert(0, str(val))
                j += 1
            i += 1

        i = 0
        while i < n:
            val = 0.0
            if n == 3:
                val = example_x0[i]
            else:
                val = 1.0 if i == 0 else 0.0
            self.x_entries[i].delete(0, tk.END)
            self.x_entries[i].insert(0, str(val))
            i += 1

    def _setup_result_columns(self, n):
        cols = ["t"]
        i = 0
        while i < n:
            cols.append(f"x[{i+1}]")
            i += 1
        cols.append("sum")

        self.tree.configure(columns=cols)

        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=110 if c != "t" else 60, anchor="center")

    def _run(self):
        try:
            A = parse_grid(self.A_entries)
            x0 = parse_vec(self.x_entries)
            if not is_square(A):
                raise ValueError("Matrix is not square.")
            if len(A) != len(x0):
                raise ValueError("Vector length must match matrix size.")
            steps = int(self.steps_var.get())
            if steps < 1:
                raise ValueError("Steps must be >= 1.")
        except Exception as e:
            messagebox.showerror("Input error", str(e))
            return

        states = predict(A, x0, steps)

        for item in self.tree.get_children():
            self.tree.delete(item)

        n = len(x0)
        t = 0
        while t < len(states):
            row = [t]
            j = 0
            s = 0.0
            while j < n:
                v = float(states[t][j])
                s += v
                row.append(f"{v:.6f}")
                j += 1
            row.append(f"{s:.6f}")
            self.tree.insert("", "end", values=row)
            t += 1

        self.last_payload = {
            "matrix_size": n,
            "steps": steps,
            "matrix_A": A,
            "initial_vector": x0,
            "state_vectors": states,
        }

        self.summary.config(state=tk.NORMAL)
        self.summary.delete("1.0", tk.END)
        self.summary.insert(tk.END, f"Computed {steps} step(s). Showing x0..x{steps}.\n")
        self.summary.insert(tk.END, "You can edit A / x0 / steps and run again.\n")
        self.summary.config(state=tk.DISABLED)

    def _save_json(self):
        if not self.last_payload:
            messagebox.showinfo("Nothing to save", "Run the prediction first.")
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
    MatrixPredictorApp().mainloop()
