"""LaTeX for the detection table and the per-plot appendix table, from extra_analyses.json.

No number in these tables is typed by hand. Usage: python3 figures/manuscript_tables.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
E = json.load(open(os.path.join(os.path.dirname(HERE), "results", "extra_analyses.json")))
MET = ("SAT", "AMS3D", "FF3D", "TreeLite3D", "WS")
NAME = {"SAT": "SegmentAnyTree", "AMS3D": "AMS3D", "FF3D": "ForestFormer3D",
        "TreeLite3D": "TreeLite3D$^{\\dagger}$", "WS": "Watershed"}
SHORT = {"SAT": "SAT", "AMS3D": "AMS3D", "FF3D": "FF3D", "TreeLite3D": "TL3D", "WS": "WS"}


def n(x):
    return f"{x:.3f}".replace("-", "$-$")


def ci(v):
    return f"[{n(v[0])}, {n(v[1])}]"


def detection_table():
    rows = []
    for g in ("Chamrousse", "Pellizzano", "Both sites"):
        P = E["pooled"][g]
        best = {c: max(MET, key=lambda m: P[m][c]) for c in ("recall", "precision", "f1", "excess")}
        b = lambda m, c: f"\\textbf{{{n(P[m][c])}}}" if best[c] == m else n(P[m][c])
        rows.append(f"\\midrule\n\\multicolumn{{8}}{{l}}{{\\textit{{{g}}}}}\\\\")
        for m in MET:
            v = P[m]
            rows.append(f"{NAME[m]} & {b(m, 'recall')} & {ci(v['recall_ci'])} & {b(m, 'precision')} & "
                        f"{b(m, 'f1')} & {ci(v['f1_ci'])} & {b(m, 'excess')} & {ci(v['excess_ci'])}\\\\")
    return "\n".join(rows)


def plot_table():
    rows = []
    for k, v in E["per_plot"].items():
        site, plot = k.split("/")
        lab = ("Cham. " if site == "cham" else "Pell. ") + plot
        f1s = " & ".join(n(v[m]["f1"]) for m in MET)
        rows.append(f"{lab} & {v['trees']} & {v['area_m2']:.0f} & {v['stems_per_ha']} & "
                    f"{v['returns_per_m2']:.1f} & {v['nearest_stem_m']:.2f} & {f1s}\\\\")
    avg = []
    for g in ("Chamrousse", "Pellizzano", "Both sites"):
        A = E["plot_average"][g]
        avg.append(f"\\multicolumn{{6}}{{l}}{{Plot average, {g if g != 'Both sites' else 'both sites'}}} & "
                   + " & ".join(n(A[m]["f1"]) for m in MET) + "\\\\")
    return "\n".join(rows) + "\n\\midrule\n" + "\n".join(avg)


if __name__ == "__main__":
    out = os.path.join(HERE, "tables")
    os.makedirs(out, exist_ok=True)
    open(os.path.join(out, "detection_rows.tex"), "w").write(detection_table() + "\n")
    open(os.path.join(out, "plot_rows.tex"), "w").write(plot_table() + "\n")
    print(detection_table()); print(); print(plot_table())
