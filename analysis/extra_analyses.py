"""Uncertainty, chance and per-plot results that complement main_evaluation.py.

1. Per-plot table: area of the scored region, trees, stems per hectare, returns per m2,
   median distance to the nearest stem, and recall, precision and F1 of each method at
   the 3.5 m tolerance.
2. Pooled recall, precision, F1 and excess over chance, with 95 % percentile intervals
   from resampling plots. For the two sites together the resampling is stratified by
   site: every resample keeps 7 Chamrousse plots and 15 Pellizzano plots.
3. Plot-averaged recall, precision and F1, where every plot weighs the same.
4. Paired differences, SegmentAnyTree minus each other method, in recall, F1 and excess,
   on the same resamples, with and without stratification by site.
5. The chance floor of each plot: as many detections as the method kept, thrown
   uniformly inside the scored region and matched in the same way (NSIM replicates).
   The floor is also broken down by crown stratum.
6. TreeLite3D with the clustering chosen on one site and applied to the other.

Usage: python3 analysis/extra_analyses.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
ROOT = os.environ.get("TREELEARN_ROOT") or os.path.dirname(BASE)
sys.path.insert(0, os.path.join(BASE, "modal"))
sys.path.insert(0, HERE)
from evaluation_core import SITES, region, results_path  # noqa: E402
from main_evaluation import MET, STRATA, count, greedy_pairs, metrics, ref_of  # noqa: E402

D = 3.5
NB = 2000
NSIM = 200
SEED = 0


def chance_plot(nd, reg, ref, dmax, rng, reps=NSIM):
    """Expected matches of a uniform detector, in total and by stratum."""
    out = {"TP": 0.0, **{s: 0.0 for s in STRATA}}
    if not nd:
        return out
    v = reg.vertices
    x0, x1, y0, y1 = v[:, 0].min(), v[:, 0].max(), v[:, 1].min(), v[:, 1].max()
    for _ in range(reps):
        got = []
        while len(got) < nd:
            c = np.c_[rng.uniform(x0, x1, nd * 3), rng.uniform(y0, y1, nd * 3)]
            got.extend(c[reg.contains_points(c)].tolist())
        c = np.array(got[:nd])
        pr = greedy_pairs(c[:, 0], c[:, 1], ref["x"], ref["y"], dmax)
        out["TP"] += len(pr)
        for _, j in pr:
            out[ref["stratum"][j]] += 1
    return {k: v / reps for k, v in out.items()}


def nearest_stem(ref):
    x, y = np.asarray(ref["x"]), np.asarray(ref["y"])
    dd = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
    np.fill_diagonal(dd, np.inf)
    return float(np.median(dd.min(axis=1)))


def pooled(rows, idx):
    tp, nd, nr, tc = (sum(rows[i][c] for i in idx) for c in ("TP", "ND", "NR", "TPc"))
    r, p, f = metrics(tp, nr, nd)
    return r, p, f, (tp - tc) / nr if nr else 0.0


def resamples(site_of, keys, stratified, rng):
    groups = {}
    for i, k in enumerate(keys):
        groups.setdefault(site_of[k], []).append(i)
    for _ in range(NB):
        if stratified:
            yield [g[j] for g in groups.values() for j in rng.integers(0, len(g), len(g))]
        else:
            yield list(rng.integers(0, len(keys), len(keys)))


def q(v):
    return [round(float(np.percentile(v, 2.5)), 3), round(float(np.percentile(v, 97.5)), 3)]


def main():
    data = json.load(open(os.path.join(ROOT, "all_detections.json")))
    keys = [f"{s}/{p}" for s, _, ps in SITES for p in ps]
    site_of = {f"{s}/{p}": n for s, n, ps in SITES for p in ps}
    rng = np.random.default_rng(SEED)

    # per plot counts and chance
    rows = {m: [] for m in MET}
    for m in MET:
        base = count(data, keys, m, D)
        for k, b in zip(keys, base):
            ref = ref_of(k)
            ch = chance_plot(b["ND"], region(ref), ref, D, rng)
            st = {s: dict(n=int((ref["stratum"] == s).sum()),
                          TP=sum(1 for j in b["ok"] if ref["stratum"][j] == s),
                          TPc=ch[s]) for s in STRATA}
            rows[m].append(dict(TP=b["TP"], ND=b["ND"], NR=b["NR"], TPc=ch["TP"], strata=st))
        print(f"chance done for {m}")

    out = {"settings": dict(tolerance=D, bootstrap=NB, chance_replicates=NSIM, seed=SEED)}

    # 1. per plot table
    out["per_plot"] = {}
    for i, k in enumerate(keys):
        ref = ref_of(k)
        dens = data["density"][k]
        reg_area = dens["area"]
        row = dict(site=site_of[k], area_m2=reg_area, trees=len(ref["x"]),
                   stems_per_ha=round(len(ref["x"]) / reg_area * 1e4),
                   returns_per_m2=dens["density"], nearest_stem_m=round(nearest_stem(ref), 2))
        for m in MET:
            r = rows[m][i]
            rr, pp, ff = metrics(r["TP"], r["NR"], r["ND"])
            row[m] = dict(recall=round(rr, 3), precision=round(pp, 3), f1=round(ff, 3),
                          detections=r["ND"], chance_recall=round(r["TPc"] / r["NR"], 3))
        out["per_plot"][k] = row

    # 2 and 3. pooled with intervals, and plot averages
    groups = [(n, [i for i, k in enumerate(keys) if site_of[k] == n]) for _, n, _ in SITES]
    groups.append(("Both sites", list(range(len(keys)))))
    out["pooled"] = {}
    out["plot_average"] = {}
    for gname, idx in groups:
        sub = [keys[i] for i in idx]
        out["pooled"][gname] = {}
        out["plot_average"][gname] = {}
        for m in MET:
            rw = [rows[m][i] for i in idx]
            r, p, f, x = pooled(rw, range(len(rw)))
            boots = {"recall": [], "precision": [], "f1": [], "excess": []}
            rr = np.random.default_rng(SEED)
            strat = gname == "Both sites"
            for s in resamples(site_of, sub, strat, rr):
                a = pooled(rw, s)
                for c, val in zip(boots, a):
                    boots[c].append(val)
            out["pooled"][gname][m] = dict(
                recall=round(r, 3), precision=round(p, 3), f1=round(f, 3), excess=round(x, 3),
                chance=round(r - x, 3), **{c + "_ci": q(v) for c, v in boots.items()},
                stratified=strat)
            per = [metrics(w["TP"], w["NR"], w["ND"]) for w in rw]
            out["plot_average"][gname][m] = {c: round(float(np.mean([t[j] for t in per])), 3)
                                            for j, c in enumerate(("recall", "precision", "f1"))}

    # 4. paired differences against SegmentAnyTree
    out["sat_minus"] = {}
    for gname, idx in groups:
        sub = [keys[i] for i in idx]
        out["sat_minus"][gname] = {}
        for m in MET:
            if m == "SAT":
                continue
            a_rows = [rows["SAT"][i] for i in idx]
            b_rows = [rows[m][i] for i in idx]
            ra, _, fa, xa = pooled(a_rows, range(len(idx)))
            rb, _, fb, xb = pooled(b_rows, range(len(idx)))
            res = dict(recall=round(ra - rb, 3), f1=round(fa - fb, 3), excess=round(xa - xb, 3))
            for strat in ((False, True) if gname == "Both sites" else (False,)):
                rr = np.random.default_rng(SEED)
                d = {"recall": [], "f1": [], "excess": []}
                for s in resamples(site_of, sub, strat, rr):
                    a, b = pooled(a_rows, s), pooled(b_rows, s)
                    d["recall"].append(a[0] - b[0]); d["f1"].append(a[2] - b[2])
                    d["excess"].append(a[3] - b[3])
                tag = "stratified" if strat else "unstratified"
                res[tag] = {c + "_ci": q(v) for c, v in d.items()}
            out["sat_minus"][gname][m] = res

    # 5. chance floor by stratum, both sites
    out["strata"] = {}
    for gname, idx in groups:
        out["strata"][gname] = {}
        for s in STRATA:
            n = sum(rows["SAT"][i]["strata"][s]["n"] for i in idx)
            out["strata"][gname][s] = {"n": n}
            for m in MET:
                tp = sum(rows[m][i]["strata"][s]["TP"] for i in idx)
                tc = sum(rows[m][i]["strata"][s]["TPc"] for i in idx)
                out["strata"][gname][s][m] = dict(recall=round(tp / n, 3), chance=round(tc / n, 3),
                                                  excess=round((tp - tc) / n, 3))

    # 6. TreeLite3D, clustering chosen on one site and applied to the other
    tl = json.load(open(results_path("treelite_evaluation.json")))
    sw = tl["sweep"]
    best_on = {s: max(sw, key=lambda t: sw[t][s]["f1"]) for s in ("Chamrousse", "Pellizzano")}
    cross = {"Chamrousse": sw[best_on["Pellizzano"]]["Chamrousse"],
             "Pellizzano": sw[best_on["Chamrousse"]]["Pellizzano"]}
    tp = sum(cross[s]["TP"] for s in cross); nd = sum(cross[s]["ND"] for s in cross)
    nr = sum(cross[s]["NR"] for s in cross)
    r, p, f = metrics(tp, nr, nd)
    out["treelite_leave_one_site_out"] = dict(
        chosen_on_chamrousse=best_on["Chamrousse"], chosen_on_pellizzano=best_on["Pellizzano"],
        both_sites=dict(recall=round(r, 3), precision=round(p, 3), f1=round(f, 3)))

    with open(os.path.join(ROOT, "extra_analyses.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("saved to extra_analyses.json")

    g = out["pooled"]["Both sites"]
    for m in MET:
        v = g[m]
        print(f"{m:11} R {v['recall']:.3f} {v['recall_ci']}  P {v['precision']:.3f} "
              f"{v['precision_ci']}  F1 {v['f1']:.3f} {v['f1_ci']}  excess {v['excess']:.3f} "
              f"{v['excess_ci']}")
    for m, v in out["sat_minus"]["Both sites"].items():
        print(f"SAT-{m:10} R {v['recall']:+.3f} F1 {v['f1']:+.3f} X {v['excess']:+.3f} | "
              f"unstrat {v['unstratified']} | strat {v['stratified']}")
    print(out["treelite_leave_one_site_out"])


if __name__ == "__main__":
    main()
