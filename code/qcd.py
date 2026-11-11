"""Benchmark on the Quantum Cluster Database, following PROTOCOL.md (tag prereg-idea10).

  python qcd.py energies {development|confirmation} MODEL   -> MACE energies (single point), saved group by group
  python qcd.py metrics {development|confirmation} MODEL    -> protocol metrics

Energies go to data/qcd/energies_{set}_{model}.csv. A run resumes where it stopped.
"""

import csv
import json
import os
import sys
import time

import numpy as np

QCD = "data/qcd/x/data_json/data.json"
OUT = "data/qcd"
SEED = 20261103
MAX_ISO = 10


def split_elements(records):
    elems = sorted({r["element_symbol"] for r in records})
    order = np.random.default_rng(SEED).permutation(len(elems))
    shuffled = [elems[i] for i in order]
    half = len(shuffled) // 2
    return set(shuffled[:half]), set(shuffled[half:])


def groups(which):
    records = list(json.load(open(QCD, encoding="utf-8")).values())
    dev, conf = split_elements(records)
    keep = dev if which == "development" else conf
    by = {}
    for r in records:
        if r["element_symbol"] in keep:
            by.setdefault((r["element_symbol"], int(r["n_atoms"])), []).append(r)
    out = {}
    for key, rs in by.items():
        rs = sorted(rs, key=lambda r: r["energy_dft"])[:MAX_ISO]
        e = [r["energy_dft"] for r in rs]
        if len(rs) >= 3 and max(e) - min(e) >= 0.05:
            out[key] = rs
    return out, sorted(keep)


def atoms_of(r):
    from ase import Atoms
    parts = r["structure_xyz"].split(";")
    n = int(parts[0])
    rows = [p.split() for p in parts[2:2 + n]]
    return Atoms([x[0] for x in rows], positions=[[float(v) for v in x[1:4]] for x in rows])


def energies(which, model):
    from mace.calculators import mace_mp

    path = f"{OUT}/energies_{which}_{model}.csv"
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            done = {(r["el"], int(r["n"])) for r in csv.DictReader(f)}
    gs, elems = groups(which)
    print(f"{which}: {len(elems)} elements, {len(gs)} groups, {sum(len(v) for v in gs.values())} isomers; "
          f"{len(done)} groups already done", flush=True)
    calc = mace_mp(model=model, default_dtype="float64", device="cpu")
    new = not os.path.exists(path)
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["el", "n", "cluster_id", "e_dft", "e_ml"])
        t0, k = time.time(), 0
        for key in sorted(gs, key=lambda x: (x[1], x[0])):
            if key in done:
                continue
            rows = []
            for r in gs[key]:
                a = atoms_of(r)
                a.calc = calc
                rows.append([key[0], key[1], r["cluster_id"], r["energy_dft"], a.get_potential_energy()])
            w.writerows(rows)
            f.flush()
            k += 1
            if k % 50 == 0:
                print(f"  {k} new groups ({time.time() - t0:.0f}s), last {key}", flush=True)
    print("done", flush=True)


def metrics(which, model):
    import pandas as pd
    from scipy.stats import spearmanr

    df = pd.read_csv(f"{OUT}/energies_{which}_{model}.csv")
    rows = []
    for (el, n), g in df.groupby(["el", "n"]):
        e_dft = g.e_dft.to_numpy() - g.e_dft.min()
        e_ml = g.e_ml.to_numpy() - g.e_ml.min()
        keep = e_dft <= 0.5
        rel = e_ml[keep] - e_ml[keep][np.argmin(e_dft[keep])]
        # (exploratory, not in the protocol) rank of the DFT minimum in the model's ordering: how many of the model's
        # candidates must be re-checked in DFT to find it
        k_needed = int(np.where(np.argsort(e_ml, kind="stable") == np.argmin(e_dft))[0][0]) + 1
        rows.append({"el": el, "n": n, "m": len(g), "k": k_needed, "top": int(np.argmin(e_ml) == np.argmin(e_dft)),
                     "chance": 1 / len(g), "rho": spearmanr(e_dft, e_ml).correlation,
                     "mae_mev_atom": 1000 * np.abs(rel - e_dft[keep]).mean() / n})
    r = pd.DataFrame(rows)
    r["size_band"] = pd.cut(r.n, [2, 10, 25, 55], labels=["3–10", "11–25", "26–55"])
    print(f"{which} {model}: {len(r)} groups, {r.el.nunique()} elements")
    print(f"  right minimum {r.top.mean():.1%} (chance {r.chance.mean():.1%}); mean Spearman {r.rho.mean():.2f}; "
          f"MAE {r.mae_mev_atom.mean():.1f} meV/atom")
    print(r.groupby("size_band", observed=True)[["top", "chance", "rho", "mae_mev_atom"]].mean().round(3).to_string())
    print("  (exploratory) DFT minimum among the model's first k candidates: " +
          ", ".join(f"k={k}: {(r.k <= k).mean():.0%}" for k in (1, 2, 3, 5)))
    per_el = r.groupby("el")[["top", "chance", "rho"]].mean().sort_values("top")
    print("  per element (worst and best):")
    print(pd.concat([per_el.head(5), per_el.tail(5)]).round(2).to_string())
    r.to_csv(f"{OUT}/metrics_{which}_{model}.csv", index=False)


if __name__ == "__main__":
    {"energies": energies, "metrics": metrics}[sys.argv[1]](sys.argv[2], sys.argv[3])
