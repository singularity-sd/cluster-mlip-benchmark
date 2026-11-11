"""Pilot (forecast T3v2-2026-0037): do "universal" interatomic potentials rank the isomers of isolated metal
clusters correctly?

Data: low-energy isomers of Ag_n, n = 3–12, optimised in DFT (B3PW91/def2-TZVP), with relative energies in
kcal/mol (Zenodo 10263658, CC0).

For every size n with at least 2 isomers:
- MACE-MP-0 energy in the DFT geometry (single point) and after relaxation with MACE;
- is the most stable isomer for MACE the same as in DFT? Spearman correlation between the orderings;
- displacement of the structure during relaxation (RMSD).

Usage (environment with mace-torch): python pilot_ag.py data/ag/datosAgn.xyz
"""

import sys

import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from scipy.stats import spearmanr

KCAL = 0.0433641  # eV per kcal/mol


def read_blocks(path):
    lines = open(path, encoding="utf-8").read().splitlines()
    i, out = 0, []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        e = float(lines[i + 1].split()[0])
        rows = [ln.split() for ln in lines[i + 2:i + 2 + n]]
        atoms = Atoms([r[0] for r in rows], positions=[[float(x) for x in r[1:4]] for r in rows])
        out.append((n, e, atoms))
        i += 2 + n
    return out


def rmsd(a, b):
    pa = a.get_positions() - a.get_positions().mean(0)
    pb = b.get_positions() - b.get_positions().mean(0)
    # optimal alignment (Kabsch), same atom order
    u, _, vt = np.linalg.svd(pb.T @ pa)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    return float(np.sqrt(((pb @ rot - pa) ** 2).sum(1).mean()))


def main():
    from mace.calculators import mace_mp

    model = sys.argv[2] if len(sys.argv) > 2 else "medium"
    calc = mace_mp(model=model, default_dtype="float64", device="cpu")
    print("model:", model)
    blocks = read_blocks(sys.argv[1])
    sizes = sorted({n for n, _, _ in blocks})
    tot_top, tot_sizes, rhos = 0, 0, []
    print("n  isomers  E_DFT(kcal/mol)  |  E_MACE single point  |  E_MACE relaxed  | relaxation RMSD (Å)")
    for n in sizes:
        iso = [(e, a) for m, e, a in blocks if m == n]
        if len(iso) < 2:
            continue
        e_dft = np.array([e for e, _ in iso])
        e_sp, e_rx, drift = [], [], []
        for _, a in iso:
            a = a.copy()
            a.calc = calc
            e_sp.append(a.get_potential_energy())
            r = a.copy()
            r.calc = calc
            BFGS(r, logfile=None).run(fmax=0.01, steps=500)
            e_rx.append(r.get_potential_energy())
            drift.append(rmsd(a, r))
        e_sp = (np.array(e_sp) - min(e_sp)) / KCAL
        e_rx = (np.array(e_rx) - min(e_rx)) / KCAL
        top = int(np.argmin(e_rx) == np.argmin(e_dft))
        rho = spearmanr(e_dft, e_rx).correlation if len(iso) > 2 else np.nan
        tot_top += top
        tot_sizes += 1
        if not np.isnan(rho):
            rhos.append(rho)
        print(f"{n:2d} {len(iso):5d}   {np.round(e_dft, 1)} | {np.round(e_sp, 1)} | {np.round(e_rx, 1)} | "
              f"{np.round(drift, 2)}  same minimum: {'yes' if top else 'NO'}", flush=True)
    print(f"\nSame most stable isomer (after MACE relaxation): {tot_top}/{tot_sizes}; "
          f"mean Spearman {np.mean(rhos):.2f} over {len(rhos)} sizes with ≥ 3 isomers")


if __name__ == "__main__":
    main()
