"""Pilot 2: MoSn_n clusters (n = 3–16), populations of isomers optimised in **PBE**/def2-TZVP
(Zenodo 11390454, CC-BY 4.0). Same functional the MACE models are trained on, so no doubt about the functional.

For every n: relative PBE energies (from the file) and MACE energies (single point in the PBE geometry, then after
relaxation with MACE); same most stable isomer? Spearman; mean absolute error on the relative energies of the isomers
within 1 eV.

Usage: python pilot_mosn.py PBE_FOLDER model
"""

import glob
import os
import sys

import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from scipy.stats import spearmanr

HARTREE = 27.211386


def read_population(path):
    lines = open(path, encoding="utf-8").read().splitlines()
    i, out = 0, []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        head = lines[i + 1].split()
        e = float(head[1]) * HARTREE
        rows = [ln.split() for ln in lines[i + 2:i + 2 + n]]
        out.append((e, Atoms([r[0] for r in rows], positions=[[float(x) for x in r[1:4]] for r in rows])))
        i += 2 + n
    return out


def main():
    from mace.calculators import mace_mp

    model = sys.argv[2] if len(sys.argv) > 2 else "medium"
    calc = mace_mp(model=model, default_dtype="float64", device="cpu")
    print("model:", model)
    files = sorted(glob.glob(os.path.join(sys.argv[1], "MoSn*_population.xyz")),
                   key=lambda p: int(os.path.basename(p)[4:].split("_")[0]))
    tops, rhos, maes, rx_tops = 0, [], [], 0
    for path in files:
        pop = read_population(path)
        e_dft = np.array([e for e, _ in pop])
        e_dft -= e_dft.min()
        keep = e_dft <= 1.0
        e_sp, e_rx = [], []
        for _, a in pop:
            a = a.copy()
            a.calc = calc
            e_sp.append(a.get_potential_energy())
            r = a.copy()
            r.calc = calc
            BFGS(r, logfile=None).run(fmax=0.02, steps=300)
            e_rx.append(r.get_potential_energy())
        e_sp = np.array(e_sp) - min(e_sp)
        e_rx = np.array(e_rx) - min(e_rx)
        top = int(np.argmin(e_sp) == np.argmin(e_dft))
        rx_top = int(np.argmin(e_rx) == np.argmin(e_dft))
        rho = spearmanr(e_dft, e_sp).correlation
        rel_dft = e_dft[keep]
        rel_sp = e_sp[keep] - e_sp[keep][np.argmin(rel_dft)]
        mae = float(np.abs(rel_sp - rel_dft).mean())
        tops += top
        rx_tops += rx_top
        rhos.append(rho)
        maes.append(mae)
        n = len(pop[0][1])
        print(f"MoSn{n - 1}: {len(pop)} isomers, PBE range {e_dft.max():.2f} eV; same minimum (single point / "
              f"relaxed): {'yes' if top else 'no'} / {'yes' if rx_top else 'no'}; Spearman {rho:.2f}; "
              f"MAE relative energies (within 1 eV) {mae:.3f} eV", flush=True)
    print(f"\nTOTAL {model}: same minimum {tops}/{len(files)} (single point), {rx_tops}/{len(files)} (relaxed); "
          f"mean Spearman {np.mean(rhos):.2f}; mean MAE {np.mean(maes):.3f} eV")


if __name__ == "__main__":
    main()
