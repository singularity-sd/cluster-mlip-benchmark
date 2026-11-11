"""Pilot 3: same check as pilot_mosn.py (MACE single point in the DFT geometries) on other isomer populations
from Zenodo, all in the format "N / name energy(Ha) ... / atoms":
- MgPt_n (n = 2–12), PBE0/def2-TZVP (Zenodo 15792001);
- PtSi_n (n = 2–11) neutral, PBE0/def2-TZVP (Zenodo 19357337);
- MB12 (M = Li, K, Be, Ca, Sc, Ti, Fe, Cu, Zn, Al, P), functional not stated in the record (Zenodo 21035557);
- NaPb12 neutral, PBE-D4/def2-SVP (Zenodo 23763138).

Usage: python pilot_pop.py LABEL model FILE1.xyz [FILE2.xyz ...]
"""

import sys

import numpy as np
from scipy.stats import spearmanr

from pilot_mosn import read_population


def main():
    from mace.calculators import mace_mp

    label, model, files = sys.argv[1], sys.argv[2], sys.argv[3:]
    calc = mace_mp(model=model, default_dtype="float64", device="cpu")
    tops, rhos, maes, n_ok = 0, [], [], 0
    for path in files:
        pop = read_population(path)
        if len(pop) < 3:
            continue
        e_dft = np.array([e for e, _ in pop])
        e_dft -= e_dft.min()
        e_ml = []
        for _, a in pop:
            a.calc = calc
            e_ml.append(a.get_potential_energy())
        e_ml = np.array(e_ml) - min(e_ml)
        keep = e_dft <= 1.0
        rel = e_ml[keep] - e_ml[keep][np.argmin(e_dft[keep])]
        top = int(np.argmin(e_ml) == np.argmin(e_dft))
        rho = spearmanr(e_dft, e_ml).correlation
        mae = float(np.abs(rel - e_dft[keep]).mean())
        tops += top
        rhos.append(rho)
        maes.append(mae)
        n_ok += 1
        name = path.replace("\\", "/").split("/")[-1].replace("_population.xyz", "")
        print(f"  {name}: {len(pop)} isomers, minimum {'yes' if top else 'no'}, Spearman {rho:.2f}, MAE {mae:.3f} eV",
              flush=True)
    print(f"TOTAL {label} {model}: same minimum {tops}/{n_ok}; mean Spearman {np.nanmean(rhos):.2f}; "
          f"mean MAE {np.mean(maes):.3f} eV")


if __name__ == "__main__":
    main()
