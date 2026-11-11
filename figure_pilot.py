"""Pilot figure: share of (element, size) groups in which the model finds the DFT minimum, per family."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

fam = ["Ag$_n$\n(B3PW91)", "MoSn$_n$\n(PBE)", "MgPt$_n$\n(PBE0)", "PtSi$_n$\n(PBE0)", "MB$_{12}$\n(n/a)"]
groups = np.array([10, 14, 10, 10, 11])
found = {"MACE-MP-0": [5, 4, 1, 3, 6], "MACE-MPA-0": [6, 3, 2, 5, 5], "MACE-OMAT-0": [6, 5, 3, 5, 4]}
chance = np.array([3.2, 1.3 - 0.5, 1.2, 1.6, 2.7])  # MoSn: the group with 2 isomers (1/2) excluded

x = np.arange(len(fam))
w = 0.2
fig, ax = plt.subplots(figsize=(7, 3.4))
for k, (name, v) in enumerate(found.items()):
    ax.bar(x + (k - 1) * w, np.array(v) / groups, w, label=name)
ax.bar(x + 2 * w, chance / groups, w, color="lightgray", hatch="//", label="random pick")
ax.set_xticks(x + 0.5 * w)
ax.set_xticklabels(fam)
ax.set_ylabel("ground state found")
ax.set_ylim(0, 1)
ax.axhline(1, color="k", lw=0.5)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig("fig_pilot.png", dpi=200)
print("wrote fig_pilot.png")
