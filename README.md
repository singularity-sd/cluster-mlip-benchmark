# Do universal machine-learned interatomic potentials find the ground state of isolated clusters? A preregistered benchmark

A preregistered benchmark of universal machine-learned interatomic potentials on isolated clusters of 55 elements.

**Paper:** [paper.pdf](paper.pdf) · **Page:** [singularitysd.com/research/clusters.html](https://singularitysd.com/research/clusters.html) · Preprint, not yet peer reviewed.

## Abstract

Universal machine-learned interatomic potentials (uMLIPs) trained on crystalline DFT data are increasingly used to search the energy landscape of nanoclusters. Their reliability for isolated clusters, whose low-coordinated atoms are absent from the training data, has not been measured systematically. We ask a simple question: given the low-energy isomers of a cluster computed with DFT, does the uMLIP identify the same ground state? We answer it with a preregistered benchmark on the Quantum Cluster Database (PBE, 55 elements, 3–55 atoms). The protocol was fixed before any comparison. We developed on 27 elements and confirmed once on 28 held-out elements (1 458 element–size groups, 13 102 isomers). \begin{itemize} \item MACE-MP-0 ranks the DFT ground state first in 33.3% of the groups (95% CI 28.5–38.0%), against 12.2% expected by chance. MACE-OMAT-0 does so in 41.8% of the groups. \item Errors on relative isomer energies are 10–17 meV/atom, comparable to the spacing between competing isomers. \item The ground state is among the five lowest-energy candidates of the model in 80% (MACE-MP-0) and 87% (MACE-OMAT-0) of the groups. Re-checking a handful of candidates with DFT is therefore a cheap and effective remedy. \end{itemize} A pilot on six independent public isomer datasets gives the same picture (34% against 18% by chance), including a case where the reference uses the training functional.

## How it was made

The plan of the test was fixed before looking at the confirmation data, and the test was run once.

## What is in this repository

| file | what it is |
|---|---|
| `code/` | full code, protocol and intermediate results |
| `fig_pilot.png` | figure |
| `fig_qcd.png` | figure |
| `figure_pilot.py` | draws the pilot figure |
| `main.tex` | LaTeX source of the paper |
| `paper.pdf` | the paper |

## Reproduce

The scripts in `code/` reproduce the numbers of the paper; each one explains at the top how to run it.

## Cite

```bibtex
@misc{singularity_cluster_mlip_benchmark_2026,
  title  = {Do universal machine-learned interatomic potentials find the ground state of isolated clusters? A preregistered benchmark},
  author = {Singularity},
  year   = {2026},
  url    = {https://github.com/singularity-sd/cluster-mlip-benchmark}
}
```

## License

Code: MIT (`LICENSE`). Paper, text and figures: CC BY 4.0 (`LICENSE-paper.md`). Third-party data keep their own licence where stated.

---

Part of [Singularity](https://singularitysd.com): forecasting where science goes next, then testing it.
