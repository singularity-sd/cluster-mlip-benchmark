# Protocol: universal potentials on the clusters of the Quantum Cluster Database

Written on 3 November 2026. Made final with the tag `prereg-idea10` after only an **inspection of the format** of
the file `data.json`. The file was downloaded from the QCD website; SHA-256 of the archive:
dc0be2a8933bf0f09670e8402ae0763541aba5b3d0484c0ea4ba6dc60d089da9.

**What was seen before the tag:**
- the fields (`cluster_id`, `element_symbol`, `n_atoms`, `energy_dft`, `energy_relative`, `structure_xyz`, …);
- 61,689 structures, 55 elements, 2,915 (element, size) groups, 2,856 with at least 3 structures;
- the fields of a single record (Ag10), including its absolute energy in the xyz header.

No energy was compared with a model.

## Question
Can "universal" interatomic potentials trained on PBE crystals rank the isomers of isolated clusters computed in PBE?
It is measured on many elements, to learn whether the error seen in the pilots is general.

## Data
- QCD (Manna et al., Sci. Data 2023, CC-BY 4.0): for every element and size from 3 to 55 atoms, the low-energy
  structures with their DFT energy.
- Only (element, size) groups with **at least 3 isomers** and an energy range of at least 0.05 eV are used.
- Reference energy: `energy_dft`, in the geometry `structure_xyz`.
- **At most 10 isomers per group:** the 10 with the lowest `energy_dft`. This is a limit on computing time, decided
  before seeing the energies.

## Split
- Elements sorted alphabetically and shuffled with `numpy.random.default_rng(20261103)`.
- **Development:** first 50%. Everything can be explored there: thresholds, models, relaxation or not, error analysis.
- **Confirmation:** second 50%, opened only once with the models and metrics fixed below.

## Models (fixed)
- **Primary:** MACE-MP-0 medium.
- **Secondary:** MACE-OMAT-0 medium.

Both in float64 on CPU. For every isomer the energy is a single point in the QCD geometry. MACE-MPA-0 and relaxation
are exploratory analyses, on the development elements only.

## Metrics
1. **Primary:** share of groups in which the model picks the same most stable isomer as the QCD ("right minimum").
2. Mean Spearman correlation between the orderings.
3. Mean absolute error on the relative energies of the isomers within 0.5 eV of the minimum, in meV/atom.
4. The same metrics per group of the periodic table and per size band (3–10, 11–25, 26–55 atoms).

## Hypotheses, written before seeing the data
- **H1:** for MACE-MP-0, the "right minimum" on the confirmation elements is below 60%. In the pilots it was around
  30%; a reliable model should exceed 90%.
- **H2:** MACE-OMAT-0 does not beat MACE-MP-0 by more than 15 percentage points on the "right minimum".
- **H3:** the error is larger for small clusters (3–10 atoms) than for large ones (26–55), which look more like
  crystals.

## Random reference
For every group, the probability of guessing the minimum at random is 1/number of isomers. The "right minimum" of a
model is compared with the mean of this probability over the groups.

## What is reported in any case
All the numbers on the confirmation, even if the hypotheses fail. No new models or thresholds after opening the
confirmation.
