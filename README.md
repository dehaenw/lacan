# LACAN: Leveraging Adjacent Co-occurrence of Atomic Neighborhoods

LACAN is a cheminformatics toolkit for **scoring, mutating, and generating drug-like molecules** using a statistical model of chemical bond environments learned from ChEMBL. It is designed as a library for generative chemistry pipelines and includes an adaptive genetic algorithm that can optimise molecules toward any user-defined scoring function.

> *"All sorts of things in this world behave like mirrors."* — Jacques Lacan

📖 **Full documentation:** https://lacan.readthedocs.io/en/latest/

---

## How it works

For every bond in a molecule, LACAN computes a pair of **ECFP2-like atom environment identifiers**, one per endpoint. Each identifier encodes atomic number, degree, hydrogen count, formal charge, and ring type (none / non-aromatic / aromatic), hashed to a 32-bit integer. The two hashes form a *bond pair*.

A **profile** (e.g. `chembl.pickle`) stores, for a large training corpus:
- `idx`: how often each atom environment appears
- `pairs`: how often each bond-pair co-occurs
- `setsize`: total number of bonds seen

The **pointwise mutual information (PMI)** for a single bond:

```
observed  = pairs[(env1, env2)] / setsize
expected  = (idx[env1] / setsize / 2) × (idx[env2] / setsize / 2)
bond_PMI  = observed / expected
```

The molecule-level score uses the minimum per-bond PMI:

```
score = min_PMI / (1 + min_PMI)
```

A score near 0 means at least one bond is chemically unusual; near 1.0 means all bond environments are well-represented in the training data. Bonds below a threshold (default PMI < 0.05) are reported as `bad_bonds`.

---

## Quick start

```python
from rdkit import Chem
from lacan import lacan, gen

profile = lacan.load_profile("chembl")

# Score a molecule
mol = Chem.MolFromSmiles("CCCc1nn(C)c2c(=O)[nH]c(-c3ccccc3)nc12")
score, info = lacan.score_mol(mol, profile)
print(f"Score: {score:.3f}  bad bonds: {info['bad_bonds']}")

# Generate drug-like molecules
mols = gen.generate_filtered_molecules(profile, n_molecules=100, n_jobs=-1)

# Optimise toward any scoring function
def my_score(mols):
    return [lacan.score_mol(m, profile)[0] for m in mols]

winners = gen.generate_optimized_molecules(my_score, profile,
                                            startN=50, generations=20)
# Returns: list of (smiles, score) sorted best-first
```

### Build a custom profile

```python
from rdkit import Chem
from lacan.lacan import get_profile_for_mols

suppl = Chem.SmilesMolSupplier("my_molecules.smi", titleLine=False)
profile = get_profile_for_mols(suppl, profile_name="my_profile", n_jobs=-1)
# Saved to lacan/data/my_profile.pickle
# Reload with: lacan.load_profile("my_profile")
```

For the full module reference, GA parameters, corpus biasing, and protection API see the **[documentation](https://lacan.readthedocs.io/en/latest/)**.

---

## Example notebooks

Worked examples are in `lacan/example_notebooks/`:

| Notebook | Contents |
|---|---|
| `generate_molecules.ipynb` | Random generation, corpus biasing, running the GA |
| `optimize_from_mol.ipynb` | Lead optimisation from a seed molecule, pharmacophore protection, `mol_cleaner` |
| `mutating_molecules.ipynb` | Atom-level mutations and score filtering |
| `evaluate_bonds.ipynb` | Per-bond PMI scoring and visualisation |
| `median_molecules.ipynb` | Molecular crossover |
| `shape_optimize_vortioxetine.ipynb` | 3D shape-guided scaffold hopping with pharmacophore locking |

---

## Installation
```bash
pip install lacan
```
Installation is done via Pip. This package requires Python ≥ 3.9 and RDKit.

For installing from source:
```bash
git clone https://github.com/wdehaen/lacan.git
cd lacan
pip install .
```

---

## Running the tests

```bash
pip install pytest
pytest                              # full suite (~155 tests)
pytest tests/test_protect.py -v    # single module
```

---

## Citation

*Preprint coming soon.*

If you use LACAN in your research, please cite:

```
Dehaen W. (2026). LACAN: Leveraging Adjacent Co-occurrence of Atomic Neighborhoods
for molecular scoring and generation. [Preprint]
```
