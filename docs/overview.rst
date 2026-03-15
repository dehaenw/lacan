Overview
========

LACAN (Leveraging Adjacent Co-occurrence of Atomic Neighborhoods) is a
statistical molecular filter and generative toolkit.  It scores molecules by
asking: *how likely is each bond, given the chemical environments on both sides
of it?*  A model trained on ~27 million bonds from ChEMBL assigns a pointwise
mutual information (PMI) score to every bond environment pair.  The
molecule-level score is derived from the minimum per-bond PMI, so a single
unusual bond is enough to flag a molecule.

Scoring model
-------------

For each bond the model computes two ECFP2-like atom environment hashes — one
per endpoint — and looks them up in a pre-built profile.  The PMI score is::

    observed  = pairs[(e1, e2)] / setsize
    expected  = (idx[e1] / setsize / 2) * (idx[e2] / setsize / 2)
    PMI       = observed / expected

The molecule-level score::

    score = min_PMI / (1 + min_PMI)

saturates toward 1.0 as the worst-bond PMI grows large, and approaches 0
when the worst bond is near zero.

Module overview
---------------

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Module
     - Purpose
   * - :mod:`lacan.lacan`
     - Core scoring: ``score_mol``, ``assess_per_bond``, profile I/O
   * - :mod:`lacan.mutate`
     - Atom-level mutations (40+ reaction SMARTS) — the EXPLOIT step
   * - :mod:`lacan.replace`
     - Coarse fragment swaps (ring / substituent / linker) — the EXPLORE step
   * - :mod:`lacan.breed`
     - Molecular crossover via fragment recombination
   * - :mod:`lacan.gen`
     - Random generation, corpus biasing, adaptive GA
   * - :mod:`lacan.protect`
     - SMARTS-based atom exclusion, bond protection, ``mol_cleaner``
   * - :mod:`lacan.decompose`
     - Molecule fragmentation; corpus building

Quick start
-----------

.. code-block:: python

    from rdkit import Chem
    from lacan import lacan, gen, mutate, replace

    profile = lacan.load_profile("chembl")

    mol = Chem.MolFromSmiles("CCCc1nn(C)c2c(=O)[nH]c(-c3ccccc3)nc12")
    score, info = lacan.score_mol(mol, profile)
    print(f"Score: {score:.3f}  bad bonds: {info['bad_bonds']}")

    # Generate drug-like molecules
    mols = gen.generate_filtered_molecules(profile, n_molecules=10, n_jobs=-1)

    # Optimise toward a scoring function
    def my_score(mols):
        return [lacan.score_mol(m, profile)[0] for m in mols]

    winners = gen.generate_optimized_molecules(my_score, profile,
                                               startN=20, generations=5)

Genetic algorithm
-----------------

:func:`~lacan.gen.generate_optimized_molecules` runs an adaptive GA that
balances exploration and exploitation each generation using two mechanisms:

**Smooth explore fraction**
    A float ``explore_fraction`` (0–1) controls the budget split between
    exploration arms (ring/substituent/linker replacement, scaffold decoration,
    crossover, random injection) and exploitation arms (atom-level mutation
    from :mod:`lacan.mutate`).  It shifts toward mutation when the population
    plateaus, and toward exploration on diversity collapse, decaying back to a
    user-set baseline otherwise.

**Per-operation Thompson Sampling bandit**
    Each operation is treated as an independent arm with a Beta posterior over
    its hit rate.  Budget is allocated proportionally to sampled weights each
    generation, so productive arms receive more budget while all arms remain
    explored.  Statistics can optionally persist across runs.

Results are collected in a :class:`~lacan.gen.HallOfFame` that retains the
all-time best diverse molecules with a Tanimoto diversity gate.

**Presets** — ``preset="ml"`` / ``"medium"`` / ``"docking"`` / ``"guacamol"``
provide sensible defaults for fast, medium, slow, and unlimited-budget scoring
functions respectively.  Individual parameters always override preset values.

See :func:`~lacan.gen.generate_optimized_molecules` for the full parameter
reference.
