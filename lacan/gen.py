"""
gen.py — Random molecule generation and adaptive genetic algorithm for LACAN.

This module has two responsibilities:

1. **Random generation** — assemble drug-like molecules from scratch by
   recursively sampling fragments from the corpus and joining them through
   their dummy-atom attachment points.

2. **Optimisation** — run an adaptive GA (:func:`generate_optimized_molecules`)
   that iteratively improves molecules toward a user-defined scoring function
   while maintaining population diversity.

Random generation
-----------------
:func:`generate_random_molecule` builds a molecule by starting from a randomly
chosen corpus fragment and recursively filling its open attachment points with
rings, linkers, or substituents sampled proportionally to their corpus
frequency.  The recursion terminates when all dummies are consumed.

:func:`generate_filtered_molecule` wraps the above with a rejection loop that
discards molecules outside the atom-count window or below the LACAN score
threshold.

:func:`generate_filtered_molecules` parallelises the above using
``multiprocessing.Pool``, giving each worker a unique seed offset so the
molecules are independent.

Corpus biasing
--------------
:func:`get_corpus_from_mols` / :func:`bias_corpus` build a merged corpus that
over-represents fragments from a reference set (e.g. known actives).  The
``ratio`` parameter controls the boost: at ratio=2 the reference fragments'
effective frequency is doubled relative to the background ChEMBL corpus.

The resulting corpus can be passed as ``fragcorpus=`` to any generation
function to steer sampling toward the desired chemotype.

Adaptive GA
-----------
See :func:`generate_optimized_molecules` for the full description.

GA design
---------
* A :class:`HallOfFame` collects the all-time best diverse molecules;
  the final result is simply the top-N from it.

* **scoring_budget** — a hard cap on molecules sent to the scoring function
  per generation.  Candidates are pre-screened with cheap property filters
  and MaxMin diversity sampling before scoring, keeping the budget tight even
  when raw operation output is large.

* **Per-operation Thompson Sampling bandit** — seven arms (ring replacement,
  linker replacement, substituent replacement, scaffold decoration, crossover,
  random injection, and atom-level mutation) each maintain a Beta posterior
  over their hit rate.  Budget is allocated proportionally to sampled weights,
  so productive arms get more budget while all arms remain explored.

* **Smooth explore_fraction** — a float in [0,1] controls the mix of
  exploration arms (ring/linker/sub-replace, decoration, crossover, random)
  vs exploitation arms (mutation).  On plateau it shifts toward mutation
  (highest empirical hit rate); on diversity collapse it shifts toward
  exploration.  It decays back toward a user-set baseline otherwise.

* **Presets** — ``preset="ml"`` / ``"medium"`` / ``"docking"`` / ``"guacamol"``
  set a coherent bundle of defaults appropriate for fast, medium, slow, and
  unlimited-budget scoring functions respectively.  Individual kwargs override
  preset values.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from lacan import lacan
from lacan import mutate, breed, decompose
from lacan import replace as replace_module
import random
import math
import csv
import os
from rdkit import DataStructs
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdFingerprintGenerator
import multiprocessing

MFPGEN = rdFingerprintGenerator.GetMorganGenerator(2, fpSize=4096)
rxn1 = rdChemReactions.ReactionFromSmarts("[*:0]-[#0:1].[*:2]-[#0:3]>>[*:0]-[*:2].[*:1]-[*:3]")
rxn2 = rdChemReactions.ReactionFromSmarts("[*:0]=[#0:1].[*:2]=[#0:3]>>[*:0]=[*:2].[*:1]=[*:3]")

_corpus_cache = {}  # keyed by min_count so different thresholds coexist

# ── Presets ───────────────────────────────────────────────────────────────────

# Each preset sets a coherent bundle of defaults.  Individual kwargs passed to
# generate_optimized_molecules() always override preset values.
_PRESETS = {
    "ml": dict(
        # Fast QSAR/ML model: 0.1–10 ms/mol. Score many molecules freely.
        # Large startN is cheap and gives a much better initial population.
        startN=1000,
        popsize=40,
        scoring_budget=200,
        generations=15,
        plateau_patience=3,
        base_explore_fraction=0.4,
        hof_size=100,
        hof_sim_threshold=0.55,
        pool_diversity_k=2,          # keep top 2*popsize, then diversity-select
        n_replacements=10,
        conservative=True,
    ),
    "medium": dict(
        # ADMET models, shape, ML ensemble: 0.1–2 s/mol.
        startN=40,
        popsize=25,
        scoring_budget=60,
        generations=12,
        plateau_patience=3,
        base_explore_fraction=0.45,
        hof_size=100,
        hof_sim_threshold=0.65,
        pool_diversity_k=2,
        n_replacements=8,
        conservative=True,
    ),
    "docking": dict(
        # Physics-based docking: 5–30 s/mol.  Tight budget is critical.
        startN=20,
        popsize=15,
        scoring_budget=20,
        generations=10,
        plateau_patience=2,
        base_explore_fraction=0.5,
        hof_size=50,
        hof_sim_threshold=0.60,
        pool_diversity_k=2,
        n_replacements=5,
        conservative=True,
    ),
    "guacamol": dict(
        # Unlimited scoring budget (GuacaMol benchmark).
        # Large population + many generations; mutation-dominant.
        # Pool is culled immediately after scoring each generation to prevent
        # unbounded growth (was causing pool=22k which slows _diverse_top_k).
        startN=1000,
        popsize=200,
        scoring_budget=30000,
        generations=25,
        plateau_patience=5,
        base_explore_fraction=0.15,
        hof_size=1000,
        hof_sim_threshold=0.55,
        pool_diversity_k=3,
        n_replacements=10,           # was 15 — reduces frag time ~30%
        conservative=True,
        maxmin_cap_factor=1,
    ),}



def load_corpus(min_count=200):
    """Load the ChEMBL fragment corpus from ``data/rls.csv``, with caching.

    The CSV is read only on the first call for a given ``min_count`` value;
    subsequent calls return the cached result immediately.

    Parameters
    ----------
    min_count : int
        Minimum occurrence count for a fragment to be included (default 200).

    Returns
    -------
    list of [smiles, count, degree, ftype, bonds]
    """
    if min_count in _corpus_cache:
        return _corpus_cache[min_count]
    _this_dir = os.path.dirname(__file__)
    _path = os.path.join(_this_dir, "data/rls.csv")
    _entries = []
    with open(_path, newline="") as _f:
        _reader = csv.reader(_f, delimiter=",", quotechar='"')
        for _e in _reader:
            if _e[0] != "smiles" and int(_e[1]) > min_count:
                _entries.append([_e[0], int(_e[1]), int(_e[2]), _e[3], _e[4]])
    _corpus_cache[min_count] = _entries
    return _entries


def _fragcount(corpus):
    """Return total occurrence counts per fragment type for *corpus*."""
    fc = {"Sub": 0, "Linker": 0, "Ring": 0}
    for e in corpus:
        fc[e[3]] += e[1]
    return fc


def generate_random_molecule(fragments=[], fragcorpus=None, needed_sample=[], bondstring="", mol=None):
    """Recursively assemble a random molecule from corpus fragments.

    Parameters
    ----------
    fragments     : list — accumulates SMILES of all fragments used so far
    fragcorpus    : list of corpus entries (default: ChEMBL rls.csv)
    needed_sample : list of lists — allowed fragment types for each open point
    bondstring    : str — bond types for the current open attachment points
    mol           : RDKit Mol — molecule built so far

    Returns
    -------
    (mol, fragments)
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()
    FRAGCOUNT = _fragcount(fragcorpus)
    if len(fragments) > 0:
        new_needed_sample = []
        new_bondstring = ""
        es = []
        for i, cns in enumerate(needed_sample):
            fragtype = random.choices(cns, [FRAGCOUNT[ft] for ft in cns])[0]
            fe = [e for e in fragcorpus if e[3] == fragtype and bondstring[i] == e[4][0]]
            es.append(random.choices([k for k in fe], [k[1]**0.666 for k in fe])[0])
        for i, e in enumerate(es):
            fragments.append(e[0])
            if bondstring[i] == "-":
                mol = rxn1.RunReactants((mol, Chem.MolFromSmiles(e[0])))[0][0]
            else:
                mol = rxn2.RunReactants((mol, Chem.MolFromSmiles(e[0])))[0][0]
            Chem.SanitizeMol(mol)
            if e[3] == "Ring" and e[2] > 1:
                new_needed_sample += (e[2] - 1) * [["Sub", "Linker"]]
                new_bondstring += e[4][1:]
            if e[3] == "Linker":
                new_needed_sample += (e[2] - 1) * [["Ring"]]
                new_bondstring += e[4][1:]
        if len(new_needed_sample) > 0:
            mol, fragments = generate_random_molecule(fragments, fragcorpus, new_needed_sample, new_bondstring, mol)
    else:
        e = random.choices([k for k in fragcorpus], [k[1]**0.666 for k in fragcorpus])[0]
        mol = Chem.MolFromSmiles(e[0])
        if e[3] == "Linker":
            mol, fragments = generate_random_molecule([e[0]], fragcorpus, e[2] * [["Ring"]], e[4], mol)
        elif e[3] == "Sub":
            mol, fragments = generate_random_molecule([e[0]], fragcorpus, [["Ring"]], e[4], mol)
        else:
            mol, fragments = generate_random_molecule([e[0]], fragcorpus, e[2] * [["Sub", "Linker", "Ring"]], e[4], mol)
    return mol, fragments


def generate_filtered_molecule(profile, fragcorpus=None, threshold=0.5, seed=None, min_atoms=14, max_atoms=35):
    """Generate one random molecule that passes the LACAN score and size filters.

    Parameters
    ----------
    profile   : LACAN profile dict
    fragcorpus: fragment corpus (default: ChEMBL rls.csv)
    threshold : minimum LACAN score to accept (default 0.5)
    seed      : int or None
    min_atoms : minimum heavy-atom count, exclusive (default 14)
    max_atoms : maximum heavy-atom count, exclusive (default 35)

    Returns
    -------
    RDKit Mol
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()
    score = 0
    if seed:
        random.seed(seed)
    while score < threshold:
        mol, _ = generate_random_molecule(fragcorpus=fragcorpus)
        score, _ = lacan.score_mol(mol, profile)
        if not (min_atoms < len(mol.GetAtoms()) < max_atoms):
            score = 0
    return mol


def generate_filtered_molecules(profile, fragcorpus=None, threshold=0.5, seed=None,
                                 min_atoms=14, max_atoms=35, n_jobs=1, n_molecules=10):
    """Generate multiple filtered molecules, optionally in parallel.

    Parameters
    ----------
    profile     : LACAN profile dict
    fragcorpus  : fragment corpus (default: ChEMBL rls.csv)
    threshold   : minimum LACAN score (default 0.5)
    seed        : int or None — base seed
    min_atoms   : minimum heavy-atom count, exclusive (default 14)
    max_atoms   : maximum heavy-atom count, exclusive (default 35)
    n_jobs      : int — parallel workers; 1 = sequential, -1 = all CPU cores
    n_molecules : int — number of molecules to generate (default 10)

    Returns
    -------
    list of RDKit Mol
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()
    if n_jobs == 1:
        mols = [generate_filtered_molecule(profile, fragcorpus, threshold, seed, min_atoms, max_atoms)
                for i in range(n_molecules)]
    else:
        if n_jobs == -1:
            n_jobs = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=n_jobs)
        mols = pool.starmap(generate_filtered_molecule,
                            [(profile, fragcorpus, threshold, seed + 823848 * i if seed else None, min_atoms, max_atoms)
                             for i in range(n_molecules)])
        pool.close()
        pool.join()
    return mols


def get_corpus_from_mols(mols, fragcorpus=None, ratio=1):
    """Merge a custom fragment corpus from reference molecules into the background corpus.

    Parameters
    ----------
    mols       : iterable of RDKit Mol objects
    fragcorpus : background corpus (default: ChEMBL rls.csv)
    ratio      : float — frequency multiplier for reference fragments (default 1)

    Returns
    -------
    list of corpus entries
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()
    custom_entries = decompose.get_corpus(mols)
    custom_fragcount = sum([e[1] for e in custom_entries])
    fragcount = sum([e[1] for e in custom_entries])
    mult = fragcount / custom_fragcount * ratio
    fragcounts = {e[0]: e[1] for e in fragcorpus}
    custom_fragcounts = {e[0]: e[1] for e in custom_entries}
    merged_entries = []
    for e in custom_entries:
        if e[0] in fragcounts:
            merged_entries.append([e[0]] + [int(e[1] * mult) + fragcounts[e[0]]] + e[2:])
        else:
            merged_entries.append([e[0]] + [int(e[1] * mult)] + e[2:])
    for e in fragcorpus:
        if e[0] not in custom_fragcounts:
            merged_entries.append(e)
    return merged_entries


def bias_corpus(mols, fragcorpus=None, ratio=2.0):
    """Build a fragment corpus biased toward the chemistry of the provided molecules.

    Parameters
    ----------
    mols       : iterable of RDKit Mol objects — reference molecules
    fragcorpus : background corpus (default: ChEMBL rls.csv)
    ratio      : float — frequency multiplier (default 2.0)

    Returns
    -------
    list of corpus entries
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()
    return get_corpus_from_mols(mols, fragcorpus=fragcorpus, ratio=ratio)


# ── GA helpers ────────────────────────────────────────────────────────────────

def _mean_diversity(smis):
    """Mean pairwise Tanimoto distance for a list of SMILES. Range [0,1]."""
    if len(smis) < 2:
        return 1.0
    fps = [MFPGEN.GetFingerprint(Chem.MolFromSmiles(s)) for s in smis]
    total, n = 0.0, 0
    for i in range(len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i+1:])
        total += sum(1 - s for s in sims)
        n += len(sims)
    return total / n if n > 0 else 1.0


def _sign(higher_is_better):
    """Return -1 if higher is better (GA minimises internally)."""
    return -1 if higher_is_better else 1


def _safe_score(scoring_function, mols):
    """Call scoring_function safely, returning 0.0 for None mols or exceptions."""
    scores = []
    valid_mols = []
    valid_indices = []
    for i, mol in enumerate(mols):
        if mol is not None:
            valid_mols.append(mol)
            valid_indices.append(i)
    if not valid_mols:
        return [0.0] * len(mols)
    try:
        raw = scoring_function(valid_mols)
    except Exception:
        raw = [0.0] * len(valid_mols)
    result = [0.0] * len(mols)
    for idx, score in zip(valid_indices, raw):
        try:
            result[idx] = float(score)
        except Exception:
            result[idx] = 0.0
    return result


def _basic_property_filter(mol, min_atoms=8, max_atoms=60, max_mw=700.0,
                             max_rotbonds=15, max_hbd=7):
    """Fast property-based pre-filter before expensive scoring.

    Returns True if the molecule passes all filters.
    """
    if mol is None:
        return False
    try:
        n = mol.GetNumAtoms()
        if not (min_atoms <= n <= max_atoms):
            return False
        mw = Descriptors.MolWt(mol)
        if mw > max_mw:
            return False
        rb = Descriptors.NumRotatableBonds(mol)
        if rb > max_rotbonds:
            return False
        hbd = Descriptors.NumHDonors(mol)
        if hbd > max_hbd:
            return False
        return True
    except Exception:
        return False


def _maxmin_sample(mols, k):
    """MaxMin diversity sampling: pick k molecules maximising pairwise diversity.

    Precomputes the full n×n Tanimoto distance matrix once, then runs the
    greedy MaxMin selection using numpy min/argmax over columns — no per-
    iteration BulkTanimotoSimilarity calls inside the loop.

    Time complexity: O(n²) for matrix build + O(n×k) for selection,
    but the selection loop is pure numpy with no Python overhead per candidate.

    Parameters
    ----------
    mols : list of RDKit Mol
    k    : int — number to select (if k >= len(mols), returns all)

    Returns
    -------
    list of RDKit Mol
    """
    import numpy as np

    n = len(mols)
    if k >= n:
        return list(mols)

    fps = [MFPGEN.GetFingerprint(m) for m in mols]

    # Build full n×n distance matrix in one pass.
    # BulkTanimotoSimilarity(fps[i], fps) returns row i of the similarity
    # matrix; we compute all rows and convert to distances in one go.
    sim_matrix = np.zeros((n, n), dtype=np.float32)
    for i, fp in enumerate(fps):
        row = DataStructs.BulkTanimotoSimilarity(fp, fps)
        sim_matrix[i] = row
    dist_matrix = 1.0 - sim_matrix  # distance = 1 - similarity

    # Seed: pick the molecule with the highest mean distance to all others
    # (the most peripheral point).
    mean_dist = dist_matrix.sum(axis=1) / (n - 1)
    selected = [int(np.argmax(mean_dist))]

    # min_dist_to_selected[i] = distance from mol i to its nearest selected mol
    min_dist_to_selected = dist_matrix[selected[0]].copy()

    selected_mask = np.zeros(n, dtype=bool)
    selected_mask[selected[0]] = True

    while len(selected) < k:
        # Zero out already-selected so they can't be picked again
        min_dist_to_selected[selected_mask] = -1.0
        best_idx = int(np.argmax(min_dist_to_selected))
        if min_dist_to_selected[best_idx] < 0:
            break
        selected.append(best_idx)
        selected_mask[best_idx] = True
        # Update: min distance to selected set = elementwise min of current
        # min-distances and distances to the newly added molecule
        np.minimum(min_dist_to_selected, dist_matrix[best_idx],
                   out=min_dist_to_selected)

    return [mols[i] for i in selected]


def _prescreened_candidates(raw_mols, budget, seen_smis,
                              min_atoms=8, max_atoms=60,
                              max_mw=700.0, max_rotbonds=15, max_hbd=7,
                              maxmin_cap_factor=10):
    """Filter and diversity-sample a raw candidate list down to scoring budget.

    Pipeline:
      1. Drop None and already-seen SMILES
      2. Apply cheap property filters
      3. Random subsample to ``budget * maxmin_cap_factor`` before MaxMin —
         caps the O(n²) MaxMin cost regardless of how large the raw pool is
      4. MaxMin diversity sampling to reach ``budget``

    Parameters
    ----------
    raw_mols          : list of RDKit Mol (may contain None)
    budget            : int — maximum number to return
    seen_smis         : set of SMILES already scored (mutated in place)
    maxmin_cap_factor : int — MaxMin input is capped at
                        ``budget * maxmin_cap_factor`` (default 10).
                        E.g. budget=200 → MaxMin sees at most 2000 candidates.
                        Raise if you want more diversity coverage at the cost
                        of speed; lower to speed up further.

    Returns
    -------
    list of (smiles, mol) tuples, length <= budget
    """
    # Step 1: deduplicate and property filter
    candidates = []
    for m in raw_mols:
        if not _basic_property_filter(m, min_atoms=min_atoms, max_atoms=max_atoms,
                                       max_mw=max_mw, max_rotbonds=max_rotbonds,
                                       max_hbd=max_hbd):
            continue
        try:
            smi = Chem.MolToSmiles(m)
        except Exception:
            continue
        if smi in seen_smis:
            continue
        candidates.append((smi, m))

    if not candidates:
        return []

    # Step 2: random subsample before MaxMin to bound O(n²) FP cost.
    # Shuffle first so the random cap is unbiased, then MaxMin picks the
    # most diverse subset from this pre-sampled pool.
    maxmin_input_cap = budget * maxmin_cap_factor
    if len(candidates) > maxmin_input_cap:
        candidates = random.sample(candidates, maxmin_input_cap)

    # Step 3: MaxMin diversity sampling down to budget
    if len(candidates) > budget:
        sampled_mols = _maxmin_sample([m for _, m in candidates], budget)
        sampled_smis = {Chem.MolToSmiles(m) for m in sampled_mols}
        candidates = [(smi, m) for smi, m in candidates if smi in sampled_smis]

    # Register all returned SMILES as seen
    for smi, _ in candidates:
        seen_smis.add(smi)

    return candidates[:budget]


def _fragment_moves_by_type_worker(args):
    """Top-level picklable worker for parallel fragment moves across parents."""
    mol_smi, profile, fragcorpus, n_replacements, conservative, min_atoms = args
    mol = Chem.MolFromSmiles(mol_smi)
    if mol is None:
        return {"ring": [], "linker": [], "sub": [], "decorate": []}
    return _fragment_moves_by_type(mol, profile, fragcorpus,
                                    n_replacements=n_replacements,
                                    conservative=conservative,
                                    min_atoms=min_atoms)


def _fragment_moves_all_parents(pool, profile, fragcorpus, n_replacements=5,
                                 conservative=True, n_jobs=1):
    """Apply fragment moves to all pool parents, in parallel if n_jobs != 1.

    Returns
    -------
    list of (parent_internal_score, by_type_dict) — one per parent
    """
    if n_jobs == 1 or len(pool) <= 1:
        out = []
        for smi, sc in pool:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                out.append((sc, {"ring": [], "linker": [], "sub": [], "decorate": []}))
            else:
                out.append((sc, _fragment_moves_by_type(mol, profile, fragcorpus,
                                                         n_replacements=n_replacements,
                                                         conservative=conservative)))
        return out

    actual_jobs = multiprocessing.cpu_count() if n_jobs == -1 else n_jobs
    worker_args = [(smi, profile, fragcorpus, n_replacements, conservative, 5)
                   for smi, _ in pool]
    with multiprocessing.Pool(processes=actual_jobs) as mp_pool:
        by_type_list = mp_pool.map(_fragment_moves_by_type_worker, worker_args)
    return [(sc, bt) for (_, sc), bt in zip(pool, by_type_list)]


def _fragment_moves_by_type(mol, profile, fragcorpus, n_replacements=5,
                             conservative=True, min_atoms=5, protect_smarts=None):
    """Apply each fragment operation separately and return per-operation outputs.

    Returns
    -------
    dict mapping operation name -> list of RDKit Mol
        Keys: "ring", "linker", "sub", "decorate"
    """
    results = {"ring": [], "linker": [], "sub": [], "decorate": []}
    op_map = {
        "ring":     replace_module.replace_ring,
        "linker":   replace_module.replace_linker,
        "sub":      replace_module.replace_substituent,
    }
    for name, op in op_map.items():
        try:
            out = op(mol, profile, score_threshold=0.0,
                     fragcorpus=fragcorpus, n_replacements=n_replacements,
                     conservative=conservative, protect_smarts=protect_smarts)
            results[name] = [m for m in out if m is not None and m.GetNumAtoms() >= min_atoms]
        except Exception:
            pass
    try:
        out = replace_module.decorate_scaffold(
            mol, profile, score_threshold=0.0,
            fragcorpus=fragcorpus, n_replacements=n_replacements,
            mode="Hydrogen", protect_smarts=protect_smarts)
        results["decorate"] = [m for m in out if m is not None and m.GetNumAtoms() >= min_atoms]
    except Exception:
        pass
    return results


def _fragment_moves(mol, profile, fragcorpus, n_replacements=5, conservative=True, min_atoms=5, protect_smarts=None):
    """Apply all fragment-level operations to a single molecule (flat list output).

    Kept for backwards compatibility with optimize_from_mol.
    """
    by_type = _fragment_moves_by_type(mol, profile, fragcorpus,
                                       n_replacements=n_replacements,
                                       conservative=conservative,
                                       min_atoms=min_atoms,
                                       protect_smarts=protect_smarts)
    new_mols = []
    for mlist in by_type.values():
        new_mols += mlist
    seen, out = set(), []
    for m in new_mols:
        smi = Chem.MolToSmiles(m)
        if smi not in seen:
            seen.add(smi)
            out.append(m)
    return out


# ── Thompson Sampling Bandit ──────────────────────────────────────────────────

class ThompsonBandit:
    """Beta-Bernoulli Thompson Sampling bandit over GA operations.

    Models each arm as a Bernoulli process: did this arm produce at least one
    improvement this generation?  Maintains Beta(α, β) posteriors over each
    arm's true hit rate.  Allocation weights are drawn from these posteriors
    each generation, so uncertain arms are still explored but arms with
    consistently higher hit rates get proportionally more budget over time.

    Thompson Sampling is well-suited here because most moves are non-productive
    (hit rate ~1–5%), so the Bernoulli posterior handles sparse signals natively.
    Budget allocation is proportional (soft) rather than winner-take-all, so all
    arms remain explored.  Statistics persist across runs (saved to JSON), letting
    the bandit accumulate meaningful signal over many optimisation campaigns.

    Prior α₀ / β₀
    --------------
    The prior encodes domain knowledge about each arm's expected hit rate,
    based on empirical observations from fragment-based drug design:

    - sub, mutate   : slightly higher prior hit rate (α=2, β=18 → ~10%)
    - ring, linker  : moderate (α=1, β=14 → ~7%)
    - decorate      : moderate (α=1, β=14 → ~7%)
    - crossover     : lower prior (α=1, β=24 → ~4%)
    - random        : lowest (α=1, β=49 → ~2%) — random injection rarely
      beats an optimised population directly

    These are weak priors (total pseudo-counts 15–25) so they are quickly
    overridden by real data once enough generations accumulate.

    Parameters
    ----------
    arm_names   : list of str
    stats_file  : str or None — path to JSON file for persistent statistics.
                  If provided, α/β counts are loaded on init and saved after
                  each update.  Pass None to disable persistence.
    prior_alpha : dict mapping arm name -> α₀  (default: built-in domain priors)
    prior_beta  : dict mapping arm name -> β₀  (default: built-in domain priors)
    """

    # Domain-knowledge priors: (alpha0, beta0) per arm
    _DEFAULT_PRIORS = {
        "ring":      (1,  14),   # ~7%  hit rate prior
        "linker":    (1,  14),   # ~7%
        "sub":       (2,  18),   # ~10%
        "decorate":  (1,  14),   # ~7%
        "crossover": (1,  24),   # ~4%
        "random":    (1,  49),   # ~2%
        "mutate":    (2,  18),   # ~10%
    }

    def __init__(self, arm_names, stats_file=None, prior_alpha=None, prior_beta=None):
        self.arms = arm_names
        self.stats_file = stats_file

        # Build priors — per-arm overrides, falling back to defaults
        pa = prior_alpha or {}
        pb = prior_beta  or {}
        self._alpha0 = {a: pa.get(a, self._DEFAULT_PRIORS.get(a, (1, 9))[0])
                        for a in arm_names}
        self._beta0  = {a: pb.get(a, self._DEFAULT_PRIORS.get(a, (1, 9))[1])
                        for a in arm_names}

        # Posterior counts (accumulated across runs if stats_file is used)
        self.alpha = {a: self._alpha0[a] for a in arm_names}  # successes + prior
        self.beta  = {a: self._beta0[a]  for a in arm_names}  # failures  + prior
        self.n_pulls = {a: 0 for a in arm_names}

        # Load persisted stats if available
        if stats_file is not None:
            self._load()

    # ── persistence ──────────────────────────────────────────────────────────

    def _load(self):
        """Load α/β counts from stats_file if it exists."""
        import json
        if not os.path.exists(self.stats_file):
            return
        try:
            with open(self.stats_file) as f:
                data = json.load(f)
            for a in self.arms:
                if a in data.get("alpha", {}):
                    self.alpha[a] = data["alpha"][a]
                if a in data.get("beta", {}):
                    self.beta[a] = data["beta"][a]
                if a in data.get("n_pulls", {}):
                    self.n_pulls[a] = data["n_pulls"][a]
        except Exception as e:
            pass  # silently ignore corrupt files

    def save(self):
        """Persist current α/β counts to stats_file."""
        import json
        if self.stats_file is None:
            return
        try:
            os.makedirs(os.path.dirname(os.path.abspath(self.stats_file)), exist_ok=True)
            with open(self.stats_file, "w") as f:
                json.dump({"alpha": self.alpha, "beta": self.beta,
                           "n_pulls": self.n_pulls}, f, indent=2)
        except Exception:
            pass

    def reset_stats(self):
        """Reset posteriors to priors and delete stats_file if present."""
        self.alpha   = {a: self._alpha0[a] for a in self.arms}
        self.beta    = {a: self._beta0[a]  for a in self.arms}
        self.n_pulls = {a: 0 for a in self.arms}
        if self.stats_file and os.path.exists(self.stats_file):
            try:
                os.remove(self.stats_file)
            except Exception:
                pass

    # ── bandit interface ──────────────────────────────────────────────────────

    def sample_weights(self):
        """Draw one sample from each arm's Beta posterior.

        Returns a dict arm -> sampled weight.  Use these as proportional
        allocation weights for the scoring budget.
        """
        return {a: random.betavariate(self.alpha[a], self.beta[a])
                for a in self.arms}

    def update(self, arm, success):
        """Update posterior for *arm*.

        Parameters
        ----------
        arm     : str
        success : bool or int — True/1 if the arm produced at least one
                  improvement this generation, False/0 otherwise.
        """
        self.n_pulls[arm] += 1
        if success:
            self.alpha[arm] += 1
        else:
            self.beta[arm] += 1

    @property
    def hit_rate(self):
        """Posterior mean hit rate per arm: α / (α + β)."""
        return {a: self.alpha[a] / (self.alpha[a] + self.beta[a])
                for a in self.arms}

    def summary(self):
        hr = self.hit_rate
        rows = []
        for a in self.arms:
            rows.append(
                f"  {a:12s}  pulls={self.n_pulls[a]:4d}"
                f"  alpha={self.alpha[a]:6.1f}  beta={self.beta[a]:6.1f}"
                f"  hit_rate={hr[a]:.3f}"
            )
        return "\n".join(rows)


# ── HallOfFame ────────────────────────────────────────────────────────────────

class HallOfFame:
    """Diversity-gated collection of the best molecules seen across all generations.

    A new molecule is accepted only if its Tanimoto similarity to every current
    member is below ``sim_threshold``.  When full, a new molecule displaces the
    worst member if it scores better.

    Parameters
    ----------
    max_size      : int — maximum number of entries (default 100)
    sim_threshold : float — max allowed Tanimoto similarity to existing members
                    (default 0.70)
    higher_is_better : bool (default True)
    """

    # Number of nearest neighbours checked for the diversity gate.
    # Checking all members becomes expensive when HoF is large and adds little
    # value — a mediocre early entry shouldn't block a great later one forever
    # just because it happens to be similar.  We check the K most similar
    # existing members instead.
    _DIVERSITY_GATE_K = 20

    def __init__(self, max_size=100, sim_threshold=0.55, higher_is_better=True):
        self.max_size = max_size
        self.sim_threshold = sim_threshold
        self.higher_is_better = higher_is_better
        self._entries = []   # list of (smiles, score, fp)

    def _is_better(self, a, b):
        return a > b if self.higher_is_better else a < b

    def offer(self, smi, score, mol):
        """Try to add (smi, score) to the HallOfFame.

        Diversity gate: the candidate is accepted only if its maximum Tanimoto
        similarity to the ``_DIVERSITY_GATE_K`` most similar existing members
        is below ``sim_threshold``.  Using the top-K rather than all members
        prevents mediocre early entries from permanently blocking better later
        ones that happen to be structurally related.

        Returns True if accepted, False if rejected.
        """
        try:
            fp = MFPGEN.GetFingerprint(mol)
        except Exception:
            return False

        if self._entries:
            existing_fps = [e[2] for e in self._entries]
            sims = DataStructs.BulkTanimotoSimilarity(fp, existing_fps)
            # Check only against the K most similar members
            top_k_sims = sorted(sims, reverse=True)[:self._DIVERSITY_GATE_K]
            if top_k_sims and top_k_sims[0] >= self.sim_threshold:
                # Similar neighbour found — only block if that neighbour is at
                # least as good as the candidate (don't block an improvement)
                most_similar_idx = int(max(range(len(sims)), key=lambda i: sims[i]))
                neighbour_score = self._entries[most_similar_idx][1]
                if not self._is_better(score, neighbour_score):
                    return False
                # Candidate is better than its nearest neighbour → replace it
                self._entries[most_similar_idx] = (smi, score, fp)
                return True

        if len(self._entries) < self.max_size:
            self._entries.append((smi, score, fp))
            return True

        # HoF full and no similar neighbour to replace: displace the worst
        # member if the candidate is better than it.
        worst_idx = min(range(len(self._entries)),
                        key=lambda i: self._entries[i][1] if self.higher_is_better
                        else -self._entries[i][1])
        if self._is_better(score, self._entries[worst_idx][1]):
            self._entries[worst_idx] = (smi, score, fp)
            return True
        return False

    def _worth_offering(self, score):
        """Quick check: is this score worth running the full offer() logic?

        Returns True if the HoF isn't full, or if the score beats the current
        worst member.  Avoids fingerprint computation for clearly hopeless cases.
        """
        if len(self._entries) < self.max_size:
            return True
        worst_score = min(e[1] for e in self._entries) if self.higher_is_better \
                      else max(e[1] for e in self._entries)
        return self._is_better(score, worst_score)

    def top(self, n=None):
        """Return top-n (smiles, score) pairs sorted best-first."""
        sorted_entries = sorted(self._entries, key=lambda e: e[1],
                                reverse=self.higher_is_better)
        result = [(smi, sc) for smi, sc, _ in sorted_entries]
        return result if n is None else result[:n]

    @property
    def best_score(self):
        if not self._entries:
            return None
        return max(e[1] for e in self._entries) if self.higher_is_better \
               else min(e[1] for e in self._entries)

    def __len__(self):
        return len(self._entries)


# ── Pool diversity-aware selection ───────────────────────────────────────────

def _diverse_top_k(pool_entries, k, diversity_weight=0.3):
    """Select k molecules from pool balancing score and diversity.

    Algorithm: iteratively pick the candidate with the best combined
    score_rank + diversity_rank (where diversity = min distance to already
    selected set), starting from the top scorer.

    Parameters
    ----------
    pool_entries    : list of (smiles, internal_score)  [lower internal = better]
    k               : int — how many to keep
    diversity_weight: float in [0,1] — weight of diversity vs score (default 0.3)

    Returns
    -------
    list of (smiles, internal_score)
    """
    if len(pool_entries) <= k:
        return list(pool_entries)

    # Pre-compute fingerprints
    fps = []
    valid = []
    for entry in pool_entries:
        try:
            fp = MFPGEN.GetFingerprint(Chem.MolFromSmiles(entry[0]))
            fps.append(fp)
            valid.append(entry)
        except Exception:
            pass

    n = len(valid)
    if n <= k:
        return valid

    # Rank by score (ascending internal score = better)
    score_rank = [0] * n
    for rank, idx in enumerate(sorted(range(n), key=lambda i: valid[i][1])):
        score_rank[idx] = rank

    selected_indices = []
    selected_fps = []
    min_dist_to_selected = [1.0] * n  # initialise: no selected yet → max distance

    # Always start with the best scorer
    best_idx = min(range(n), key=lambda i: valid[i][1])
    selected_indices.append(best_idx)
    selected_fps.append(fps[best_idx])
    sims = DataStructs.BulkTanimotoSimilarity(fps[best_idx], fps)
    min_dist_to_selected = [min(1 - s, d) for s, d in zip(sims, min_dist_to_selected)]

    while len(selected_indices) < k:
        # Rank by diversity (higher min-distance = more diverse = better)
        remaining = [i for i in range(n) if i not in selected_indices]
        div_scores = [min_dist_to_selected[i] for i in remaining]
        max_div = max(div_scores) if max(div_scores) > 0 else 1.0
        div_rank = {i: div_scores[j] / max_div for j, i in enumerate(remaining)}

        # Combined score: lower score_rank_normalised + higher div is better
        max_srank = n - 1
        best = max(remaining,
                   key=lambda i: (1 - diversity_weight) * (1 - score_rank[i] / max_srank)
                                 + diversity_weight * div_rank[i])

        selected_indices.append(best)
        selected_fps.append(fps[best])
        sims = DataStructs.BulkTanimotoSimilarity(fps[best], fps)
        min_dist_to_selected = [min(1 - s, d) for s, d in zip(sims, min_dist_to_selected)]

    return [valid[i] for i in selected_indices]


def generate_optimized_molecules(
        scoring_function,
        profile,
        seed=123,
        startN=None,
        generations=None,
        popsize=None,
        scoring_budget=None,
        higher_is_better=True,
        diversity_threshold=0.35,
        plateau_patience=None,
        base_explore_fraction=None,
        explore_ratio=None,          # legacy alias for base_explore_fraction
        hof_size=None,
        hof_sim_threshold=None,
        sim_threshold=None,          # legacy alias for hof_sim_threshold
        pool_diversity_k=2,
        n_replacements=None,
        conservative=None,
        fragcorpus=None,
        n_jobs=-1,
        quiet=False,
        callback=None,
        seed_mols=None,
        preset="ml",
        bandit_stats_file=None,
        maxmin_cap_factor=None,      # cap for MaxMin prescreening; 1 = skip MaxMin (fast scorers)
        # Legacy / ignored parameters kept for API compatibility
        win_threshold=None,          # silently ignored — use HallOfFame instead
        **kwargs):
    """
    Adaptive GA with per-operation Thompson Sampling bandit, scoring budget, and HallOfFame.

    Quick start
    -----------
    For most use cases, set ``preset=`` and override individual parameters as
    needed::

        # Fast ML model
        winners = gen.generate_optimized_molecules(rf_score, profile, preset="ml")

        # Slow docking function — tight budget
        winners = gen.generate_optimized_molecules(dock, profile, preset="docking",
                                                    scoring_budget=15)

        # GuacaMol benchmark — unlimited scoring budget
        winners = gen.generate_optimized_molecules(guac_score, profile, preset="guacamol")

    Presets
    -------
    ``preset="ml"``        — fast QSAR/ML models (0.1–10 ms/mol)
    ``preset="medium"``    — ADMET, shape, ML ensembles (0.1–2 s/mol)
    ``preset="docking"``   — physics-based docking (5–30 s/mol)
    ``preset="guacamol"``  — unlimited scoring budget, large population,
    mutation-dominant allocation

    All preset values can be overridden by passing explicit kwargs.

    Key parameters
    --------------
    scoring_budget : int
        Hard cap on molecules sent to the scoring function per generation.
        Candidates are pre-screened with cheap property filters and MaxMin
        diversity sampling before scoring.  This is the most important knob
        for slow scoring functions.

    base_explore_fraction : float in [0,1]
        Baseline fraction of the scoring budget allocated to exploration arms
        (ring/linker/sub replacement, decoration, crossover, random injection)
        vs exploitation arms (atom-level mutation).  On scoring plateau this
        shifts toward mutation (highest empirical hit rate); on diversity
        collapse it shifts toward exploration.

    hof_size : int
        Maximum size of the HallOfFame.  The HoF collects the best diverse
        molecules found across all generations; the final result is the top-N
        from it sorted best-first.

    hof_sim_threshold : float
        Maximum Tanimoto similarity allowed between two HallOfFame members
        (diversity gate).  Lower = more diverse but smaller effective HoF.

    diversity_threshold : float
        If mean pairwise Tanimoto distance of the pool falls below this,
        explore_fraction is boosted.

    plateau_patience : int
        After this many generations without improvement, explore_fraction is
        shifted toward mutation to exploit the arm with the highest hit rate.

    Returns
    -------
    list of (smiles, score) tuples
        All HallOfFame members sorted best-first.
    """
    # ── resolve preset ────────────────────────────────────────────────────────
    p = _PRESETS.get(preset, _PRESETS["ml"]).copy()

    def _resolve(val, key):
        """Return val if not None, else preset default."""
        return val if val is not None else p[key]

    startN             = _resolve(startN,             "startN")
    popsize            = _resolve(popsize,             "popsize")
    scoring_budget     = _resolve(scoring_budget,      "scoring_budget")
    generations        = _resolve(generations,         "generations")
    plateau_patience   = _resolve(plateau_patience,    "plateau_patience")
    n_replacements     = _resolve(n_replacements,      "n_replacements")
    conservative       = _resolve(conservative,        "conservative")
    hof_size           = _resolve(hof_size,            "hof_size")
    pool_diversity_k   = _resolve(pool_diversity_k,    "pool_diversity_k")
    # maxmin_cap_factor: preset may define it; kwarg overrides; fallback=10 (original default)
    _mmcf_default = p.get("maxmin_cap_factor", 10)
    _maxmin_cap_factor = maxmin_cap_factor if maxmin_cap_factor is not None else _mmcf_default

    # Handle legacy aliases
    _base_ef = base_explore_fraction if base_explore_fraction is not None \
               else (explore_ratio if explore_ratio is not None
                     else p["base_explore_fraction"])
    _hof_sim = hof_sim_threshold if hof_sim_threshold is not None \
               else (sim_threshold if sim_threshold is not None
                     else p["hof_sim_threshold"])

    if fragcorpus is None:
        fragcorpus = load_corpus()
    random.seed(seed)
    sign = _sign(higher_is_better)

    if not quiet:
        print(f"[LACAN GA] preset={preset!r}  startN={startN}  popsize={popsize}"
              f"  scoring_budget={scoring_budget}/gen  generations={generations}"
              f"  base_explore_fraction={_base_ef:.2f}  maxmin_cap_factor={_maxmin_cap_factor}")
        if win_threshold is not None:
            print(f"[LACAN GA] Note: win_threshold={win_threshold} is ignored in this version. "
                  f"Results are returned as all HallOfFame entries (top {hof_size} diverse molecules).")

    # ── data structures ───────────────────────────────────────────────────────
    hof = HallOfFame(max_size=hof_size, sim_threshold=_hof_sim,
                     higher_is_better=higher_is_better)
    pool = []           # list of (smiles, internal_score)
    seen_smis = set()   # all SMILES ever scored

    # Bandit over 7 arms — Thompson Sampling with optional persistent stats
    EXPLORE_ARMS = ["ring", "linker", "sub", "decorate", "crossover", "random"]
    EXPLOIT_ARMS = ["mutate"]
    bandit = ThompsonBandit(EXPLORE_ARMS + EXPLOIT_ARMS,
                            stats_file=bandit_stats_file)

    # ── initial population ────────────────────────────────────────────────────
    if seed_mols is not None:
        if not quiet:
            print(f"  Seeding from {len(seed_mols)} provided molecules...")
        start_mols = [m for m in seed_mols if m is not None]
        for m in start_mols:
            try:
                seen_smis.add(Chem.MolToSmiles(m))
            except Exception:
                pass
    else:
        if not quiet:
            print(f"  Generating {startN} initial molecules...")
        start_mols = generate_filtered_molecules(profile, fragcorpus=fragcorpus,
                                                  min_atoms=15, n_molecules=startN,
                                                  n_jobs=n_jobs, seed=seed)
        for m in start_mols:
            try:
                seen_smis.add(Chem.MolToSmiles(m))
            except Exception:
                pass

    raw_scores = _safe_score(scoring_function, start_mols)
    for mol, rs in zip(start_mols, raw_scores):
        try:
            smi = Chem.MolToSmiles(mol)
        except Exception:
            continue
        if hof._worth_offering(rs):
            hof.offer(smi, rs, mol)
        pool.append((smi, rs * sign))

    # For the initial pool, keep a larger slice proportional to startN so that
    # a large initial sweep (e.g. startN=1000) benefits the first few GA
    # generations rather than being immediately discarded.  We use diversity-
    # aware selection so the larger pool is also well-spread, then decay to
    # popsize at the start of gen 1 via the normal per-generation cull.
    initial_pool_size = min(len(pool), max(popsize * pool_diversity_k,
                                           startN // 10))
    pool.sort(key=lambda x: x[1])
    pool = pool[:initial_pool_size * pool_diversity_k]
    pool = _diverse_top_k(pool, initial_pool_size, diversity_weight=0.3)

    plateau_count = 0
    explore_fraction = _base_ef
    # Initialise best_internal from the culled pool so the plateau tracker
    # has a consistent reference point from the start.
    best_internal = pool[0][1] if pool else 0.0

    if not quiet:
        print(f"  Initial pool: {len(pool)} molecules  |  HoF: {len(hof)}"
              f"  |  best={hof.best_score}")

    # ── generational loop ─────────────────────────────────────────────────────
    import time as _time
    _gen_t0 = _time.perf_counter()
    _viable_for_crossover = None  # cache: recomputed only when pool changes significantly

    for gen_idx in range(generations):
        _t_gen_start = _time.perf_counter()

        if seed:
            seed = (seed * 6543 + gen_idx) % (2**31)

        # ── phase: pool cull + diversity-select ───────────────────────────────
        _t0 = _time.perf_counter()
        pool.sort(key=lambda x: x[1])
        pool = pool[:popsize * pool_diversity_k]
        pool = _diverse_top_k(pool, popsize, diversity_weight=0.3)
        _t_pool_cull = _time.perf_counter() - _t0

        culled_best = pool[0][1] if pool else best_internal
        if culled_best < best_internal - 1e-6:
            best_internal = culled_best
            plateau_count = 0
        else:
            plateau_count += 1

        # ── phase: mean diversity ─────────────────────────────────────────────
        _t0 = _time.perf_counter()
        diversity = _mean_diversity([smi for smi, _ in pool])
        _t_diversity = _time.perf_counter() - _t0

        force_explore = plateau_count >= plateau_patience
        low_diversity = diversity < diversity_threshold

        if low_diversity and not force_explore:
            explore_fraction = min(explore_fraction + 0.10, 0.60)
        elif force_explore:
            explore_fraction = max(explore_fraction - 0.05, 0.05)
        else:
            explore_fraction = max(explore_fraction - 0.03, _base_ef)

        n_exploit_budget = max(1, int(scoring_budget * (1 - explore_fraction)))
        n_explore_budget = max(1, scoring_budget - n_exploit_budget)

        parent_mols = [Chem.MolFromSmiles(smi) for smi, _ in pool]
        parent_score_map = {smi: sc for smi, sc in pool}
        arm_candidates = {a: [] for a in EXPLORE_ARMS + EXPLOIT_ARMS}

        # ── phase: fragment moves (ring/linker/sub/decorate) ──────────────────
        _t0 = _time.perf_counter()
        parent_results = _fragment_moves_all_parents(
            pool, profile, fragcorpus,
            n_replacements=n_replacements,
            conservative=conservative,
            n_jobs=n_jobs,
        )
        for (parent_sc, by_type) in parent_results:
            for arm in ["ring", "linker", "sub", "decorate"]:
                for m in by_type.get(arm, []):
                    arm_candidates[arm].append((m, parent_sc))
        _t_frag = _time.perf_counter() - _t0

        # ── Budget allocation via Thompson Sampling ───────────────────────────
        ts_weights = bandit.sample_weights()

        def _allocate_budget(arms, total_budget):
            weights = {a: ts_weights[a] for a in arms}
            total_w = sum(weights.values())
            per_arm = {}
            remaining = total_budget
            arm_list = list(arms)
            for i, a in enumerate(arm_list):
                if i == len(arm_list) - 1:
                    per_arm[a] = max(1, remaining)
                else:
                    if total_w > 0:
                        alloc = max(1, int(total_budget * weights[a] / total_w))
                    else:
                        alloc = max(1, total_budget // len(arm_list))
                    per_arm[a] = alloc
                    remaining = max(0, remaining - alloc)
            return per_arm

        explore_alloc = _allocate_budget(EXPLORE_ARMS, n_explore_budget)
        exploit_alloc = _allocate_budget(EXPLOIT_ARMS, n_exploit_budget)
        arm_alloc = {**explore_alloc, **exploit_alloc}

        # ── phase: crossover ──────────────────────────────────────────────────
        # We call breed.breed directly (not cross_breed_mols) to control
        # hacrange. Viability check (fragment_molecule) is cached across
        # generations — recomputed only when the pool population changes
        # (i.e. when plateau_count resets, meaning new best molecules entered).
        _t0 = _time.perf_counter()
        _crossover_err = None
        pool_mean_sc = sum(sc for _, sc in pool) / len(pool) if pool else 0.0
        if len(parent_mols) >= 2:
            try:
                # Recompute viable parents only on improvement (plateau reset)
                # or if cache is empty. Saves ~2-4s of fragment_molecule calls.
                if _viable_for_crossover is None or plateau_count == 0:
                    _viable_for_crossover = []
                    for m in parent_mols:
                        if m is None:
                            continue
                        try:
                            s, c = breed.fragment_molecule(m, 3)
                            if not c:
                                s, c = breed.fragment_molecule(m, 2)
                            if s and c:
                                _viable_for_crossover.append(m)
                        except Exception:
                            pass

                if len(_viable_for_crossover) >= 2:
                    _nmols_per_pair = max(1, explore_alloc["crossover"] // max(1, len(_viable_for_crossover)))
                    _actual_jobs = multiprocessing.cpu_count() if n_jobs == -1 else max(1, n_jobs)
                    _breed_inputs = [
                        (m, random.choice(_viable_for_crossover), profile,
                         _nmols_per_pair, 3, (0.6, 1.5), (0.25, 0.75), False)
                        for m in _viable_for_crossover
                    ]
                    with multiprocessing.Pool(processes=_actual_jobs) as _xpool:
                        _xresults = _xpool.starmap(breed.breed, _breed_inputs)
                    for _xmols in _xresults:
                        for m in _xmols:
                            arm_candidates["crossover"].append((m, pool_mean_sc))
                else:
                    _crossover_err = f"too few fragmentable parents ({len(_viable_for_crossover) if _viable_for_crossover else 0}/{len(parent_mols)})"
            except Exception as _e:
                _crossover_err = repr(_e)
        _t_crossover = _time.perf_counter() - _t0

        # ── phase: random injection ───────────────────────────────────────────
        _t0 = _time.perf_counter()
        n_rand = max(2, explore_alloc["random"])
        try:
            randoms = generate_filtered_molecules(
                profile, fragcorpus=fragcorpus, min_atoms=14,
                n_molecules=n_rand, n_jobs=n_jobs, seed=seed)
            for m in randoms:
                arm_candidates["random"].append((m, pool_mean_sc))
        except Exception:
            pass
        _t_random = _time.perf_counter() - _t0

        # ── phase: mutation ───────────────────────────────────────────────────
        _t0 = _time.perf_counter()
        try:
            mutants = mutate.apply_mutations_mols(
                parent_mols, profile, score_threshold=0.0, n_jobs=n_jobs)
            # Cap raw mutants: scoring budget for mutate arm * small multiplier.
            # apply_mutations_mols on 200 parents produces ~18k candidates but
            # we only need ~425. A cap of 5x leaves enough diversity without
            # running 18k through the property filter loop.
            _mutant_cap = exploit_alloc["mutate"] * 5
            if len(mutants) > _mutant_cap:
                random.shuffle(mutants)
                mutants = mutants[:_mutant_cap]
            for m in mutants:
                arm_candidates["mutate"].append((m, pool_mean_sc))
        except Exception as _e:
            if not quiet:
                print(f"    MUTATE_ERR: {repr(_e)}")
        _t_mutate = _time.perf_counter() - _t0

        # ── phase: prescreening ───────────────────────────────────────────────
        _t0 = _time.perf_counter()
        prescreened = {}
        _prescreen_counts = {}  # raw -> screened per arm, for diagnostics
        for arm, alloc in arm_alloc.items():
            raw_pairs = arm_candidates[arm]
            if not raw_pairs:
                prescreened[arm] = []
                bandit.update(arm, False)
                _prescreen_counts[arm] = (0, 0)
                continue
            raw_mols = [m for m, _ in raw_pairs]
            raw_smi_to_parent_sc = {}
            for m, psc in raw_pairs:
                try:
                    raw_smi_to_parent_sc[Chem.MolToSmiles(m)] = psc
                except Exception:
                    pass
            screened = _prescreened_candidates(raw_mols, alloc, seen_smis,
                                                maxmin_cap_factor=_maxmin_cap_factor)
            prescreened[arm] = [(smi, mol, raw_smi_to_parent_sc.get(smi, pool_mean_sc))
                                 for smi, mol in screened]
            _prescreen_counts[arm] = (len(raw_pairs), len(screened))
        _t_prescreen = _time.perf_counter() - _t0

        # ── phase: scoring ────────────────────────────────────────────────────
        all_candidates = []
        for arm in EXPLORE_ARMS + EXPLOIT_ARMS:
            for item in prescreened[arm]:
                all_candidates.append((arm, item[0], item[1], item[2]))

        if not all_candidates:
            plateau_count += 1
            if not quiet:
                print(f"  Gen {gen_idx+1}/{generations}  [no candidates — skipping]")
            continue

        _t0 = _time.perf_counter()
        score_mols = [mol for _, _, mol, _ in all_candidates]
        raw_scores = _safe_score(scoring_function, score_mols)
        _t_score = _time.perf_counter() - _t0

        arm_yields = {a: [] for a in EXPLORE_ARMS + EXPLOIT_ARMS}
        for (arm, smi, mol, parent_sc), rs in zip(all_candidates, raw_scores):
            internal = rs * sign
            pool.append((smi, internal))
            if hof._worth_offering(rs):
                hof.offer(smi, rs, mol)
            arm_yields[arm].append(max(0.0, parent_sc - internal))

        # ── immediate pool cull ───────────────────────────────────────────────
        # Cull right after scoring so the pool never grows to 20k+ entries.
        # Without this, _diverse_top_k at the next gen start gets a huge list.
        _pool_hard_cap = popsize * pool_diversity_k * 4  # generous but bounded
        if len(pool) > _pool_hard_cap:
            pool.sort(key=lambda x: x[1])
            pool = pool[:_pool_hard_cap]

        # ── phase: HoF offers (counted above, separate timing) ────────────────
        for arm, yields in arm_yields.items():
            success = any(y > 1e-6 for y in yields) if yields else False
            bandit.update(arm, success)
        bandit.save()

        _t_gen_total = _time.perf_counter() - _t_gen_start

        # ── per-generation report ─────────────────────────────────────────────
        if not quiet:
            best_raw = pool[0][1] / sign if pool else float("nan")

            # Candidate counts per arm: "mutate:2134->423"
            arm_counts = "  ".join(
                f"{a}:{_prescreen_counts.get(a,(0,0))[0]}->{_prescreen_counts.get(a,(0,0))[1]}"
                for a in EXPLORE_ARMS + EXPLOIT_ARMS
            )

            print(
                f"  Gen {gen_idx+1}/{generations}"
                f"  pool={len(pool)}"
                f"  explore_f={explore_fraction:.2f}"
                f"  diversity={diversity:.3f}"
                f"  plateau={plateau_count}"
                f"  best_pool={best_raw:.4f}"
                f"  HoF={len(hof)}  best_hof={hof.best_score}"
                f"\n    TIMING(s): total={_t_gen_total:.1f}"
                f"  pool_cull={_t_pool_cull:.2f}"
                f"  diversity={_t_diversity:.2f}"
                f"  frag={_t_frag:.2f}"
                f"  crossover={_t_crossover:.2f}"
                f"  random={_t_random:.2f}"
                f"  mutate={_t_mutate:.2f}"
                f"  prescreen={_t_prescreen:.2f}"
                f"  score={_t_score:.2f}"
                f"  scored={len(all_candidates)}"
                f"\n    CANDIDATES: {arm_counts}"
                + (f"\n    CROSSOVER_ERR: {_crossover_err}" if _crossover_err else "")
            )

        if callback is not None:
            best_pool_raw = (min(s for _, s in pool) / sign) if pool else None
            callback({
                "generation":       gen_idx + 1,
                "explore_fraction": explore_fraction,
                "pool_size":        len(pool),
                "hof_size":         len(hof),
                "diversity":        diversity,
                "plateau":          plateau_count,
                "best_pool":        best_pool_raw,
                "best_hof":         hof.best_score,
                "bandit":           bandit.hit_rate,
                "timing": {
                    "total":      _t_gen_total,
                    "pool_cull":  _t_pool_cull,
                    "diversity":  _t_diversity,
                    "frag":       _t_frag,
                    "crossover":  _t_crossover,
                    "random":     _t_random,
                    "mutate":     _t_mutate,
                    "prescreen":  _t_prescreen,
                    "score":      _t_score,
                    "n_scored":   len(all_candidates),
                },
            })

    if not quiet:
        print(f"\n[LACAN GA] Done. HallOfFame: {len(hof)} molecules.")
        print("[LACAN GA] Per-arm bandit summary (posterior hit rates):")
        print(bandit.summary())
        if bandit_stats_file:
            print(f"[LACAN GA] Bandit stats saved to {bandit_stats_file}")

    return hof.top()


class GAReporter:
    """Collects per-generation statistics from :func:`generate_optimized_molecules` and plots them.

    Pass an instance as the ``callback=`` argument to the GA.  After the run,
    call :meth:`plot` to visualise the results, or compare multiple runs by
    passing a list of reporters to :func:`compare`.

    Parameters
    ----------
    label : str
        Name for this run, shown in plot legends (default ``"run"``).

    Example
    -------
    ::

        reporter = GAReporter(label="docking_preset")
        winners = gen.generate_optimized_molecules(
            my_score_fn, profile,
            preset="docking",
            callback=reporter,
        )
        reporter.plot()
    """

    def __init__(self, label="run"):
        self.label = label
        self.history = []

    def __call__(self, stats):
        self.history.append(stats)

    def plot(self, ax=None, show=True):
        """Plot score, diversity, explore fraction, and per-arm hit rates over time."""
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        if not self.history:
            print("No data recorded yet.")
            return None

        gens         = [h["generation"]       for h in self.history]
        best_pool    = [h["best_pool"]         for h in self.history]
        best_hof     = [h["best_hof"]          for h in self.history]
        diversity    = [h["diversity"]         for h in self.history]
        explore_f    = [h["explore_fraction"]  for h in self.history]
        hof_sizes    = [h["hof_size"]          for h in self.history]

        fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
        fig.suptitle(f"GA run: {self.label}", fontsize=13, fontweight="bold")

        # Panel 1: scores
        ax1 = axes[0]
        ax1.plot(gens, best_pool, "o-", color="steelblue", label="best in pool", lw=1.5)
        ax1.plot(gens, best_hof,  "s--", color="darkorange", label="best HoF", lw=1.5)
        ax1.set_ylabel("Score")
        ax1.legend(loc="lower right", fontsize=9)
        ax1.set_title("Score per generation")

        # Panel 2: diversity + explore fraction + HoF size
        ax2 = axes[1]
        ax2.plot(gens, diversity, "^-", color="mediumpurple", label="diversity", lw=1.5)
        ax2.plot(gens, explore_f, "v--", color="teal", label="explore_fraction", lw=1.5)
        ax2_r = ax2.twinx()
        ax2_r.bar(gens, hof_sizes, color="salmon", alpha=0.4, label="HoF size")
        ax2_r.set_ylabel("HoF size", color="salmon")
        ax2_r.tick_params(axis="y", labelcolor="salmon")
        ax2.set_ylabel("Fraction / Distance")
        ax2.set_ylim(0, 1)
        ax2.legend(loc="upper left", fontsize=9)
        ax2_r.legend(loc="upper right", fontsize=9)
        ax2.set_title("Diversity, explore fraction, and HallOfFame size")

        # Panel 3: bandit arm hit rates over time
        ax3 = axes[2]
        colors3 = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        arm_names = list(self.history[0]["bandit"].keys()) if self.history else []
        for i, arm in enumerate(arm_names):
            rates = [h["bandit"].get(arm, 0) for h in self.history]
            ax3.plot(gens, rates, "o-", label=arm,
                     color=colors3[i % len(colors3)], lw=1.5, alpha=0.85)
        ax3.set_ylabel("Posterior hit rate")
        ax3.set_xlabel("Generation")
        ax3.set_ylim(0, None)
        ax3.set_title("Bandit arm posterior hit rates (Thompson Sampling)")
        ax3.legend(fontsize=8, loc="upper right", ncol=4)

        plt.tight_layout()
        if show:
            plt.show()
        return fig

    @staticmethod
    def compare(reporters, metric="best_pool", show=True):
        """Overlay score curves from multiple runs."""
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(9, 4))
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

        for i, rep in enumerate(reporters):
            if not rep.history:
                continue
            gens = [h["generation"] for h in rep.history]
            vals = [h[metric] for h in rep.history]
            if metric == "best_hof":
                vals = [v if v is not None else float("nan") for v in vals]
            ax.plot(gens, vals, "o-", color=colors[i % len(colors)],
                    label=rep.label, lw=1.8)

        ylabel = {
            "best_pool": "Score (best in pool)",
            "best_hof":  "Score (best HoF)",
            "diversity": "Mean Tanimoto distance",
        }.get(metric, metric)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Generation")
        ax.set_title(f"GA comparison — {metric}")
        ax.legend(fontsize=9)
        plt.tight_layout()
        if show:
            plt.show()
        return fig

    def summary(self):
        """Print a compact summary table of the run."""
        if not self.history:
            print("No data.")
            return
        print(f"\n{'Gen':>4}  {'ExpF':>5}  {'Pool':>5}  {'HoF':>5}  "
              f"{'Div':>6}  {'BestPool':>9}  {'BestHoF':>9}")
        print("-" * 58)
        for h in self.history:
            bh = f"{h['best_hof']:.4f}" if h["best_hof"] is not None else "    —   "
            bp = f"{h['best_pool']:.4f}" if h["best_pool"] is not None else "    —   "
            print(f"{h['generation']:>4}  {h['explore_fraction']:>5.2f}"
                  f"  {h['pool_size']:>5}  {h['hof_size']:>5}  "
                  f"{h['diversity']:>6.3f}  {bp:>9}  {bh:>9}")
        print()


def optimize_from_mol(
        mol,
        scoring_function,
        profile,
        generations=20,
        beam_width=10,
        n_replacements=15,
        win_threshold=0.8,
        higher_is_better=True,
        conservative=True,
        fragcorpus=None,
        plateau_patience=5,
        quiet=False,
        callback=None,
        protect_smarts=None):
    """Optimise a single seed molecule using fragment-level operations.

    Unlike :func:`generate_optimized_molecules`, which explores chemical space
    from random starting structures, this function starts from *one specific
    molecule* and iteratively improves it by applying ring, linker, substituent,
    and decoration operations — then retaining the best results as parents for
    the next round (beam search).

    Parameters
    ----------
    mol              : RDKit Mol — the seed / lead molecule
    scoring_function : callable — ``fn(list[Mol]) -> list[float]``
    profile          : LACAN profile dict
    generations      : int — maximum number of rounds (default 20)
    beam_width       : int — top molecules kept as parents each round (default 10)
    n_replacements   : int — fragment operation attempts per parent (default 15)
    win_threshold    : float — score at which a molecule is reported as winner
                       (default 0.8).  If no winners found, best beam returned.
    higher_is_better : bool (default True)
    conservative     : bool — prefer similar fragment replacements (default True)
    fragcorpus       : fragment corpus (default: ChEMBL rls.csv)
    plateau_patience : int — stop early after this many generations without
                       improvement (default 5)
    quiet            : bool (default False)
    callback         : callable or None — called each generation with stats dict
    protect_smarts   : str or None — SMARTS pattern; atoms matching it are skipped
                       by all fragment operations every generation.  The protected
                       set is re-derived from each parent mol on every call via
                       :func:`~lacan.protect.get_protected_atoms`, so it survives
                       the SMILES round-trip in the beam transparently.  ``None``
                       disables atom exclusion entirely (default).

    Returns
    -------
    list of (smiles, score) tuples, sorted best-first
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()

    sign = _sign(higher_is_better)
    _win = win_threshold * sign

    seed_smi = Chem.MolToSmiles(mol)
    seed_score_raw = _safe_score(scoring_function, [mol])[0]
    beam = [(seed_smi, seed_score_raw * sign)]
    winners = []
    all_seen_smis = {seed_smi}

    if seed_score_raw * sign < _win:
        winners.append((seed_smi, seed_score_raw))

    best_internal = beam[0][1]
    plateau_count = 0

    if not quiet:
        print(f"optimize_from_mol: seed score = {seed_score_raw:.4f}")

    for gen_idx in range(generations):
        parent_mols = [Chem.MolFromSmiles(smi) for smi, _ in beam]
        candidates = []
        for parent in parent_mols:
            candidates += _fragment_moves(parent, profile, fragcorpus,
                                          n_replacements=n_replacements,
                                          conservative=conservative,
                                          protect_smarts=protect_smarts)

        novel = []
        for m in candidates:
            try:
                smi = Chem.MolToSmiles(m)
            except Exception:
                continue
            if smi not in all_seen_smis:
                all_seen_smis.add(smi)
                novel.append((smi, m))

        if not novel:
            plateau_count += 1
            if not quiet:
                print(f"  Gen {gen_idx+1}: no novel candidates — plateau {plateau_count}/{plateau_patience}")
            if plateau_count >= plateau_patience:
                break
            continue

        novel_smis, novel_mols = zip(*novel)
        raw_scores = _safe_score(scoring_function, list(novel_mols))

        new_beam_entries = list(beam)
        gen_best = best_internal

        for smi, rs in zip(novel_smis, raw_scores):
            internal = rs * sign
            if internal < _win:
                winners.append((smi, rs))
            new_beam_entries.append((smi, internal))
            if internal < gen_best:
                gen_best = internal

        new_beam_entries = sorted(new_beam_entries, key=lambda x: x[1])[:beam_width]
        beam = new_beam_entries

        if gen_best < best_internal - 1e-6:
            best_internal = gen_best
            plateau_count = 0
        else:
            plateau_count += 1

        best_raw = best_internal / sign
        if not quiet:
            print(f"  Gen {gen_idx+1}/{generations} | candidates={len(novel)} "
                  f"| winners={len(winners)} | best={best_raw:.4f} "
                  f"| plateau={plateau_count}/{plateau_patience}")

        if callback is not None:
            callback({
                "generation":   gen_idx + 1,
                "n_candidates": len(novel),
                "n_winners":    len(winners),
                "best_score":   best_raw,
                "plateau":      plateau_count,
            })

        if plateau_count >= plateau_patience:
            if not quiet:
                print("  Plateau reached — stopping early.")
            break

    if not winners:
        if not quiet:
            print("  No winners above threshold — returning best beam members.")
        winners = [(smi, sc / sign) for smi, sc in beam]

    seen, deduped = set(), []
    for smi, sc in sorted(winners, key=lambda x: x[1], reverse=higher_is_better):
        if smi not in seen:
            seen.add(smi)
            deduped.append((smi, sc))
    return deduped


def next_population(population):
    """Placeholder for a future incremental population-advance API."""
    raise NotImplementedError("next_population is not yet implemented")
