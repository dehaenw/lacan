"""
gen.py — Random molecule generation and adaptive genetic algorithm for LACAN.

This module has two responsibilities:

1. **Random generation** — assemble drug-like molecules from scratch by
   recursively sampling fragments from the corpus and joining them through
   their dummy-atom attachment points.

2. **Optimisation** — run an adaptive GA (:func:`generate_optimized_molecules`)
   that drives any external scoring function toward molecules that score above
   a user-defined threshold while maintaining population diversity.

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
"""

from rdkit import Chem
from lacan import lacan
from lacan import mutate, breed, decompose
from lacan import replace as replace_module
import random
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


def load_corpus(min_count=200):
    """Load the ChEMBL fragment corpus from ``data/rls.csv``, with caching.

    The CSV is read only on the first call for a given ``min_count`` value;
    subsequent calls return the cached result immediately.  This makes the
    function safe to use as a default-argument resolver inside other functions
    without causing module-level I/O at import time (which would break Sphinx
    autodoc and slow down any code that merely imports the module).

    Parameters
    ----------
    min_count : int
        Minimum occurrence count for a fragment to be included (default 200).
        Fragments appearing fewer than this many times in ChEMBL are excluded
        to avoid noise.  Lower values include rarer fragments; higher values
        restrict to the most common ones.

    Returns
    -------
    list of [smiles, count, degree, ftype, bonds]
        Corpus entries sorted by count descending.
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
    """Return total occurrence counts per fragment type for *corpus*.

    Parameters
    ----------
    corpus : list of corpus entries

    Returns
    -------
    dict with keys ``"Sub"``, ``"Linker"``, ``"Ring"``
    """
    fc = {"Sub": 0, "Linker": 0, "Ring": 0}
    for e in corpus:
        fc[e[3]] += e[1]
    return fc


def generate_random_molecule(fragments=[], fragcorpus=None, needed_sample=[], bondstring="", mol=None):
    """Recursively assemble a random molecule from corpus fragments.

    This is a recursive builder.  On the first call (``fragments=[]``) it picks
    a seed fragment from the corpus, then calls itself to fill each of the
    seed's open attachment points with appropriately typed fragments.

    **Attachment-point grammar** — every corpus fragment carries a ``bonds``
    string (e.g. ``"--"`` or ``"-="``).  Each character encodes the bond type
    expected at one dummy atom.  When a Ring is placed, its remaining dummies
    are filled with Subs or Linkers; when a Linker is placed, its remaining
    dummy must connect to a Ring; when a Sub is placed it has no further
    dummies and the branch terminates.

    Sampling weights are ``count^0.666`` (slightly more uniform than
    frequency-proportional to avoid the most common fragments dominating
    completely).

    Parameters
    ----------
    fragments     : list — accumulates SMILES of all fragments used so far
                    (internal state; pass ``[]`` or omit on first call)
    fragcorpus    : list of corpus entries (default: ChEMBL rls.csv)
    needed_sample : list of lists — each entry is a list of allowed fragment
                    types for the next open attachment point
                    (internal state; omit on first call)
    bondstring    : str — bond types for the current open attachment points
                    (internal state; omit on first call)
    mol           : RDKit Mol — molecule built so far
                    (internal state; omit on first call)

    Returns
    -------
    (mol, fragments) : RDKit Mol and list of fragment SMILES used
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

    Repeatedly calls :func:`generate_random_molecule` until it produces a
    molecule whose LACAN score ≥ ``threshold`` and whose heavy-atom count is
    in ``(min_atoms, max_atoms)``.  There is no iteration cap — this will loop
    indefinitely if the constraints are too tight for the corpus.

    Parameters
    ----------
    profile   : LACAN profile dict
    fragcorpus: fragment corpus (default: ChEMBL rls.csv)
    threshold : minimum LACAN score to accept (default 0.5)
    seed      : int or None — if provided, seeds the random state before the
                first attempt (useful for reproducible parallel generation)
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

    When ``n_jobs=1`` molecules are generated sequentially.  For ``n_jobs=-1``
    or any value > 1 a ``multiprocessing.Pool`` is used; each worker receives a
    unique seed (``seed + 823848 * i``) so they explore different regions of
    chemical space independently.

    Parameters
    ----------
    profile     : LACAN profile dict
    fragcorpus  : fragment corpus (default: ChEMBL rls.csv)
    threshold   : minimum LACAN score (default 0.5)
    seed        : int or None — base seed; each worker's seed is derived from
                  this by an offset
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
                            [(profile, fragcorpus, threshold, seed + 823848 * i, min_atoms, max_atoms)
                             for i in range(n_molecules)])
        pool.close()
        pool.join()
    return mols


def get_corpus_from_mols(mols, fragcorpus=None, ratio=1):
    """Merge a custom fragment corpus from reference molecules into the background corpus.

    Decomposes each molecule in *mols* via :func:`~lacan.decompose.get_corpus`,
    then merges the resulting fragment counts into *fragcorpus* with a
    multiplier of *ratio*.  Fragments present in both corpora have their counts
    added; fragments only in the custom set are inserted at their boosted count;
    fragments only in the background are kept unchanged.

    The net effect is that sampling from the returned corpus will draw reference
    fragments with approximately ``ratio`` times higher probability than their
    background frequency alone would give.

    This is the low-level implementation; the public interface is
    :func:`bias_corpus`.

    Parameters
    ----------
    mols       : iterable of RDKit Mol objects — reference molecules
    fragcorpus : background corpus to merge into (default: ChEMBL rls.csv)
    ratio      : float — frequency multiplier for reference fragments (default 1,
                 which gives equal weighting; 2 doubles the effective frequency)

    Returns
    -------
    list of corpus entries (same format as *fragcorpus*)
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()
    custom_entries = decompose.get_corpus(mols)  # this can take a while ...
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


def _fragment_moves(mol, profile, fragcorpus, n_replacements=5, conservative=True, min_atoms=5):
    """Apply all fragment-level operations to a single molecule.

    This is the per-molecule workhorse shared by :func:`_explore` and
    :func:`optimize_from_mol`.  It runs :func:`~lacan.replace.replace_ring`,
    :func:`~lacan.replace.replace_substituent`,
    :func:`~lacan.replace.replace_linker`, and
    :func:`~lacan.replace.decorate_scaffold` on *mol* and returns the union of
    their outputs.

    All operations use ``score_threshold=0.0`` so that partially-valid
    intermediates are not discarded here — the caller applies any
    objective-function filter.

    Parameters
    ----------
    mol            : RDKit Mol
    profile        : LACAN profile dict
    fragcorpus     : fragment corpus list
    n_replacements : attempts per operation (default 5)
    conservative   : if True, prefer structurally similar replacements
    min_atoms      : discard products with fewer heavy atoms (default 5)

    Returns
    -------
    list of RDKit Mol — deduplicated by SMILES within this call.
    """
    new_mols = []
    for op in [replace_module.replace_ring,
               replace_module.replace_substituent,
               replace_module.replace_linker]:
        try:
            new_mols += op(mol, profile, score_threshold=0.0,
                           fragcorpus=fragcorpus, n_replacements=n_replacements,
                           conservative=conservative)
        except Exception:
            pass
    try:
        new_mols += replace_module.decorate_scaffold(
            mol, profile, score_threshold=0.0,
            fragcorpus=fragcorpus, n_replacements=n_replacements,
            mode="Hydrogen")
    except Exception:
        pass
    seen, out = set(), []
    for m in new_mols:
        if m is None or m.GetNumAtoms() < min_atoms:
            continue
        smi = Chem.MolToSmiles(m)
        if smi not in seen:
            seen.add(smi)
            out.append(m)
    return out


def _explore(mols, profile, fragcorpus, n_random, seed, min_atoms=5, conservative=True, n_jobs=-1):
    """
    Exploration step: coarse fragment operations + crossover + fresh random mols.
    Returns a flat list of new candidate molecules (filtered to min_atoms).
    Each operation is wrapped in try/except so a bad molecule never kills the step.
    """
    new_mols = []
    for mol in mols:
        new_mols += _fragment_moves(mol, profile, fragcorpus,
                                    n_replacements=5, conservative=conservative,
                                    min_atoms=min_atoms)
    # crossover between population members
    if len(mols) >= 2:
        try:
            new_mols += breed.cross_breed_mols(mols, profile, score_threshold=0.0, nmols=2, n_jobs=n_jobs)
        except Exception:
            pass
    # fresh random molecules to inject diversity
    try:
        new_mols += generate_filtered_molecules(profile, fragcorpus=fragcorpus,
                                                min_atoms=14, n_molecules=n_random, n_jobs=n_jobs, seed=seed)
    except Exception:
        pass
    # Filter out very small molecules that can crash 3D-based scoring functions
    return [m for m in new_mols if m.GetNumAtoms() >= min_atoms]


def _exploit(mols, profile, min_atoms=5, n_jobs=-1):
    """
    Exploitation step: fine-grained atom-level mutations.
    Returns a flat list of new candidate molecules (filtered to min_atoms).
    """
    try:
        result = mutate.apply_mutations_mols(mols, profile, score_threshold=0.0, n_jobs=n_jobs)
        return [m for m in result if m.GetNumAtoms() >= min_atoms]
    except Exception:
        return []


def _safe_score(scoring_function, mols):
    """
    Call scoring_function and return a list of scores, replacing any
    exception (e.g. from 3D embedding failures in shape-align) with 0.0.
    Also filters out None mols before calling the function and fills
    back 0.0 for those positions.
    """
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
    # Map back, filling 0.0 for None mols
    result = [0.0] * len(mols)
    for idx, score in zip(valid_indices, raw):
        try:
            result[idx] = float(score)
        except Exception:
            result[idx] = 0.0
    return result


def _sign(higher_is_better):
    """Return -1 if higher scores are better, +1 if lower scores are better.

    The GA internally minimises all scores.  Multiplying raw scores by this
    value converts maximisation objectives to minimisation objectives so the
    rest of the GA logic is uniform.
    """
    return -1 if higher_is_better else 1


def bias_corpus(mols, fragcorpus=None, ratio=2.0):
    """
    Build a fragment corpus biased toward the chemistry of the provided molecules.

    Fragments found in *mols* have their occurrence counts boosted relative to
    the background corpus. *ratio* controls how strongly: at ratio=1 the custom
    fragments are weighted equally to their frequency in the background; at
    ratio=2 (default) they are weighted twice as heavily, so the GA will
    preferentially sample fragments that appear in your reference molecules.

    Typical use: pass a set of known actives, scaffold set, or reference
    structures to steer the GA toward chemotypes you care about.

    Returns a merged corpus list that can be passed as fragcorpus= to
    generate_filtered_molecules() or generate_optimized_molecules().

    Example
    -------
        actives = [Chem.MolFromSmiles(s) for s in active_smiles]
        biased = gen.bias_corpus(actives, ratio=3.0)
        winners = gen.generate_optimized_molecules(score_fn, profile,
                                                   fragcorpus=biased, ...)
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()
    return get_corpus_from_mols(mols, fragcorpus=fragcorpus, ratio=ratio)


def generate_optimized_molecules(
        scoring_function,
        profile,
        seed=123,
        startN=50,
        generations=10,
        popsize=20,
        win_threshold=0.8,
        sim_threshold=0.45,
        higher_is_better=True,
        # Diversity threshold: below this mean pairwise distance, switch to explore
        diversity_threshold=0.4,
        # After this many generations without improvement, force an explore step
        plateau_patience=3,
        # Fraction of a generation budget spent on exploration vs exploitation
        # when the GA freely chooses (0 = all exploit, 1 = all explore)
        explore_ratio=0.5,
        quiet=False,
        conservative=True,
        fragcorpus=None,
        n_jobs=-1,
        callback=None,
        seed_mols=None):
    """
    GA with adaptive exploration / exploitation.

    Decision logic per generation
    ──────────────────────────────
    The GA tracks two signals:

      diversity  = mean pairwise Tanimoto distance of the current population
      plateau    = number of consecutive generations with no improvement

    From those it picks a mode:

      EXPLORE  if  diversity < diversity_threshold  OR  plateau >= plateau_patience
      EXPLOIT  otherwise

    EXPLORE uses coarse fragment operations (ring/linker/substituent replacement,
    scaffold decoration) plus crossover and fresh random molecules.
    These make big structural jumps and are good for finding new chemotypes.

    EXPLOIT uses fine-grained atom mutations. These refine good leads and are
    fast, so they suit both slow (docking) and fast (QSAR) objectives.

    For slow scoring functions, reduce startN/popsize and set n_jobs=-1.
    For fast ones, you can increase both freely.

    conservative : bool (default True)
        If True, fragment replacement operations sample structurally similar
        replacements with higher probability. Set to False for more random
        exploration (useful if the search is stuck in a chemical series).

    n_jobs : int (default -1)
        Number of parallel worker processes for mutation and generation steps.
        -1 uses all available CPU cores. Set to 1 to disable parallelism
        (useful for debugging or when the scoring function itself is already
        parallelised and spawning additional pools would be counterproductive).

    callback : callable or None (default None)
        If provided, called at the end of each generation with a dict of
        per-generation statistics.  Use :class:`GAReporter` for built-in
        plotting and comparison.  The dict contains:

        * ``generation``  — 1-based generation number
        * ``mode``        — ``"EXPLORE"`` or ``"EXPLOIT"``
        * ``pool_size``   — number of molecules in the active pool
        * ``n_winners``   — cumulative winner count
        * ``diversity``   — mean pairwise Tanimoto distance of the pool
        * ``plateau``     — consecutive generations without improvement
        * ``best_pool``   — best raw score in the current pool (None if empty)
        * ``best_winner`` — best raw score among all winners so far (None if none yet)

    seed_mols : list of RDKit Mol or None (default None)
        If provided, the initial population is seeded with these molecules
        (scored and split into winners/pool) instead of generating *startN*
        random structures from scratch.  Useful when you have a known active
        series or reference structures you want to optimise from.  The GA will
        still inject fresh random molecules in later explore steps if the pool
        runs dry.

    ┌─────────────────────────────────────────────────┐
    │  START: generate startN random molecules        │
    │         score all → split into winners/pool     │
    └──────────────────┬──────────────────────────────┘
                       │
              ┌────────▼────────┐
              │  for each gen   │◄────────────────────────┐
              └────────┬────────┘                         │
                       │                                  │
           ┌───────────▼────────────┐                     │
           │  diversity < threshold │                     │
           │    OR plateau hit?     │                     │
           └──────┬─────────┬───────┘                     │
                  │ YES     │ NO                           │
           EXPLORE▼         ▼EXPLOIT                      │
      ring/sub/linker    mutations                         │
      decoration         (fine-grained)                   │
      crossover                                           │
      +fresh randoms                                      │
                  │         │                             │
                  └────┬────┘                             │
                       │                                  │
              ┌────────▼────────────┐                     │
              │  score new mols     │                     │
              │  update winners     │                     │
              │  cull pool to       │                     │
              │  popsize by score   │                     │
              └────────┬────────────┘                     │
                       │                                  │
              ┌────────▼────────┐                         │
              │  more gens?     ├─ YES ───────────────────┘
              └────────┬────────┘
                       │ NO
              ┌────────▼────────┐
              │  return winners │
              └─────────────────┘
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()
    random.seed(seed)
    sign = _sign(higher_is_better)  # multiply scores so lower is always better internally
    _win = win_threshold * sign

    winners = []       # list of (smiles, raw_score) — molecules that beat win_threshold
    winnerfps = []     # fps of winners for diversity gating
    pool = []          # list of (smiles, internal_score) — active population

    # ── initial population ────────────────────────────────────────────────────
    if seed_mols is not None:
        if not quiet:
            print(f"Seeding initial population from {len(seed_mols)} provided molecules...")
        start_mols = [m for m in seed_mols if m is not None]
    else:
        if not quiet:
            print("Generating initial population...")
        start_mols = generate_filtered_molecules(profile, fragcorpus=fragcorpus,
                                                 min_atoms=15, n_molecules=startN,
                                                 n_jobs=n_jobs, seed=seed)
    raw_scores = _safe_score(scoring_function, start_mols)
    for mol, rs in zip(start_mols, raw_scores):
        smi = Chem.MolToSmiles(mol)
        internal = rs * sign
        if internal < _win:
            winners.append((smi, rs))
            winnerfps.append(MFPGEN.GetFingerprint(mol))
        else:
            pool.append((smi, internal))

    if not quiet:
        print(f"Initial population: {len(pool)} in pool, {len(winners)} winners")

    plateau_count = 0
    best_internal = min((s for _, s in pool), default=0.0)

    # ── generational loop ─────────────────────────────────────────────────────
    for gen in range(generations):
        if seed:
            seed = (seed * 6543 + gen) % (2**31)

        # Sort pool and trim to popsize
        pool = sorted(pool, key=lambda x: x[1])[:popsize]

        # Remove molecules too similar to existing winners
        if winnerfps:
            pool = [(smi, sc) for smi, sc in pool
                    if max(DataStructs.BulkTanimotoSimilarity(
                        MFPGEN.GetFingerprint(Chem.MolFromSmiles(smi)), winnerfps)) < sim_threshold]

        # Decide: explore or exploit?
        diversity = _mean_diversity([smi for smi, _ in pool])
        force_explore = plateau_count >= plateau_patience
        do_explore = force_explore or (diversity < diversity_threshold)
        mode = "EXPLORE" if do_explore else "EXPLOIT"

        if not quiet:
            best_raw = pool[0][1] / sign if pool else float('nan')
            print(f"Gen {gen+1}/{generations} | mode={mode} | pool={len(pool)} "
                  f"| winners={len(winners)} | diversity={diversity:.2f} "
                  f"| plateau={plateau_count} | best_pool={best_raw:.3f}")

        # If pool is exhausted, inject fresh mols and force explore
        if len(pool) == 0:
            if not quiet:
                print("  Pool empty — injecting fresh molecules")
            fresh = generate_filtered_molecules(profile, fragcorpus=fragcorpus,
                                                min_atoms=15, n_molecules=startN,
                                                n_jobs=n_jobs, seed=seed)
            raw_scores = _safe_score(scoring_function, fresh)
            for mol, rs in zip(fresh, raw_scores):
                pool.append((Chem.MolToSmiles(mol), rs * sign))
            pool = sorted(pool, key=lambda x: x[1])[:popsize]
            mode = "EXPLORE"

        parent_mols = [Chem.MolFromSmiles(smi) for smi, _ in pool]
        n_explore_random = max(2, int(startN * explore_ratio * 0.3))

        # ── generate new candidates ───────────────────────────────────────────
        if do_explore:
            new_mols = _explore(parent_mols, profile, fragcorpus, n_explore_random, seed,
                                conservative=conservative, n_jobs=n_jobs)
        else:
            new_mols = _exploit(parent_mols, profile, n_jobs=n_jobs)
            # Always add a small explore component even in exploit mode
            # so we never fully converge
            n_extra = max(1, int(len(parent_mols) * explore_ratio * 0.5))
            new_mols += _explore(parent_mols[:n_extra], profile, fragcorpus, 1, seed,
                                 conservative=conservative, n_jobs=n_jobs)

        if not new_mols:
            if not quiet:
                print("  No new molecules generated this generation")
            plateau_count += 1
            continue

        # ── score and sort into winners / pool ────────────────────────────────
        raw_scores = _safe_score(scoring_function, new_mols)
        existing_smis = {smi for smi, _ in pool} | {smi for smi, _ in winners}
        improved = False

        for mol, rs in zip(new_mols, raw_scores):
            smi = Chem.MolToSmiles(mol)
            if smi in existing_smis:
                continue
            existing_smis.add(smi)
            internal = rs * sign
            if internal < _win:
                fp = MFPGEN.GetFingerprint(mol)
                sim_to_winners = (max(DataStructs.BulkTanimotoSimilarity(fp, winnerfps))
                                  if winnerfps else 0.0)
                if sim_to_winners < sim_threshold:
                    winners.append((smi, rs))
                    winnerfps.append(fp)
                    improved = True
            else:
                pool.append((smi, internal))

        # Update plateau tracker using best score in pool
        new_best = min((s for _, s in pool), default=best_internal)
        if new_best < best_internal - 1e-6:
            best_internal = new_best
            plateau_count = 0
        else:
            plateau_count += 1 if not improved else 0

        # Fire callback with per-generation stats
        if callback is not None:
            best_pool_raw = (min(s for _, s in pool) / sign) if pool else None
            best_winner_raw = (sorted(winners, key=lambda x: x[1] * sign)[0][1]
                               if winners else None)
            callback({
                "generation":  gen + 1,
                "mode":        mode,
                "pool_size":   len(pool),
                "n_winners":   len(winners),
                "diversity":   diversity,
                "plateau":     plateau_count,
                "best_pool":   best_pool_raw,
                "best_winner": best_winner_raw,
            })

    winners = sorted(winners, key=lambda x: x[1] * sign)
    return winners


class GAReporter:
    """Collects per-generation statistics from :func:`generate_optimized_molecules` and plots them.

    Pass an instance as the ``callback=`` argument to the GA.  After the run,
    call :meth:`plot` to visualise the results, or compare multiple runs by
    passing a list of reporters to :func:`plot_runs`.

    Parameters
    ----------
    label : str
        Name for this run, shown in plot legends (default ``"run"``).

    Example
    -------
    ::

        reporter = GAReporter(label="explore_ratio=0.7")
        winners = gen.generate_optimized_molecules(
            my_score_fn, profile,
            generations=10, startN=30,
            callback=reporter,
        )
        reporter.plot()

    To compare settings::

        r1 = GAReporter("default")
        gen.generate_optimized_molecules(score_fn, p, callback=r1)

        r2 = GAReporter("high explore")
        gen.generate_optimized_molecules(score_fn, p, explore_ratio=0.8, callback=r2)

        GAReporter.compare([r1, r2])
    """

    def __init__(self, label="run"):
        self.label = label
        self.history = []   # list of stat dicts, one per generation

    def __call__(self, stats):
        """Called by the GA at the end of each generation with a stats dict.

        The dict contains:
        ``generation``, ``mode`` (``"EXPLORE"`` or ``"EXPLOIT"``),
        ``pool_size``, ``n_winners``, ``diversity``, ``plateau``,
        ``best_pool`` (best score in current pool),
        ``best_winner`` (best score among all winners so far, or None).
        """
        self.history.append(stats)

    # ── single-run plot ──────────────────────────────────────────────────────

    def plot(self, ax=None, show=True):
        """Plot score and mode history for this run.

        Parameters
        ----------
        ax   : matplotlib Axes or None — if None, a new figure is created
        show : bool — call ``plt.show()`` at the end (default True)

        Returns
        -------
        matplotlib Figure
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        if not self.history:
            print("No data recorded yet.")
            return None

        gens = [h["generation"] for h in self.history]
        best_pool = [h["best_pool"] for h in self.history]
        best_winner = [h["best_winner"] for h in self.history]
        diversity = [h["diversity"] for h in self.history]
        modes = [h["mode"] for h in self.history]
        n_winners = [h["n_winners"] for h in self.history]

        fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
        fig.suptitle(f"GA run: {self.label}", fontsize=13, fontweight="bold")

        # ── panel 1: scores ─────────────────────────────────────────────────
        ax1 = axes[0]
        ax1.plot(gens, best_pool, "o-", color="steelblue", label="best in pool", linewidth=1.5)
        # best_winner may be None for early gens
        bw_x = [g for g, bw in zip(gens, best_winner) if bw is not None]
        bw_y = [bw for bw in best_winner if bw is not None]
        if bw_x:
            ax1.plot(bw_x, bw_y, "s--", color="darkorange", label="best winner", linewidth=1.5)
        ax1.set_ylabel("Score")
        ax1.legend(loc="lower right", fontsize=9)
        ax1.set_title("Score per generation")
        # shade explore/exploit
        for i, (g, mode) in enumerate(zip(gens, modes)):
            color = "#d4eaf7" if mode == "EXPLORE" else "#fdebd0"
            ax1.axvspan(g - 0.5, g + 0.5, color=color, alpha=0.5, zorder=0)

        # ── panel 2: diversity & winners ────────────────────────────────────
        ax2 = axes[1]
        ax2.plot(gens, diversity, "^-", color="mediumpurple", label="diversity", linewidth=1.5)
        ax2_r = ax2.twinx()
        ax2_r.bar(gens, n_winners, color="salmon", alpha=0.5, label="winners")
        ax2_r.set_ylabel("# winners", color="salmon")
        ax2_r.tick_params(axis="y", labelcolor="salmon")
        ax2.set_ylabel("Mean Tanimoto distance")
        ax2.legend(loc="upper left", fontsize=9)
        ax2_r.legend(loc="upper right", fontsize=9)
        ax2.set_title("Population diversity & winner count")
        for i, (g, mode) in enumerate(zip(gens, modes)):
            color = "#d4eaf7" if mode == "EXPLORE" else "#fdebd0"
            ax2.axvspan(g - 0.5, g + 0.5, color=color, alpha=0.5, zorder=0)

        # ── panel 3: mode timeline ───────────────────────────────────────────
        ax3 = axes[2]
        mode_vals = [1 if m == "EXPLORE" else 0 for m in modes]
        ax3.step(gens, mode_vals, where="mid", color="gray", linewidth=1)
        ax3.fill_between(gens, mode_vals, step="mid", alpha=0.3, color="steelblue")
        ax3.set_yticks([0, 1])
        ax3.set_yticklabels(["EXPLOIT", "EXPLORE"])
        ax3.set_xlabel("Generation")
        ax3.set_title("Mode per generation")

        # legend patch for shading
        exp_patch = mpatches.Patch(color="#d4eaf7", alpha=0.8, label="EXPLORE")
        expl_patch = mpatches.Patch(color="#fdebd0", alpha=0.8, label="EXPLOIT")
        axes[0].legend(
            handles=[ax1.lines[0], ax1.lines[1] if bw_x else ax1.lines[0],
                     exp_patch, expl_patch],
            labels=["best in pool", "best winner", "EXPLORE", "EXPLOIT"],
            loc="lower right", fontsize=8
        )

        plt.tight_layout()
        if show:
            plt.show()
        return fig

    # ── multi-run comparison ─────────────────────────────────────────────────

    @staticmethod
    def compare(reporters, metric="best_pool", show=True):
        """Overlay score curves from multiple runs on a single plot.

        Parameters
        ----------
        reporters : list of GAReporter
        metric    : ``"best_pool"`` or ``"best_winner"`` (default ``"best_pool"``)
        show      : bool

        Returns
        -------
        matplotlib Figure
        """
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(9, 4))
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

        for i, rep in enumerate(reporters):
            if not rep.history:
                continue
            gens = [h["generation"] for h in rep.history]
            vals = [h[metric] for h in rep.history]
            if metric == "best_winner":
                vals = [v if v is not None else float("nan") for v in vals]
            color = colors[i % len(colors)]
            ax.plot(gens, vals, "o-", color=color, label=rep.label, linewidth=1.8)
            # shade modes
            modes = [h["mode"] for h in rep.history]
            for g, mode in zip(gens, modes):
                hatch = "/" if mode == "EXPLORE" else None
                ax.axvspan(g - 0.5, g + 0.5, color=color, alpha=0.06,
                           hatch=hatch, edgecolor=color, linewidth=0)

        ylabel = "Score (best in pool)" if metric == "best_pool" else "Score (best winner)"
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
        print(f"\n{'Gen':>4}  {'Mode':>7}  {'Pool':>5}  {'Win':>5}  {'Div':>6}  {'BestPool':>9}  {'BestWin':>9}")
        print("-" * 62)
        for h in self.history:
            bw = f"{h['best_winner']:.4f}" if h["best_winner"] is not None else "    —   "
            print(f"{h['generation']:>4}  {h['mode']:>7}  {h['pool_size']:>5}  "
                  f"{h['n_winners']:>5}  {h['diversity']:>6.3f}  "
                  f"{h['best_pool']:>9.4f}  {bw:>9}")
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
        callback=None):
    """Optimise a single seed molecule using fragment-level operations.

    Unlike :func:`generate_optimized_molecules`, which explores chemical space
    from random starting structures, this function starts from *one specific
    molecule* and iteratively improves it by applying ring, linker, substituent,
    and decoration operations — then retaining the best results as parents for
    the next round.

    This is the right tool when you have a **lead compound** (e.g. a known
    active, a docking hit, or a manually designed scaffold) and want to
    explore structural variations that improve a property while staying
    chemically reasonable.  The scoring function can be anything: an ML model
    from :mod:`lacan.gen`, a docking wrapper, a property calculator, etc.

    Algorithm
    ---------
    Each generation:

    1. Apply all four fragment operations (``replace_ring``,
       ``replace_substituent``, ``replace_linker``, ``decorate_scaffold``) to
       every molecule in the current beam — using ``_fragment_moves``.
    2. Score all candidates with ``scoring_function``.
    3. Separate winners (above ``win_threshold``) from non-winners.
    4. Keep the top ``beam_width`` non-winners as parents for the next round
       (beam search).
    5. Track plateau and stop early if no improvement for
       ``plateau_patience`` generations.

    Note: unlike the full GA there is no crossover (single-mol context), no
    random injection, and no EXPLOIT/EXPLORE switching — fragment moves alone
    are sufficient for lead optimisation.  For larger chemical-space sweeps
    use :func:`generate_optimized_molecules` with ``seed_mols=``.

    Parameters
    ----------
    mol              : RDKit Mol — the seed / lead molecule
    scoring_function : callable — ``fn(list[Mol]) -> list[float]``.
                       Higher or lower is better depending on
                       ``higher_is_better``.
    profile          : LACAN profile dict (from :func:`~lacan.lacan.load_profile`)
    generations      : int — maximum number of optimisation rounds (default 20)
    beam_width       : int — how many top molecules to keep as parents each
                       round (default 10)
    n_replacements   : int — fragment operation attempts per parent per round
                       (default 15)
    win_threshold    : float — score at which a molecule is considered a winner
                       and reported (default 0.8)
    higher_is_better : bool — True if higher scores are better (default True)
    conservative     : bool — prefer structurally similar fragment replacements
                       (default True)
    fragcorpus       : fragment corpus; loads ChEMBL default if None
    plateau_patience : int — stop early if best score hasn't improved for this
                       many consecutive generations (default 5)
    quiet            : bool — suppress progress output (default False)
    callback         : callable or None — called each generation with a stats
                       dict containing ``generation``, ``n_candidates``,
                       ``n_winners``, ``best_score``, ``plateau``

    Returns
    -------
    list of (smiles, score) tuples
        All winner molecules found during the run, sorted best-first.
        If no molecule exceeded ``win_threshold``, returns the best
        ``beam_width`` molecules seen across all generations.

    Example
    -------
    ::

        from rdkit import Chem
        from lacan import gen, lacan

        profile = lacan.load_profile("chembl")

        # Use the built-in LACAN score as the objective
        def lacan_score(mols):
            return [lacan.score_mol(m, profile)[0] for m in mols]

        lead = Chem.MolFromSmiles("CCCc1nn(C)c2c(=O)[nH]c(-c3ccccc3)nc12")
        results = gen.optimize_from_mol(lead, lacan_score, profile,
                                        generations=15, win_threshold=0.7)
        for smi, score in results[:5]:
            print(f"{score:.3f}  {smi}")

    Using an ML model from the gen module::

        # Assume you have trained an RF model on D3 actives:
        model_fn = gen.get_scoring_function(model, scaler)
        results = gen.optimize_from_mol(lead, model_fn, profile,
                                        win_threshold=0.6, higher_is_better=True)
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()

    sign = _sign(higher_is_better)
    _win = win_threshold * sign

    # Initialise beam with the seed molecule
    seed_smi = Chem.MolToSmiles(mol)
    seed_score_raw = _safe_score(scoring_function, [mol])[0]
    beam = [(seed_smi, seed_score_raw * sign)]   # (smiles, internal_score)
    winners = []
    all_seen_smis = {seed_smi}

    # Check if seed itself is already a winner
    if seed_score_raw * sign < _win:
        winners.append((seed_smi, seed_score_raw))

    best_internal = beam[0][1]
    plateau_count = 0

    if not quiet:
        print(f"optimize_from_mol: seed score = {seed_score_raw:.4f}")

    for gen_idx in range(generations):
        # Generate candidates from all beam members
        parent_mols = [Chem.MolFromSmiles(smi) for smi, _ in beam]
        candidates = []
        for parent in parent_mols:
            candidates += _fragment_moves(parent, profile, fragcorpus,
                                          n_replacements=n_replacements,
                                          conservative=conservative)

        # Deduplicate against everything seen so far
        novel = []
        for m in candidates:
            smi = Chem.MolToSmiles(m)
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

        # Separate winners from pool; update beam
        new_beam_entries = list(beam)  # carry forward existing beam members
        gen_best = best_internal

        for smi, rs in zip(novel_smis, raw_scores):
            internal = rs * sign
            if internal < _win:
                winners.append((smi, rs))
            new_beam_entries.append((smi, internal))
            if internal < gen_best:
                gen_best = internal

        # Sort beam by internal score (lower = better) and trim to beam_width
        new_beam_entries = sorted(new_beam_entries, key=lambda x: x[1])[:beam_width]
        beam = new_beam_entries

        # Plateau detection
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
                print(f"  Plateau reached — stopping early.")
            break

    # If no winners found, return the best molecules from beam as consolation
    if not winners:
        if not quiet:
            print("  No winners above threshold — returning best beam members.")
        winners = [(smi, sc / sign) for smi, sc in beam]

    # Sort winners best-first and deduplicate by SMILES
    seen, deduped = set(), []
    for smi, sc in sorted(winners, key=lambda x: x[1] * sign, reverse=higher_is_better):
        if smi not in seen:
            seen.add(smi)
            deduped.append((smi, sc))
    return deduped


def next_population(population):
    """Placeholder for a future incremental population-advance API.

    Not yet implemented.  Raises :exc:`NotImplementedError`.
    """
    raise NotImplementedError("next_population is not yet implemented")
