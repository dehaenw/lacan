"""
replace.py — Coarse-grained fragment replacement operations for LACAN.

This module provides four scaffold-editing operations that swap entire fragments
(rings, substituents, linkers) for alternatives drawn from the fragment corpus
(rls.csv).  Together they form the EXPLORE step of the adaptive GA in gen.py.

All four operations share these behaviours:

* **Random site selection** — when a molecule has multiple rings, linkers, or
  substituents that could be operated on, one is chosen at random each call.
  This ensures diverse output across repeated calls and prevents the same site
  from being targeted every time.
* **Conservative sampling** (``conservative=True``, the default) — replacement
  fragments are drawn with probability proportional to both their corpus
  frequency *and* their structural similarity to the fragment being replaced
  (Morgan fp Tanimoto, radius 2, including dummy atoms).  Set
  ``conservative=False`` for purely frequency-weighted sampling.
* **Protection** — bonds touching protected atoms are never cut; ring systems
  containing protected atoms are skipped entirely.
* **Deduplication** — output molecules are deduplicated by InChIKey and scored
  by the LACAN profile; molecules scoring 0 are discarded.

Reaction conventions
--------------------
The ``_breakbond`` SMARTS ``[!#0:0]!@[R:1]`` matches exo-bonds as the pair
``(non-ring atom, ring atom)``, so in every match tuple ``b[0]`` is the
**exo** atom and ``b[1]`` is the **ring** atom.  This is important for
:func:`replace_ring` and :func:`replace_linker` which use these indices to
identify which atoms belong to the scaffold vs the side chains.
"""

from rdkit import Chem
from lacan import decompose as decompose_module
from lacan import lacan
from lacan.protect import (
    atom_is_protected,
    score_mol_ignoring_protected_bonds,
)
import random
import csv
import os
from rdkit import DataStructs
from rdkit.Chem import rdChemReactions, inchi, rdFingerprintGenerator

# ---------------------------------------------------------------------------
# Reactions
# ---------------------------------------------------------------------------

rxn1 = rdChemReactions.ReactionFromSmarts(
    "[*:0]-[#0:1].[*:2]-[#0:3]>>[*:0]-[*:2].[*:1]-[*:3]"
)
"""Connect two single-bond dummy attachment points."""

rxn2 = rdChemReactions.ReactionFromSmarts(
    "[*:0]=[#0:1].[*:2]=[#0:3]>>[*:0]=[*:2].[*:1]=[*:3]"
)
"""Connect two double-bond dummy attachment points."""

htodummy = rdChemReactions.ReactionFromSmarts("[H1,H2,H3:0]>>[*:0]-[#0]")
"""Convert one hydrogen to a single-bond dummy attachment point."""

decompose1 = rdChemReactions.ReactionFromSmarts(
    "[!#0:0]!@;-[R:1]>>[!R:0]-[#0].[R:1]-[#0]"
)
"""Cut a single-bond ring-exo bond, yielding (exo-fragment, ring-fragment)."""

singledummy = Chem.MolFromSmarts("[#0;$([*]-[*])]")
"""Matches a dummy atom attached by a single bond."""

doubledummy = Chem.MolFromSmarts("[#0;$([*]=[*])]")
"""Matches a dummy atom attached by a double bond."""

# _breakbond matches (non-ring-atom, ring-atom) pairs.
# In every match b=(b[0], b[1]): b[0] = exo atom, b[1] = ring atom.
_breakbond = Chem.MolFromSmarts("[!#0:0]!@[R:1]")

# ---------------------------------------------------------------------------
# Fragment corpus
# ---------------------------------------------------------------------------

_FRAGFPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
"""Morgan generator (r=2) used for fragment similarity scoring.

Including dummy atoms at radius 2 captures the immediate attachment
environment, giving meaningful similarity values for fragments that
differ mainly in the context around their attachment points.
"""

_corpus_cache = {}  # keyed by min_count so different thresholds coexist


def load_corpus(min_count=200, path=None):
    """Load a fragment corpus CSV, with caching.

    By default loads the built-in ChEMBL corpus from ``data/rls.csv``.  Pass
    ``path`` to load a custom corpus built with ``python -m lacan.decompose``
    (e.g. from COCONUT or your own dataset).

    The CSV is read only on the first call for a given ``(path, min_count)``
    combination; subsequent calls return the cached list immediately.

    Parameters
    ----------
    min_count : int
        Minimum occurrence count for a fragment to be included (default 200).
        Lower values allow rarer fragments into the replacement pool.
    path : str or None
        Path to a corpus CSV produced by :mod:`lacan.decompose`.  If ``None``
        (default) the built-in ChEMBL corpus is used.

    Returns
    -------
    list of [smiles, count, degree, ftype, bonds]
    """
    cache_key = (path, min_count)
    if cache_key in _corpus_cache:
        return _corpus_cache[cache_key]
    if path is None:
        _this_dir = os.path.dirname(__file__)
        path = os.path.join(_this_dir, "data/rls.csv")
    _entries = []
    with open(path, newline="") as _f:
        _reader = csv.reader(_f, delimiter=",", quotechar='"')
        for _e in _reader:
            if _e[0] != "smiles" and int(_e[1]) > min_count:
                _entries.append([_e[0], int(_e[1]), int(_e[2]), _e[3], _e[4]])
    _corpus_cache[cache_key] = _entries
    return _entries

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _frag_fp(smi):
    """Return a Morgan fingerprint for a fragment SMILES string, or None on failure.

    Dummy atoms (atomic number 0) are treated as regular heavy atoms, so the
    fingerprint captures the attachment-point chemical environment — exactly
    what we want when comparing fragments by their connection context.
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return _FRAGFPGEN.GetFingerprint(mol)


def _similarity_weights(query_smi, candidates, base_weight_power=0.5, sim_boost=2.0):
    """Compute conservative sampling weights for a list of replacement fragment entries.

    The weight assigned to each candidate is::

        weight = count ^ base_weight_power * (1 + sim_boost * tanimoto(query, candidate))

    This blends two signals:

    * **Frequency** (``count ^ base_weight_power``) — common corpus fragments
      are preferred, avoiding rare exotic scaffolds.
    * **Similarity** (Tanimoto of Morgan fps) — fragments resembling the one
      being replaced get a bonus, producing conservative bioisosteric swaps.

    Parameters
    ----------
    query_smi         : SMILES of the fragment currently in the molecule
    candidates        : list of corpus entries ``[smiles, count, degree, ftype, bonds]``
    base_weight_power : exponent applied to occurrence count (default 0.5)
    sim_boost         : multiplier for the similarity bonus term (default 2.0)

    Returns a list of floats parallel to ``candidates``.
    """
    qfp = _frag_fp(query_smi)
    if qfp is None:
        return [e[1] ** base_weight_power for e in candidates]
    cfps = [_frag_fp(e[0]) for e in candidates]
    weights = []
    for i, e in enumerate(candidates):
        occ_w = e[1] ** base_weight_power
        sim = DataStructs.TanimotoSimilarity(qfp, cfps[i]) if cfps[i] is not None else 0.0
        weights.append(occ_w * (1.0 + sim_boost * sim))
    return weights


def _filter_and_dedup(new_mols, p):
    """Sanitize, LACAN-score, and deduplicate a list of product molecules.

    A molecule is discarded if:

    * RDKit sanitization fails,
    * the LACAN score is exactly 0 (at least one bond is a hard alert), or
    * it shares an InChIKey with a molecule already in the output list.

    Protected bonds on the products are honoured by
    :func:`~lacan.protect.score_mol_ignoring_protected_bonds`.

    Returns a list of valid, unique RDKit Mol objects.
    """
    filtered = []
    iks = []
    for m in new_mols:
        try:
            Chem.SanitizeMol(m)
            score, _ = score_mol_ignoring_protected_bonds(m, p)
            if score > 0:
                ik = inchi.MolToInchiKey(m)
                if ik not in iks:
                    filtered.append(m)
                    iks.append(ik)
        except Exception:
            pass
    return filtered


def _get_ring_systems(mol):
    """Return fused ring systems as a list of frozensets of atom indices.

    Simple rings that share at least one atom (fused or bridged systems) are
    merged into a single frozenset, so naphthalene yields one system of 10
    atoms rather than two overlapping 6-membered rings.
    """
    ring_atom_sets = [set(r) for r in mol.GetRingInfo().AtomRings()]
    merged = True
    while merged:
        merged = False
        new_sets = []
        while ring_atom_sets:
            current = ring_atom_sets.pop(0)
            fused = False
            for i, other in enumerate(ring_atom_sets):
                if current & other:
                    ring_atom_sets[i] = current | other
                    fused = True
                    merged = True
                    break
            if not fused:
                new_sets.append(frozenset(current))
        ring_atom_sets = new_sets
    return ring_atom_sets


def _restore_bond_protection(product, orig_mol):
    """Re-apply ``_lp`` bond protection from *orig_mol* onto *product*.

    RDKit reactions do not carry bond properties from reactants to products, so
    ``_lp`` marks are lost after stitching.  This function restores them by
    comparing atom-index pairs: if atoms *i* and *j* are bonded in *product*
    and the same pair was a protected bond in *orig_mol*, the product bond is
    re-marked ``_lp=True``.

    ``FragmentOnBonds`` preserves original atom indices for all non-dummy atoms,
    and the stitching reactions only form new bonds at dummy-atom sites, so
    heavy-atom indices are stable throughout.

    Parameters
    ----------
    product  : RDKit Mol (assembled replacement product, may be RWMol)
    orig_mol : RDKit Mol with ``_lp`` bond properties set

    Returns a new RDKit Mol with protection restored (or *product* unchanged if
    *orig_mol* has no protected bonds).
    """
    # Build set of atom-pair tuples for protected bonds in the original molecule
    protected_pairs = set()
    for b in orig_mol.GetBonds():
        if b.HasProp("_lp") and b.GetBoolProp("_lp"):
            protected_pairs.add(
                tuple(sorted([b.GetBeginAtomIdx(), b.GetEndAtomIdx()]))
            )
    if not protected_pairs:
        return product

    rwmol = Chem.RWMol(product)
    for b in rwmol.GetBonds():
        pair = tuple(sorted([b.GetBeginAtomIdx(), b.GetEndAtomIdx()]))
        if pair in protected_pairs:
            b.SetBoolProp("_lp", True)
    return rwmol.GetMol()



    """Return fused ring systems as a list of frozensets of atom indices.

    Simple rings that share at least one atom (fused or bridged systems) are
    merged into a single frozenset, so naphthalene yields one system of 10
    atoms rather than two overlapping 6-membered rings.

    Using RDKit's ``GetRingInfo().AtomRings()`` directly (rather than matching
    decompose-derived SMARTS back onto the molecule) avoids the ambiguity where
    the dummy ``*`` atom in a decompose fragment SMARTS can match non-ring atoms
    in the original molecule.
    """
    ring_atom_sets = [set(r) for r in mol.GetRingInfo().AtomRings()]
    merged = True
    while merged:
        merged = False
        new_sets = []
        while ring_atom_sets:
            current = ring_atom_sets.pop(0)
            fused = False
            for i, other in enumerate(ring_atom_sets):
                if current & other:
                    ring_atom_sets[i] = current | other
                    fused = True
                    merged = True
                    break
            if not fused:
                new_sets.append(frozenset(current))
        ring_atom_sets = new_sets
    return ring_atom_sets


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def decorate_scaffold(mol, p, score_threshold, fragcorpus=None, n_replacements=100,
                      mode="Hydrogen", replacements_per_mol=1, chance_of_linker=0.5):
    """Add substituents to a scaffold by decorating free positions.

    Two modes are supported:

    **Hydrogen mode** (default)
        For each of the ``n_replacements`` attempts, one H-bearing atom is
        chosen at random and its hydrogen replaced with a dummy, then a
        single-attachment corpus fragment is joined at that position.
        ``replacements_per_mol`` controls how many H → fragment substitutions
        are made per attempt.  Protected atoms are skipped when selecting the
        H position.

    **Dummy mode**
        The scaffold must already contain ``*`` dummy atoms marking desired
        decoration sites (e.g. ``c1c(*)c(*)co1``).  Each dummy is filled with
        a corpus fragment.  With probability ``chance_of_linker`` a two-
        attachment linker fragment is inserted before the terminal substituent,
        extending the chain length.

    Parameters
    ----------
    mol               : RDKit Mol (scaffold to decorate)
    p                 : LACAN profile dict
    score_threshold   : minimum LACAN score required for output molecules;
                        set 0.0 in the GA explore phase to accept all structures
    fragcorpus        : list of corpus entries (default: ChEMBL rls.csv)
    n_replacements    : number of decorated molecules to attempt generating
    mode              : ``"Hydrogen"`` or ``"Dummy"``
    replacements_per_mol : (Hydrogen mode only) H → fragment substitutions per attempt
    chance_of_linker  : (Dummy mode only) probability of prepending a linker fragment

    Returns
    -------
    list of RDKit Mol
        Deduplicated molecules that pass ``score_threshold``.
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()
    fs = [e for e in fragcorpus if e[2] == 1 and e[4] == "-"]
    fd = [e for e in fragcorpus if e[2] == 1 and e[4] == "="]
    fl = [e for e in fragcorpus if e[2] == 2 and e[3] == "Linker" and e[4] == "--"]
    new_mols = []

    if mode == "Dummy":
        sdb_matches = [m for m in mol.GetSubstructMatches(singledummy)
                       if not atom_is_protected(mol.GetAtomWithIdx(m[0]))]
        ddb_matches = [m for m in mol.GetSubstructMatches(doubledummy)
                       if not atom_is_protected(mol.GetAtomWithIdx(m[0]))]
        sdb = len(sdb_matches)
        ddb = len(ddb_matches)
        for i in range(n_replacements):
            cmol = Chem.Mol(mol)
            fe1 = [Chem.MolFromSmiles(e[0]) for e in random.choices(fs, [e[1] ** 0.5 for e in fs], k=sdb)]
            fe2 = [Chem.MolFromSmiles(e[0]) for e in random.choices(fd, [e[1] ** 0.5 for e in fd], k=ddb)]
            fe3 = [Chem.MolFromSmiles(e[0]) for e in random.choices(fl, [e[1] ** 0.5 for e in fl], k=sdb)]
            for j in range(sdb):
                if random.random() > chance_of_linker:
                    cmol = random.choice(rxn1.RunReactants((cmol, fe3[j])))[0]
                    Chem.SanitizeMol(cmol)
                cmol = random.choice(rxn1.RunReactants((cmol, fe1[j])))[0]
                Chem.SanitizeMol(cmol)
            for j in range(ddb):
                cmol = random.choice(rxn2.RunReactants((cmol, fe2[j])))[0]
                Chem.SanitizeMol(cmol)
            new_mols.append(_restore_bond_protection(cmol, mol))

    elif mode == "Hydrogen":
        for i in range(n_replacements):
            cmol = Chem.Mol(mol)
            for j in range(replacements_per_mol):
                candidates = [a.GetIdx() for a in cmol.GetAtoms()
                              if a.GetTotalNumHs() > 0 and not atom_is_protected(a)]
                if not candidates:
                    break
                chosen = random.choice(candidates)
                prods = htodummy.RunReactants((cmol,))
                valid = []
                for prod in prods:
                    try:
                        Chem.SanitizeMol(prod[0])
                        for a in prod[0].GetAtoms():
                            if a.GetAtomicNum() == 0:
                                if any(n.GetIdx() == chosen for n in a.GetNeighbors()):
                                    valid.append(prod[0])
                                    break
                    except Exception:
                        pass
                if valid:
                    cmol = random.choice(valid)
                else:
                    if prods:
                        cmol = random.choice(prods)[0]
                        Chem.SanitizeMol(cmol)
            fe = [Chem.MolFromSmiles(e[0]) for e in random.choices(
                fs, [e[1] ** 0.5 for e in fs], k=replacements_per_mol)]
            for j in range(replacements_per_mol):
                cmol = random.choice(rxn1.RunReactants((cmol, fe[j])))[0]
                Chem.SanitizeMol(cmol)
            new_mols.append(_restore_bond_protection(cmol, mol))
    else:
        print("mode not supported")

    filtered_prods = []
    iks = []
    for m in new_mols:
        try:
            Chem.SanitizeMol(m)
            score, _ = score_mol_ignoring_protected_bonds(m, p)
            if score > score_threshold:
                ik = inchi.MolToInchiKey(m)
                if ik not in iks:
                    filtered_prods.append(m)
                    iks.append(ik)
        except Exception:
            pass
    return filtered_prods


def replace_substituent(mol, p, score_threshold, fragcorpus=None, n_replacements=100,
                        conservative=True):
    """Replace one non-ring substituent with a randomly sampled corpus fragment.

    A substituent is defined here strictly: a fragment produced by cutting a
    single ring-exo bond that (a) has exactly one attachment point (one ``*``
    dummy) and (b) contains no ring atoms from the original molecule.  This
    prevents entire ring systems from being mistaken for substituents when the
    SMARTS fires on inter-ring bonds.

    Algorithm
    ---------
    1. Run ``decompose1`` (``[!#0]!@;-[R]>>[!R]-[*].[R]-[*]``) on the
       molecule to enumerate all possible ring-exo single-bond cuts.
    2. For each cut, classify the two fragments: smaller = substituent,
       larger = scaffold.
    3. **Filter**: discard any pair where the smaller fragment contains ring
       atoms from the original molecule (it is a ring system, not a sub).
    4. Also skip if the scaffold's attachment atom is protected.
    5. **Randomly pick** one valid substituent site for each of the
       ``n_replacements`` attempts, so all sites are sampled across repeated
       calls.
    6. Draw a replacement from corpus substituents (degree 1, single bond),
       optionally similarity-weighted, and join it to the scaffold dummy.

    Parameters
    ----------
    mol             : RDKit Mol
    p               : LACAN profile dict
    score_threshold : minimum LACAN score for output molecules
    fragcorpus      : fragment corpus
    n_replacements  : number of replacement attempts
    conservative    : if True, bias sampling toward structurally similar fragments

    Returns
    -------
    list of RDKit Mol
    """
    # Pre-compute ring atom indices for the original molecule so we can reject
    # "substituents" that actually contain ring atoms.
    if fragcorpus is None:
        fragcorpus = load_corpus()
    ring_atoms = set()
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms.update(ring)

    prods = decompose1.RunReactants((mol,))
    pairs = []      # (substituent_mol, scaffold_mol) tuples
    sub_smis = []   # SMILES of each substituent for similarity weighting

    for prod in prods:
        if len(prod) != 2:
            continue
        try:
            [Chem.SanitizeMol(pp) for pp in prod]
        except Exception:
            continue
        sorted_prod = sorted(prod, key=lambda m: m.GetNumAtoms())
        substituent, scaffold = sorted_prod[0], sorted_prod[1]

        # Reject if the smaller fragment has ring atoms in the original molecule
        # (it is an inter-ring bond cut, not a true ring-exo substituent cut).
        sub_match = mol.GetSubstructMatch(substituent)
        if sub_match and any(idx in ring_atoms for idx in sub_match
                             if mol.GetAtomWithIdx(idx).GetAtomicNum() != 0):
            continue

        # Exactly one dummy and no rings in fragment
        sub_smi = Chem.MolToSmiles(substituent)
        if sub_smi.count("*") != 1:
            continue
        if substituent.GetRingInfo().NumRings() > 0:
            continue

        # Skip if the scaffold attachment point (neighbour of dummy) is protected
        dummy_neighbor = None
        for a in scaffold.GetAtoms():
            if a.GetAtomicNum() == 0:
                for n in a.GetNeighbors():
                    dummy_neighbor = n
        if dummy_neighbor and atom_is_protected(dummy_neighbor):
            continue

        pairs.append((substituent, scaffold))
        sub_smis.append(sub_smi)

    if not pairs:
        return []

    fe = [e for e in fragcorpus if e[2] == 1 and e[4] == "-"]
    if not fe:
        return []

    new_mols = []
    for i in range(n_replacements):
        # Randomly pick one substituent site per attempt
        pair_idx = random.randrange(len(pairs))
        scaffold = pairs[pair_idx][1]
        if conservative:
            weights = _similarity_weights(sub_smis[pair_idx], fe)
        else:
            weights = [e[1] ** 0.5 for e in fe]
        replacement = random.choices(fe, weights=weights, k=1)[0]
        frag_mol = Chem.MolFromSmiles(replacement[0])
        try:
            new_mol = rxn1.RunReactants((frag_mol, scaffold))[0][0]
            new_mols.append(_restore_bond_protection(new_mol, mol))
        except Exception:
            pass
    return _filter_and_dedup(new_mols, p)


def replace_linker(mol, p, score_threshold, fragcorpus=None, n_replacements=100,
                   conservative=True):
    """Replace one non-ring linker connecting two ring systems.

    A linker is a non-ring fragment with exactly two attachment points (degree 2)
    that bridges two ring systems.  The ``decompose_molecule`` utility identifies
    them as SMARTS strings with two ``*`` dummies.

    Algorithm
    ---------
    1. Decompose the molecule to get linker SMARTS strings.
    2. For each linker SMARTS, find all substructure matches in the original mol.
    3. For each match, identify the two exo-bonds using ``_breakbond``: atoms in
       the linker match set whose ``b[0]`` (the exo/linker atom) appears in the
       match are the cut points.  Exactly two such bonds are required.
    4. Collect all valid sites (linker SMARTS + two ring fragments after cutting).
    5. **Randomly pick one site** from all valid sites so that molecules with
       multiple linkers sample them uniformly across calls.
    6. Draw ``n_replacements`` replacement linkers from the corpus (degree 2,
       type "Linker", bond string "``--``") and stitch the two ring fragments
       back through each.

    ``_breakbond`` convention: ``b[0]`` = non-ring (linker) atom,
    ``b[1]`` = ring atom.  Only ``b[0]`` is checked against the linker match
    set to locate the two linker-to-ring attachment bonds.

    Parameters
    ----------
    mol             : RDKit Mol
    p               : LACAN profile dict
    score_threshold : minimum LACAN score for output molecules
    fragcorpus      : fragment corpus
    n_replacements  : number of replacement linkers to try at the chosen site
    conservative    : if True, similarity-bias the linker sampling

    Returns
    -------
    list of RDKit Mol, or ``[]`` if no replaceable linker is found.
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()
    dcmol = decompose_module.decompose_molecule(Chem.MolToSmiles(mol))
    linkers_smarts = dcmol[1]
    if not linkers_smarts:
        return []

    breakable = mol.GetSubstructMatches(_breakbond)
    # _breakbond: b[0] = exo/linker atom, b[1] = ring atom

    # Collect all valid linker sites across all linker types and matches
    valid_sites = []  # list of (linker_smi, [frag1, frag2])

    for linker_smi in linkers_smarts:
        linker_pattern = Chem.MolFromSmarts(linker_smi)
        if linker_pattern is None:
            continue
        for linker_match in mol.GetSubstructMatches(linker_pattern):
            linker_set = set(linker_match)
            # b[0] is the linker (non-ring) atom — find breakable bonds where
            # the exo atom belongs to this linker
            to_break = [b for b in breakable if b[0] in linker_set]
            if len(to_break) != 2:
                continue
            if any(atom_is_protected(mol.GetAtomWithIdx(idx))
                   for pair in to_break for idx in pair):
                continue
            bond_indices = [mol.GetBondBetweenAtoms(*b).GetIdx() for b in to_break]
            cmol = Chem.FragmentOnBonds(mol, bond_indices)
            # Keep the two ring-side fragments (each has exactly 1 dummy)
            frags = [Chem.MolFromSmiles(f) for f in Chem.MolToSmiles(cmol).split(".")
                     if f.count("*") == 1]
            if len(frags) != 2:
                continue
            valid_sites.append((linker_smi, frags))

    if not valid_sites:
        return []

    # Randomly pick one linker site so multi-linker molecules sample uniformly
    linker_smi, frags = random.choice(valid_sites)

    fe = [e for e in fragcorpus if e[2] == 2 and e[3] == "Linker" and e[4] == "--"]
    if not fe:
        return []

    if conservative:
        weights = _similarity_weights(linker_smi, fe)
    else:
        weights = [e[1] ** 0.5 for e in fe]
    fr = [Chem.MolFromSmiles(e[0]) for e in random.choices(fe, weights=weights, k=n_replacements)]

    new_mols = []
    for new_linker in fr:
        try:
            combined = random.choice(rxn1.RunReactants((new_linker, frags[0])))[0]
            Chem.SanitizeMol(combined)
            combined = random.choice(rxn1.RunReactants((combined, frags[1])))[0]
            Chem.SanitizeMol(combined)
            new_mols.append(_restore_bond_protection(combined, mol))
        except Exception:
            pass
    return _filter_and_dedup(new_mols, p)


def replace_ring(mol, p, score_threshold, fragcorpus=None, n_replacements=100,
                 conservative=True):
    """Replace one ring system with a different ring drawn from the corpus.

    Algorithm
    ---------
    1. Find fused ring systems via :func:`_get_ring_systems`.
    2. For each ring system walk bonds to find exo-bonds (leaving the system,
       not themselves in a ring).
    3. Skip protected systems; randomly pick one.
    4. Cut the exo-bonds with ``FragmentOnBonds``.  Use ``GetMolFrags()``
       (without ``asMols``) to get tuples of *original* atom indices for each
       fragment — these are stable and comparable to ``ring_system``.  The
       fragment whose non-dummy atoms are all in ``ring_system`` is the old ring;
       everything else is a side chain to reattach.  Build side-chain Mol objects
       with ``GetMolFrags(asMols=True)`` in the same order.
    5. Sample corpus rings; stitch side chains back.
    """
    if fragcorpus is None:
        fragcorpus = load_corpus()

    ring_systems = _get_ring_systems(mol)
    if not ring_systems:
        return []

    valid = []
    for ring_system in ring_systems:
        if any(atom_is_protected(mol.GetAtomWithIdx(i)) for i in ring_system):
            continue
        seen = set()
        exo_bonds = []
        for idx in ring_system:
            for bond in mol.GetAtomWithIdx(idx).GetBonds():
                bid = bond.GetIdx()
                if bid in seen:
                    continue
                seen.add(bid)
                if bond.IsInRing():
                    continue
                other = bond.GetOtherAtomIdx(idx)
                if other in ring_system:
                    continue
                btype = "-" if bond.GetBondTypeAsDouble() == 1.0 else "="
                exo_bonds.append((bid, btype))
        if not exo_bonds:
            continue
        ring_smi = Chem.MolFragmentToSmiles(mol, list(ring_system))
        valid.append((ring_system, exo_bonds, ring_smi))

    if not valid:
        return []

    ring_system, exo_bonds, ring_smi = random.choice(valid)
    degree     = len(exo_bonds)
    bond_idxs  = [e[0] for e in exo_bonds]
    bond_types = [e[1] for e in exo_bonds]

    # FragmentOnBonds adds new dummy atoms at indices >= mol.GetNumAtoms().
    # GetMolFrags() (no asMols) returns tuples of original-mol atom indices,
    # so we can compare directly against ring_system.
    cmol = Chem.FragmentOnBonds(mol, bond_idxs)
    orig_n = mol.GetNumAtoms()

    frag_idx_tuples = Chem.GetMolFrags(cmol)          # tuples of atom indices in cmol
    frag_mols       = Chem.GetMolFrags(cmol, asMols=True, sanitizeFrags=False)

    ring_frag_pos = None
    for i, idx_tuple in enumerate(frag_idx_tuples):
        # Non-dummy atoms in this fragment: those with index < orig_n
        heavy = [idx for idx in idx_tuple if idx < orig_n]
        if all(idx in ring_system for idx in heavy):
            ring_frag_pos = i
            break

    if ring_frag_pos is None:
        return []

    side_chains = [m for i, m in enumerate(frag_mols) if i != ring_frag_pos]

    if len(side_chains) != degree:
        return []

    sbc = bond_types.count("-")
    dbc = bond_types.count("=")
    fe = [e for e in fragcorpus
          if e[2] == degree and e[3] == "Ring"
          and e[4].count("-") == sbc and e[4].count("=") == dbc]
    if not fe:
        return []

    weights = _similarity_weights(ring_smi, fe) if conservative else [e[1] ** 0.5 for e in fe]

    new_mols = []
    for entry in random.choices(fe, weights=weights, k=n_replacements):
        new_ring = Chem.MolFromSmiles(entry[0])
        if new_ring is None:
            continue
        try:
            ok = True
            for sc, btype in zip(side_chains, bond_types):
                rxn = rxn1 if btype == "-" else rxn2
                products = rxn.RunReactants((new_ring, sc))
                if not products:
                    ok = False
                    break
                new_ring = random.choice(products)[0]
                Chem.SanitizeMol(new_ring)
            if ok:
                new_mols.append(_restore_bond_protection(new_ring, mol))
        except Exception:
            pass

    return _filter_and_dedup(new_mols, p)


if __name__ == "__main__":
    random.seed(123)
    decorate_scaffold(Chem.MolFromSmiles("c1c(*)c(*)co1"), mode="Dummy")
