"""protect.py — Bond protection and SMARTS-based atom exclusion for LACAN.

This module provides two mechanisms:

**SMARTS-based atom exclusion** (stateless)
    Pass a ``protect_smarts`` SMARTS string to any operation in
    :mod:`lacan.replace`, :mod:`lacan.mutate`, or
    :func:`~lacan.gen.optimize_from_mol`.  Atoms matching the pattern are
    skipped by all mutation and fragment-replacement operations each call.
    Because the protected set is re-derived from the molecule on every call
    via :func:`get_protected_atoms`, it survives SMILES round-trips
    transparently, no atom properties are stored on the molecule.

**Bond protection** (stateful, stored as ``_lp`` bond property)
    Protected bonds are excluded from the LACAN score computation.  Use this
    when a molecule contains a required but chemically unusual motif (e.g. a
    Michael-acceptor covalent warhead) whose LACAN score would unfairly
    penalise the whole molecule.  Bond protection is stored as the ``_lp``
    boolean property on RDKit bond objects and is carried through all
    replacement operations via the internal ``_restore_bond_protection``
    helper in :mod:`lacan.replace`.

All functions return *new* molecules — the originals are never modified in place.

The :func:`mol_cleaner` utility uses bond protection internally to iteratively
repair LACAN-failing molecules while preserving the parts that already pass.
"""

from rdkit import Chem
from lacan import lacan

PROP = "_lp"
"""Name of the RDKit **bond** boolean property used to mark bond protection.

Note: atom protection is no longer stored as a mol property — use
``protect_smarts`` parameters instead (see :func:`get_protected_atoms`).
"""


# ---------------------------------------------------------------------------
# SMARTS-based atom exclusion helpers
# ---------------------------------------------------------------------------

def get_protected_atoms(mol, protect_smarts):
    """Return the set of atom indices matched by *protect_smarts* in *mol*.

    Call it once at the top of any function that needs to skip protected atoms,
    then pass the resulting frozenset to the per-atom checks.

    Parameters
    ----------
    mol            : RDKit Mol
    protect_smarts : str or None.  If None, returns an empty frozenset.

    Returns
    -------
    frozenset of int
        Atom indices matched by the SMARTS.  Empty frozenset if
        *protect_smarts* is ``None`` or ``""``.

    Raises
    ------
    ValueError
        If *protect_smarts* is not ``None`` and cannot be parsed.
    """
    if not protect_smarts:
        return frozenset()
    pattern = Chem.MolFromSmarts(protect_smarts)
    if pattern is None:
        raise ValueError(f"Invalid SMARTS: {protect_smarts!r}")
    return frozenset(
        idx for match in mol.GetSubstructMatches(pattern) for idx in match
    )


def reaction_touches_protected(mol, rxn, protect_smarts):
    """Return True if any match of *rxn* overlaps an atom matched by *protect_smarts*.

    Called by :func:`~lacan.mutate.apply_mutations` before each reaction to
    skip operations that would modify a protected atom.  Returns ``False``
    immediately when *protect_smarts* is ``None`` or ``""``, adding no
    overhead for unprotected molecules.

    Parameters
    ----------
    mol            : RDKit Mol
    rxn            : RDKit ChemicalReaction
    protect_smarts : str or None — SMARTS identifying atoms to exclude;
                     ``None`` disables the check entirely.

    Returns
    -------
    bool
        ``True`` if the reaction should be skipped; ``False`` if safe.
    """
    protected = get_protected_atoms(mol, protect_smarts)
    if not protected:
        return False
    for i in range(rxn.GetNumReactantTemplates()):
        template = rxn.GetReactantTemplate(i)
        for match in mol.GetSubstructMatches(template):
            if set(match) & protected:
                return True
    return False


# ---------------------------------------------------------------------------
# Bond protection — setting
# ---------------------------------------------------------------------------

def protect_bonds_for_idx(mol, bond_indices):
    """Mark specific bonds as protected.

    Protected bonds are excluded from the LACAN score (see
    :func:`score_mol_ignoring_protected_bonds`).  They do **not** block
    mutations — use ``protect_smarts`` in the operation functions for that.

    Parameters
    ----------
    mol          : RDKit Mol
    bond_indices : iterable of bond indices to protect

    Returns a new RDKit Mol.
    """
    rwmol = Chem.RWMol(mol)
    for idx in bond_indices:
        rwmol.GetBondWithIdx(idx).SetBoolProp(PROP, True)
    return rwmol.GetMol()


def protect_rejected_bonds(mol, profile=None, t=0.05):
    """Protect all bonds that currently fail the LACAN score threshold.

    After calling this, :func:`score_mol_ignoring_protected_bonds` will score
    the molecule as if those bonds do not exist.  This is useful when a
    molecule contains a structural motif that is required (e.g. a reactive
    warhead) but would otherwise cause the whole molecule to score 0.

    Parameters
    ----------
    mol     : RDKit Mol
    profile : LACAN profile dict; loads the default ChEMBL profile if None
    t       : bond score threshold (default 0.05)

    Returns a new RDKit Mol with failing bonds protected.
    """
    if profile is None:
        profile = lacan.load_profile("chembl")
    apb = lacan.assess_per_bond(mol, profile)
    bad_bond_indices = [i for i, score in enumerate(apb) if score < t]
    return protect_bonds_for_idx(mol, bad_bond_indices)


# ---------------------------------------------------------------------------
# Bond protection — querying
# ---------------------------------------------------------------------------

def get_protected_bond_indices(mol):
    """Return a list of indices of all protected bonds in *mol*."""
    return [b.GetIdx() for b in mol.GetBonds()
            if b.HasProp(PROP) and b.GetBoolProp(PROP)]


def bond_is_protected(bond):
    """Return True if the RDKit Bond object has the ``_lp`` protection mark."""
    return bond.HasProp(PROP) and bond.GetBoolProp(PROP)


# ---------------------------------------------------------------------------
# Scoring with bond exclusion
# ---------------------------------------------------------------------------

def score_mol_ignoring_protected_bonds(mol, profile=None, mode="score", t=0.05):
    """Score a molecule while ignoring any bonds marked as protected.

    This is a drop-in replacement for :func:`lacan.lacan.score_mol` that
    omits protected bonds from both the minimum-score calculation and the
    ``bad_bonds`` list.  If all bonds are protected the function returns 1.0
    (trivially passes).

    Parameters
    ----------
    mol     : RDKit Mol (may have protected bonds)
    profile : LACAN profile dict; loads ChEMBL default if None
    mode    : ``"score"`` (continuous, 0–1) or ``"threshold"`` (0 or 1)
    t       : bond score threshold (default 0.05)

    Returns
    -------
    (score, info) where ``info["bad_bonds"]`` lists unprotected failing bond indices.
    """
    if profile is None:
        profile = lacan.load_profile("chembl")
    apb = lacan.assess_per_bond(mol, profile)
    protected_bonds = set(get_protected_bond_indices(mol))
    apb_active = [score for i, score in enumerate(apb) if i not in protected_bonds]
    if not apb_active:
        return 1.0, {"bad_bonds": []}  # all bonds protected — trivially passes
    info = {}
    info["bad_bonds"] = [i for i, b in enumerate(apb) if b < t and i not in protected_bonds]
    if mode == "threshold":
        score = 0 if min(apb_active) < t else 1
    else:
        score = min(apb_active) / (1 + min(apb_active))
    return score, info


# ---------------------------------------------------------------------------
# mol_cleaner
# ---------------------------------------------------------------------------

def _raw_mutations(mol, profile):
    """Generate all valid single-step mutation products without score filtering.

    Unlike :func:`~lacan.mutate.apply_mutations`, this function does **not**
    score or filter the products.  It is used internally by :func:`mol_cleaner`
    so that partially-fixed intermediates (which still have some bad bonds) are
    not discarded before we get the chance to evaluate them ourselves.

    Degenerate products are rejected:
    * Disconnected molecules (contain ``"."`` in their canonical SMILES)
    * Molecules with zero bonds
    * Molecules with fewer than 4 heavy atoms.

    Parameters
    ----------
    mol     : RDKit Mol (typically has protected bonds from :func:`mol_cleaner`)
    profile : LACAN profile dict (unused here, kept for API consistency)

    Returns a list of sanitized, connected, non-degenerate RDKit Mol objects.
    """
    from lacan.mutate import mutate_ops

    products = []
    for rxn in mutate_ops.values():
        for prod_tuple in rxn.RunReactants((mol,)):
            try:
                m = prod_tuple[0]
                Chem.SanitizeMol(m)
                if m.GetNumBonds() == 0:
                    continue
                if m.GetNumHeavyAtoms() < 4:
                    continue
                if "." in Chem.MolToSmiles(m):
                    continue
                products.append(m)
            except Exception:
                pass
    return products


def mol_cleaner(mol, profile=None, score_threshold=0.5, t=0.05, max_iter=100,
                lateral_patience=5):
    """Iteratively mutate a molecule to eliminate all LACAN bond violations.

    This function is designed for molecules that mostly pass the LACAN profile
    but have one or more bad bond environments that need to be fixed.  It
    freezes the parts that already pass (using bond protection) and mutates
    only the failing regions.

    Strategy
    --------
    Each iteration the cleaner:

    1. Generates **all** single-step mutation products via :func:`_raw_mutations`
       — crucially *without* score-filtering the products.  Score-filtering at
       this stage would reject all partially-fixed intermediates (those that
       still have some bad bonds), preventing multi-step repair paths.

    2. Evaluates each product by counting its unprotected LACAN violations.
       The parent molecule's protected bond mask is re-derived from the
       *product's own* bond scores (via :func:`_reprotect`), so each candidate
       is assessed on its own bond landscape.

    3. Picks the best candidate:

       * **Improvement step** — a candidate with fewer violations than the
         current molecule.  Resets the lateral counter.
       * **Lateral move** — if no improvement is available and the lateral
         patience budget allows, accept the highest-scoring candidate with the
         *same* violation count.

    4. Re-protects the accepted candidate's newly-passing bonds so subsequent
       mutations stay focused on the remaining bad regions.

    Parameters
    ----------
    mol              : RDKit Mol to clean
    profile          : LACAN profile dict (loads ChEMBL default if None)
    score_threshold  : final LACAN score required to consider the molecule clean
                       (default 0.5)
    t                : bond score threshold used throughout (default 0.05)
    max_iter         : hard cap on total iterations (default 100)
    lateral_patience : consecutive lateral steps allowed before giving up
                       (default 5)

    Returns
    -------
    RDKit Mol if a clean version is found, else None.
    """
    if profile is None:
        profile = lacan.load_profile("chembl")

    min_heavy = max(4, mol.GetNumHeavyAtoms() - 6)

    def _reprotect(m):
        apb = lacan.assess_per_bond(m, profile)
        good = [i for i, sc in enumerate(apb) if sc >= t]
        return protect_bonds_for_idx(m, good)

    def _eval(m):
        if m.GetNumBonds() == 0:
            return 999, 0.0
        _, info = score_mol_ignoring_protected_bonds(m, profile, t=t)
        s, _ = score_mol_ignoring_protected_bonds(m, profile, t=t)
        return len(info["bad_bonds"]), s

    working_mol = _reprotect(mol)
    current_v, current_s = _eval(working_mol)

    lateral_count = 0
    visited_smis = {Chem.MolToSmiles(working_mol)}

    for _ in range(max_iter):
        if current_v == 0 and current_s >= score_threshold:
            return working_mol

        candidates = _raw_mutations(working_mol, profile)
        if not candidates:
            break

        best_improve_v = current_v
        best_improve_s = -1.0
        best_improve_mol = None
        best_lateral_s = -1.0
        best_lateral_mol = None

        for c in candidates:
            if c.GetNumHeavyAtoms() < min_heavy:
                continue
            smi = Chem.MolToSmiles(c)
            if smi in visited_smis:
                continue
            c_prot = _reprotect(c)
            v, s = _eval(c_prot)

            if v < best_improve_v or (v == best_improve_v and s > best_improve_s):
                best_improve_v = v
                best_improve_s = s
                best_improve_mol = c_prot

            if v == current_v and s > best_lateral_s:
                best_lateral_s = s
                best_lateral_mol = c_prot

        if best_improve_mol is not None and best_improve_v < current_v:
            visited_smis.add(Chem.MolToSmiles(best_improve_mol))
            working_mol = best_improve_mol
            current_v = best_improve_v
            current_s = best_improve_s
            lateral_count = 0
        elif best_lateral_mol is not None and lateral_count < lateral_patience:
            visited_smis.add(Chem.MolToSmiles(best_lateral_mol))
            working_mol = best_lateral_mol
            current_s = best_lateral_s
            lateral_count += 1
        else:
            break

    final_score, _ = score_mol_ignoring_protected_bonds(working_mol, profile, t=t)
    return working_mol if final_score >= score_threshold else None
