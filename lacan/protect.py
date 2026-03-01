"""
protect.py — Atom and bond protection for LACAN.

This module provides two orthogonal protection mechanisms:

**Atom protection**
    Protected atoms are skipped by all mutation (mutate.py) and fragment
    replacement (replace.py) operations.  Any reaction template that would
    match a protected atom is skipped entirely.  Use this to freeze a
    pharmacophore, warhead, or other structural feature while allowing the
    rest of the molecule to be explored.

**Bond protection**
    Protected bonds are excluded from the LACAN score computation.  Use this
    when a molecule contains a required but chemically unusual motif (e.g. a
    Michael-acceptor covalent warhead) whose LACAN score would unfairly
    penalise the whole molecule.

Both use the same RDKit atom/bond property ``_lp`` (a bool) so the footprint
on the molecule is minimal.  All functions return *new* molecules — the
originals are never modified in place.

The :func:`mol_cleaner` utility combines both mechanisms to iteratively repair
LACAN-failing molecules while preserving the parts that already pass.
"""

from rdkit import Chem
from lacan import lacan

PROP = "_lp"
"""Name of the RDKit atom/bond boolean property used to mark protection."""


# ---------------------------------------------------------------------------
# Setting protection
# ---------------------------------------------------------------------------

def protect_atoms_for_idx(mol, indices):
    """Mark specific atoms as protected.

    Protected atoms are ignored by all mutation and replacement reactions —
    any reaction whose SMARTS template would match a protected atom is
    skipped.  Downstream operations read the ``_lp`` property to enforce this.

    Parameters
    ----------
    mol     : RDKit Mol (not modified in place)
    indices : iterable of atom indices to protect

    Returns a new RDKit Mol with the ``_lp`` property set on the listed atoms.
    """
    rwmol = Chem.RWMol(mol)
    for idx in indices:
        rwmol.GetAtomWithIdx(idx).SetBoolProp(PROP, True)
    return rwmol.GetMol()


def protect_atoms_matching_smarts(mol, smarts):
    """Mark all atoms matching a SMARTS pattern as protected.

    This is the most convenient entry point when you want to lock a chemotype
    (e.g. a phenyl ring, an amide, a CF₃ group) without knowing the atom
    indices in advance.

    Parameters
    ----------
    mol    : RDKit Mol
    smarts : valid SMARTS string; raises ``ValueError`` if it cannot be parsed

    Returns a new RDKit Mol.
    """
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError(f"Invalid SMARTS: {smarts!r}")
    indices = [idx for match in mol.GetSubstructMatches(pattern) for idx in match]
    return protect_atoms_for_idx(mol, indices)


def unprotect_atoms_for_idx(mol, indices):
    """Remove protection from specific atoms.

    Only the listed atoms are affected; other protected atoms keep their mark.

    Parameters
    ----------
    mol     : RDKit Mol
    indices : iterable of atom indices to unprotect

    Returns a new RDKit Mol.
    """
    rwmol = Chem.RWMol(mol)
    for idx in indices:
        atom = rwmol.GetAtomWithIdx(idx)
        if atom.HasProp(PROP):
            atom.ClearProp(PROP)
    return rwmol.GetMol()


def unprotect_atoms_all(mol):
    """Remove protection from all atoms.

    Parameters
    ----------
    mol : RDKit Mol (not modified in place)

    Returns a new RDKit Mol with no atom protection.
    """
    rwmol = Chem.RWMol(mol)
    for atom in rwmol.GetAtoms():
        if atom.HasProp(PROP):
            atom.ClearProp(PROP)
    return rwmol.GetMol()


def protect_bonds_for_idx(mol, bond_indices):
    """Mark specific bonds as protected.

    Protected bonds are excluded from the LACAN score (see
    :func:`score_mol_ignoring_protected_bonds`).  They do **not** block
    mutations — atom protection is used for that.

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
# Querying protection
# ---------------------------------------------------------------------------

def get_protected_atom_indices(mol):
    """Return a list of indices of all protected atoms in *mol*."""
    return [a.GetIdx() for a in mol.GetAtoms()
            if a.HasProp(PROP) and a.GetBoolProp(PROP)]


def get_protected_bond_indices(mol):
    """Return a list of indices of all protected bonds in *mol*."""
    return [b.GetIdx() for b in mol.GetBonds()
            if b.HasProp(PROP) and b.GetBoolProp(PROP)]


def atom_is_protected(atom):
    """Return True if the RDKit Atom object has the ``_lp`` protection mark."""
    return atom.HasProp(PROP) and atom.GetBoolProp(PROP)


def bond_is_protected(bond):
    """Return True if the RDKit Bond object has the ``_lp`` protection mark."""
    return bond.HasProp(PROP) and bond.GetBoolProp(PROP)


# ---------------------------------------------------------------------------
# Reaction safety check
# ---------------------------------------------------------------------------

def reaction_touches_protected(mol, rxn):
    """Return True if any match of *rxn*'s reactant templates overlaps a protected atom.

    This is called by :func:`~lacan.mutate.apply_mutations` before each
    reaction to skip operations that would modify a protected atom.  The check
    iterates over all reactant templates and all their substructure matches; if
    any match atom index is in the protected set the reaction is skipped.

    Parameters
    ----------
    mol : RDKit Mol (may have protected atoms)
    rxn : RDKit ChemicalReaction

    Returns True (skip this reaction) or False (safe to apply).
    """
    protected = set(get_protected_atom_indices(mol))
    if not protected:
        return False
    for i in range(rxn.GetNumReactantTemplates()):
        template = rxn.GetReactantTemplate(i)
        for match in mol.GetSubstructMatches(template):
            if set(match) & protected:
                return True
    return False


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
        apb_active = [1.0]  # all bonds protected — trivially passes
    info = {}
    info["bad_bonds"] = [i for i, b in enumerate(apb) if b < t and i not in protected_bonds]
    if mode == "threshold":
        score = 0 if min(apb_active) < t else 1
    else:
        score = min(0.5 * (min(apb_active) / t) ** 0.5, 1.0)
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

    Protected atoms are still respected: any reaction whose template matches a
    protected atom is skipped.

    Degenerate products are rejected:
    * Disconnected molecules (contain ``"."`` in their canonical SMILES) — these
      arise from ring-opening and atom-deletion reactions and would collapse the
      search into tiny fragments.
    * Molecules with zero bonds — these score trivially as perfect (no bond
      violations) and would be incorrectly accepted.
    * Molecules with fewer than 4 heavy atoms.

    Parameters
    ----------
    mol     : RDKit Mol (typically has protected bonds from :func:`mol_cleaner`)
    profile : LACAN profile dict (unused here, kept for API consistency)

    Returns a list of sanitized, connected, non-degenerate RDKit Mol objects.
    """
    from rdkit.Chem import rdChemReactions
    from lacan.mutate import mutate_ops

    products = []
    for rxn in mutate_ops.values():
        if reaction_touches_protected(mol, rxn):
            continue
        for prod_tuple in rxn.RunReactants((mol,)):
            try:
                m = prod_tuple[0]
                Chem.SanitizeMol(m)
                # Reject disconnected, empty, or degenerate products
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
         *same* violation count.  This is essential for molecules with multiple
         independent violations: a lateral move on one region can reposition
         bad bonds to a state where the *next* mutation fixes one of them.

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

    # Minimum heavy-atom count: never accept a product smaller than this.
    # Atom-deletion mutations are legitimate but shouldn't collapse the molecule.
    min_heavy = max(4, mol.GetNumHeavyAtoms() - 6)

    def _reprotect(m):
        """Freeze currently-passing bonds so only bad regions stay mutable."""
        apb = lacan.assess_per_bond(m, profile)
        good = [i for i, sc in enumerate(apb) if sc >= t]
        return protect_bonds_for_idx(m, good)

    def _eval(m):
        """Return (n_violations, score) for a molecule with bond protection set.
        Returns (999, 0.0) for degenerate molecules (no bonds)."""
        if m.GetNumBonds() == 0:
            return 999, 0.0
        _, info = score_mol_ignoring_protected_bonds(m, profile, t=t)
        s, _ = score_mol_ignoring_protected_bonds(m, profile, t=t)
        return len(info["bad_bonds"]), s

    # Initialise: freeze all currently-passing bonds
    working_mol = _reprotect(mol)
    current_v, current_s = _eval(working_mol)

    lateral_count = 0
    visited_smis = {Chem.MolToSmiles(working_mol)}

    for _ in range(max_iter):
        if current_v == 0 and current_s >= score_threshold:
            return working_mol

        # Generate raw products — NO score filtering so partials are kept
        candidates = _raw_mutations(working_mol, profile)
        if not candidates:
            break

        best_improve_v = current_v
        best_improve_s = -1.0
        best_improve_mol = None
        best_lateral_s = -1.0
        best_lateral_mol = None

        for c in candidates:
            # Skip molecules that have shrunk too much
            if c.GetNumHeavyAtoms() < min_heavy:
                continue
            smi = Chem.MolToSmiles(c)
            if smi in visited_smis:
                continue
            # Re-protect the candidate on its own bond scores
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
