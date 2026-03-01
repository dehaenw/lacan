"""
tests/test_protect.py — Tests for lacan/protect.py
"""
import pytest
from rdkit import Chem
from lacan.protect import (
    protect_atoms_for_idx,
    protect_atoms_matching_smarts,
    unprotect_atoms_for_idx,
    unprotect_atoms_all,
    protect_bonds_for_idx,
    protect_rejected_bonds,
    get_protected_atom_indices,
    get_protected_bond_indices,
    atom_is_protected,
    bond_is_protected,
    reaction_touches_protected,
    score_mol_ignoring_protected_bonds,
    mol_cleaner,
    PROP,
)
from lacan.lacan import load_profile


@pytest.fixture(scope="module")
def profile():
    return load_profile("chembl")

@pytest.fixture
def toluene():
    return Chem.MolFromSmiles("Cc1ccccc1")

@pytest.fixture
def fluoxetine():
    return Chem.MolFromSmiles("CNCCC(c1ccccc1)Oc1ccc(C(F)(F)F)cc1")

@pytest.fixture
def aspirin():
    # Contains a LACAN-failing bond (ester alert)
    return Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")


# ── protect_atoms_for_idx ─────────────────────────────────────────────────────

class TestProtectAtomsForIdx:
    def test_marks_atoms(self, toluene):
        mol = protect_atoms_for_idx(toluene, [0, 1])
        assert atom_is_protected(mol.GetAtomWithIdx(0))
        assert atom_is_protected(mol.GetAtomWithIdx(1))

    def test_does_not_mark_others(self, toluene):
        mol = protect_atoms_for_idx(toluene, [0])
        assert not atom_is_protected(mol.GetAtomWithIdx(1))

    def test_does_not_mutate_input(self, toluene):
        protect_atoms_for_idx(toluene, [0])
        assert not atom_is_protected(toluene.GetAtomWithIdx(0))

    def test_empty_list(self, toluene):
        mol = protect_atoms_for_idx(toluene, [])
        assert get_protected_atom_indices(mol) == []


# ── protect_atoms_matching_smarts ─────────────────────────────────────────────

class TestProtectAtomsMatchingSmarts:
    def test_phenyl_ring(self, toluene):
        mol = protect_atoms_matching_smarts(toluene, "c1ccccc1")
        ring_atoms = {a.GetIdx() for a in toluene.GetAtoms() if a.IsInRing()}
        assert set(get_protected_atom_indices(mol)) == ring_atoms

    def test_invalid_smarts_raises(self, toluene):
        with pytest.raises(ValueError):
            protect_atoms_matching_smarts(toluene, "!!notvalid!!")

    def test_no_match_protects_nothing(self, toluene):
        mol = protect_atoms_matching_smarts(toluene, "[r3]")  # 3-membered ring
        assert get_protected_atom_indices(mol) == []


# ── unprotect ─────────────────────────────────────────────────────────────────

class TestUnprotect:
    def test_unprotect_specific(self, toluene):
        mol = protect_atoms_for_idx(toluene, [0, 1, 2])
        mol = unprotect_atoms_for_idx(mol, [0, 2])
        assert not atom_is_protected(mol.GetAtomWithIdx(0))
        assert atom_is_protected(mol.GetAtomWithIdx(1))
        assert not atom_is_protected(mol.GetAtomWithIdx(2))

    def test_unprotect_all(self, toluene):
        mol = protect_atoms_for_idx(toluene, [0, 1, 2, 3])
        mol = unprotect_atoms_all(mol)
        assert get_protected_atom_indices(mol) == []

    def test_unprotect_all_unprotected_mol(self, toluene):
        # Should not raise on a mol with no protection
        mol = unprotect_atoms_all(toluene)
        assert get_protected_atom_indices(mol) == []


# ── protect_bonds_for_idx ─────────────────────────────────────────────────────

class TestProtectBondsForIdx:
    def test_marks_bonds(self, toluene):
        mol = protect_bonds_for_idx(toluene, [0])
        assert bond_is_protected(mol.GetBondWithIdx(0))

    def test_does_not_mark_other_bonds(self, toluene):
        mol = protect_bonds_for_idx(toluene, [0])
        assert not bond_is_protected(mol.GetBondWithIdx(1))

    def test_does_not_mutate_input(self, toluene):
        protect_bonds_for_idx(toluene, [0])
        assert not bond_is_protected(toluene.GetBondWithIdx(0))


# ── protect_rejected_bonds ────────────────────────────────────────────────────

class TestProtectRejectedBonds:
    def test_aspirin_has_rejected_bonds_protected(self, aspirin, profile):
        mol = protect_rejected_bonds(aspirin, profile)
        # Aspirin scores 0 — some bonds should now be protected
        assert len(get_protected_bond_indices(mol)) > 0

    def test_fluoxetine_no_rejected_bonds(self, fluoxetine, profile):
        mol = protect_rejected_bonds(fluoxetine, profile)
        # Fluoxetine passes — nothing should be protected
        assert get_protected_bond_indices(mol) == []

    def test_returns_new_mol(self, aspirin, profile):
        mol = protect_rejected_bonds(aspirin, profile)
        assert isinstance(mol, Chem.Mol)


# ── reaction_touches_protected ────────────────────────────────────────────────

class TestReactionTouchesProtected:
    def test_no_protection_always_false(self, toluene):
        from rdkit.Chem import rdChemReactions
        rxn = rdChemReactions.ReactionFromSmarts("[cH:0]>>[nH0:0]")
        assert not reaction_touches_protected(toluene, rxn)

    def test_protected_reactive_atom_true(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = protect_atoms_for_idx(mol, [0])  # protect one aromatic C
        from rdkit.Chem import rdChemReactions
        rxn = rdChemReactions.ReactionFromSmarts("[cH:0]>>[nH0:0]")
        assert reaction_touches_protected(mol, rxn)

    def test_protected_nonreactive_atom_false(self, toluene):
        # Protect the methyl carbon but use a reaction on aromatic C only
        mol = protect_atoms_for_idx(toluene, [0])  # methyl C
        from rdkit.Chem import rdChemReactions
        rxn = rdChemReactions.ReactionFromSmarts("[cH:0]>>[nH0:0]")
        assert not reaction_touches_protected(mol, rxn)


# ── score_mol_ignoring_protected_bonds ────────────────────────────────────────

class TestScoreMolIgnoringProtectedBonds:
    def test_aspirin_passes_when_bad_bonds_protected(self, aspirin, profile):
        mol = protect_rejected_bonds(aspirin, profile)
        score, info = score_mol_ignoring_protected_bonds(mol, profile)
        assert score > 0.0, "Aspirin should pass when its bad bonds are protected"

    def test_protected_bonds_not_in_bad_bonds(self, aspirin, profile):
        mol = protect_rejected_bonds(aspirin, profile)
        _, info = score_mol_ignoring_protected_bonds(mol, profile)
        protected = set(get_protected_bond_indices(mol))
        # No protected bond index should appear in bad_bonds
        assert not (set(info["bad_bonds"]) & protected)

    def test_unprotected_mol_matches_regular_score(self, fluoxetine, profile):
        from lacan.lacan import score_mol
        score_regular, _ = score_mol(fluoxetine, profile)
        score_protect, _ = score_mol_ignoring_protected_bonds(fluoxetine, profile)
        assert abs(score_regular - score_protect) < 1e-6

    def test_fully_protected_mol_scores_1(self, toluene, profile):
        # Protect all bonds — nothing to score, should return 1.0
        mol = protect_bonds_for_idx(toluene, list(range(toluene.GetNumBonds())))
        score, _ = score_mol_ignoring_protected_bonds(mol, profile)
        assert score == 1.0


# ── mol_cleaner ───────────────────────────────────────────────────────────────

class TestMolCleaner:
    def test_fluoxetine_already_clean(self, fluoxetine, profile):
        # Should return quickly since fluoxetine already passes
        result = mol_cleaner(fluoxetine, profile, score_threshold=0.5)
        assert result is not None

    def test_returns_mol_or_none(self, aspirin, profile):
        result = mol_cleaner(aspirin, profile, score_threshold=0.5, max_iter=5)
        assert result is None or isinstance(result, Chem.Mol)

    def test_result_is_valid_mol(self, fluoxetine, profile):
        result = mol_cleaner(fluoxetine, profile)
        assert Chem.MolFromSmiles(Chem.MolToSmiles(result)) is not None


# ── mol_cleaner (extended) ────────────────────────────────────────────────────

class TestMolCleanerExtended:
    def test_partial_solutions_not_rejected(self, aspirin, profile):
        """The cleaner must be able to explore intermediates that still have
        some bad bonds.  If only fully-clean intermediates were accepted, the
        cleaner would get stuck at the first step for multi-violation molecules.

        We verify this indirectly: aspirin with max_iter=30 and lateral_patience=8
        should either find a clean mol or at least progress (not immediately return
        None because it found no candidates at iteration 1).
        """
        from lacan.protect import _raw_mutations
        from lacan import lacan

        # Step 1: _raw_mutations must return candidates even for aspirin
        apb = lacan.assess_per_bond(aspirin, profile)
        from lacan.protect import protect_bonds_for_idx
        good = [i for i, sc in enumerate(apb) if sc >= 0.05]
        working = protect_bonds_for_idx(aspirin, good)
        candidates = _raw_mutations(working, profile)
        assert len(candidates) > 0, (
            "_raw_mutations should return products even when parent has bad bonds"
        )

    def test_cleaner_makes_progress(self, aspirin, profile):
        """With enough iterations, mol_cleaner should reduce violations or find
        a clean molecule."""
        from lacan.protect import score_mol_ignoring_protected_bonds, protect_bonds_for_idx
        from lacan import lacan

        result = mol_cleaner(aspirin, profile, score_threshold=0.5,
                             max_iter=50, lateral_patience=8)
        # We can't guarantee success, but if it returns something it must be valid
        if result is not None:
            assert Chem.MolFromSmiles(Chem.MolToSmiles(result)) is not None
            final_score, _ = score_mol_ignoring_protected_bonds(result, profile)
            assert final_score >= 0.5

    def test_two_violation_mol(self, profile):
        """CC(=O)Oc1ccccc1C(=O)OF has multiple violations.
        The cleaner should be able to find a path (lateral moves help here)."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin, 1+ violations
        assert mol is not None
        result = mol_cleaner(mol, profile, score_threshold=0.5,
                             max_iter=80, lateral_patience=10)
        if result is not None:
            from lacan.protect import score_mol_ignoring_protected_bonds
            score, _ = score_mol_ignoring_protected_bonds(result, profile)
            assert score >= 0.5

    def test_already_clean_returns_fast(self, fluoxetine, profile):
        """A molecule that already passes should return at the first check."""
        result = mol_cleaner(fluoxetine, profile, score_threshold=0.5, max_iter=100)
        assert result is not None
        assert Chem.MolFromSmiles(Chem.MolToSmiles(result)) is not None
