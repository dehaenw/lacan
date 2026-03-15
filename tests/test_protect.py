"""
tests/test_protect.py — Tests for lacan/protect.py
"""
import pytest
from rdkit import Chem
from lacan.protect import (
    get_protected_atoms,
    protect_bonds_for_idx,
    protect_rejected_bonds,
    get_protected_bond_indices,
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
    return Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")


# ── get_protected_atoms ───────────────────────────────────────────────────────

class TestGetProtectedAtoms:
    def test_phenyl_ring(self, toluene):
        ring_atoms = frozenset(a.GetIdx() for a in toluene.GetAtoms() if a.IsInRing())
        assert get_protected_atoms(toluene, "c1ccccc1") == ring_atoms

    def test_no_smarts_returns_empty(self, toluene):
        assert get_protected_atoms(toluene, None) == frozenset()

    def test_empty_smarts_returns_empty(self, toluene):
        assert get_protected_atoms(toluene, "") == frozenset()

    def test_invalid_smarts_raises(self, toluene):
        with pytest.raises(ValueError):
            get_protected_atoms(toluene, "!!notvalid!!")

    def test_no_match_returns_empty(self, toluene):
        assert get_protected_atoms(toluene, "[r3]") == frozenset()

    def test_multiple_matches_all_covered(self):
        # Naphthalene has two rings — all aromatic atoms should be found
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        pa = get_protected_atoms(mol, "c1ccccc1")
        assert len(pa) == mol.GetNumAtoms()


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
        assert len(get_protected_bond_indices(mol)) > 0

    def test_fluoxetine_no_rejected_bonds(self, fluoxetine, profile):
        mol = protect_rejected_bonds(fluoxetine, profile)
        assert get_protected_bond_indices(mol) == []

    def test_returns_new_mol(self, aspirin, profile):
        mol = protect_rejected_bonds(aspirin, profile)
        assert isinstance(mol, Chem.Mol)


# ── reaction_touches_protected ────────────────────────────────────────────────

class TestReactionTouchesProtected:
    def test_no_protection_always_false(self, toluene):
        from rdkit.Chem import rdChemReactions
        rxn = rdChemReactions.ReactionFromSmarts("[cH:0]>>[nH0:0]")
        assert not reaction_touches_protected(toluene, rxn, None)

    def test_protected_reactive_atom_true(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        from rdkit.Chem import rdChemReactions
        rxn = rdChemReactions.ReactionFromSmarts("[cH:0]>>[nH0:0]")
        # Protect one aromatic C via SMARTS
        assert reaction_touches_protected(mol, rxn, "c1ccccc1")

    def test_protected_nonreactive_atom_false(self, toluene):
        # Protect only the methyl carbon; reaction targets aromatic C
        from rdkit.Chem import rdChemReactions
        rxn = rdChemReactions.ReactionFromSmarts("[cH:0]>>[nH0:0]")
        assert not reaction_touches_protected(toluene, rxn, "[CH3]")


# ── score_mol_ignoring_protected_bonds ────────────────────────────────────────

class TestScoreMolIgnoringProtectedBonds:
    def test_aspirin_passes_when_bad_bonds_protected(self, aspirin, profile):
        mol = protect_rejected_bonds(aspirin, profile)
        score, info = score_mol_ignoring_protected_bonds(mol, profile)
        assert score > 0.0

    def test_protected_bonds_not_in_bad_bonds(self, aspirin, profile):
        mol = protect_rejected_bonds(aspirin, profile)
        _, info = score_mol_ignoring_protected_bonds(mol, profile)
        protected = set(get_protected_bond_indices(mol))
        assert not (set(info["bad_bonds"]) & protected)

    def test_unprotected_mol_matches_regular_score(self, fluoxetine, profile):
        from lacan.lacan import score_mol
        score_regular, _ = score_mol(fluoxetine, profile)
        score_protect, _ = score_mol_ignoring_protected_bonds(fluoxetine, profile)
        assert abs(score_regular - score_protect) < 1e-6

    def test_fully_protected_mol_scores_1(self, toluene, profile):
        mol = protect_bonds_for_idx(toluene, list(range(toluene.GetNumBonds())))
        score, _ = score_mol_ignoring_protected_bonds(mol, profile)
        assert score == 1.0


# ── mol_cleaner ───────────────────────────────────────────────────────────────

class TestMolCleaner:
    def test_fluoxetine_already_clean(self, fluoxetine, profile):
        result = mol_cleaner(fluoxetine, profile, score_threshold=0.5)
        assert result is not None

    def test_returns_mol_or_none(self, aspirin, profile):
        result = mol_cleaner(aspirin, profile, score_threshold=0.5, max_iter=5)
        assert result is None or isinstance(result, Chem.Mol)

    def test_result_is_valid_mol(self, fluoxetine, profile):
        result = mol_cleaner(fluoxetine, profile)
        assert Chem.MolFromSmiles(Chem.MolToSmiles(result)) is not None


class TestMolCleanerExtended:
    def test_partial_solutions_not_rejected(self, aspirin, profile):
        from lacan.protect import _raw_mutations
        from lacan import lacan
        apb = lacan.assess_per_bond(aspirin, profile)
        good = [i for i, sc in enumerate(apb) if sc >= 0.05]
        working = protect_bonds_for_idx(aspirin, good)
        candidates = _raw_mutations(working, profile)
        assert len(candidates) > 0

    def test_cleaner_makes_progress(self, aspirin, profile):
        result = mol_cleaner(aspirin, profile, score_threshold=0.5,
                             max_iter=50, lateral_patience=8)
        if result is not None:
            assert Chem.MolFromSmiles(Chem.MolToSmiles(result)) is not None
            final_score, _ = score_mol_ignoring_protected_bonds(result, profile)
            assert final_score >= 0.5

    def test_two_violation_mol(self, profile):
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        assert mol is not None
        result = mol_cleaner(mol, profile, score_threshold=0.5,
                             max_iter=80, lateral_patience=10)
        if result is not None:
            score, _ = score_mol_ignoring_protected_bonds(result, profile)
            assert score >= 0.5

    def test_already_clean_returns_fast(self, fluoxetine, profile):
        result = mol_cleaner(fluoxetine, profile, score_threshold=0.5, max_iter=100)
        assert result is not None
        assert Chem.MolFromSmiles(Chem.MolToSmiles(result)) is not None
