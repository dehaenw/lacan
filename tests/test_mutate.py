"""
tests/test_mutate.py — Tests for lacan/mutate.py
"""
import pytest
from rdkit import Chem
from lacan.mutate import apply_mutations, mutate_ops, mutate_smarts


class TestMutateOps:
    def test_all_smarts_compile(self):
        """Every SMARTS in mutate_smarts must produce a valid reaction."""
        from rdkit.Chem import rdChemReactions
        for name, smarts in mutate_smarts.items():
            rxn = rdChemReactions.ReactionFromSmarts(smarts)
            assert rxn is not None, f"SMARTS for {name!r} failed to compile"

    def test_mutate_ops_keys_match_smarts(self):
        assert set(mutate_ops.keys()) == set(mutate_smarts.keys())


class TestApplyMutations:
    def test_returns_list(self, toluene, chembl_profile):
        result = apply_mutations(toluene, chembl_profile, score_threshold=0.0)
        assert isinstance(result, list)

    def test_all_results_are_mols(self, toluene, chembl_profile):
        result = apply_mutations(toluene, chembl_profile, score_threshold=0.0)
        assert all(isinstance(m, Chem.Mol) for m in result)

    def test_all_results_are_valid(self, toluene, chembl_profile):
        result = apply_mutations(toluene, chembl_profile, score_threshold=0.0)
        for m in result:
            smi = Chem.MolToSmiles(m)
            assert Chem.MolFromSmiles(smi) is not None

    def test_nonempty_for_drug(self, ibuprofen, chembl_profile):
        result = apply_mutations(ibuprofen, chembl_profile, score_threshold=0.0)
        assert len(result) > 0

    def test_strict_threshold_reduces_results(self, ibuprofen, chembl_profile):
        loose = apply_mutations(ibuprofen, chembl_profile, score_threshold=0.0)
        strict = apply_mutations(ibuprofen, chembl_profile, score_threshold=0.99)
        assert len(strict) <= len(loose)

    def test_no_duplicates_by_inchikey(self, toluene, chembl_profile):
        from rdkit.Chem.inchi import MolToInchiKey
        result = apply_mutations(toluene, chembl_profile, score_threshold=0.0)
        inchikeys = [MolToInchiKey(m) for m in result]
        assert len(inchikeys) == len(set(inchikeys)), "Duplicate molecules in output"

    def test_mode_all(self, toluene, chembl_profile):
        result = apply_mutations(toluene, chembl_profile, score_threshold=0.0, mode="all")
        assert len(result) > 0

    def test_add_methyl(self):
        """addC mutation should add CH3 to a molecule with H."""
        from rdkit.Chem import rdChemReactions
        mol = Chem.MolFromSmiles("c1ccccc1")
        rxn = mutate_ops["addC"]
        prods = rxn.RunReactants((mol,))
        assert len(prods) > 0

    def test_aroc_to_n(self):
        """aroCtoN should replace aromatic CH with N in a 6-ring."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        rxn = mutate_ops["aroCtoN"]
        prods = rxn.RunReactants((mol,))
        assert len(prods) > 0

    # Regression: mode="random" bug
    # The original code assigned `ops` in random mode but still iterated over
    # `mutate_ops`. This test checks that random mode produces fewer ops than all.
    @pytest.mark.xfail(reason="Known bug: random mode iterates all ops, not just one")
    def test_random_mode_applies_single_op(self, toluene, chembl_profile):
        # In a fixed seed, "random" should produce a subset of "all" outputs
        import random
        random.seed(42)
        result_all = apply_mutations(toluene, chembl_profile, score_threshold=0.0, mode="all")
        random.seed(42)
        result_random = apply_mutations(toluene, chembl_profile, score_threshold=0.0, mode="random")
        assert len(result_random) < len(result_all)
