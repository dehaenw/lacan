"""
tests/test_decompose.py — Tests for lacan/decompose.py
"""
import pytest
from rdkit import Chem
from lacan.decompose import decompose_molecule, get_bonds_string, get_corpus


class TestDecomposeMolecule:
    def test_returns_three_lists(self):
        rings, linkers, subs = decompose_molecule("Cc1ccccc1")
        assert isinstance(rings, list)
        assert isinstance(linkers, list)
        assert isinstance(subs, list)

    def test_toluene_has_ring_and_sub(self):
        rings, linkers, subs = decompose_molecule("Cc1ccccc1")
        assert len(rings) >= 1
        assert len(subs) >= 1
        assert len(linkers) == 0

    def test_fragments_contain_dummy_atoms(self):
        rings, linkers, subs = decompose_molecule("Cc1ccccc1")
        for frag in rings + linkers + subs:
            assert "*" in frag, f"Fragment {frag!r} should contain dummy atom"

    def test_benzene_no_substituents(self):
        rings, linkers, subs = decompose_molecule("c1ccccc1")
        # Benzene has no exocyclic bonds to break
        assert linkers == [] and subs == []

    def test_diphenylmethane_has_linker(self):
        rings, linkers, subs = decompose_molecule("c1ccc(Cc2ccccc2)cc1")
        assert len(rings) >= 2
        assert len(linkers) >= 1

    def test_invalid_smiles_returns_empty(self):
        rings, linkers, subs = decompose_molecule("not_valid_smiles")
        assert rings == [] and linkers == [] and subs == []

    def test_fragments_are_valid_smiles(self):
        rings, linkers, subs = decompose_molecule("Cc1ccccc1C(=O)O")
        for frag in rings + linkers + subs:
            mol = Chem.MolFromSmiles(frag)
            assert mol is not None, f"Fragment {frag!r} is not valid SMILES"

    def test_asmol_false_accepts_mol(self):
        mol = Chem.MolFromSmiles("Cc1ccccc1")
        rings, linkers, subs = decompose_molecule(mol, asSmi=False)
        assert len(rings) >= 1

    def test_multi_ring_molecule(self):
        # Ibuprofen: two rings, several subs
        rings, linkers, subs = decompose_molecule("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        assert len(rings) >= 1

    def test_reproducible(self):
        smi = "Cc1ccccc1"
        r1, l1, s1 = decompose_molecule(smi)
        r2, l2, s2 = decompose_molecule(smi)
        assert sorted(r1) == sorted(r2)
        assert sorted(l1) == sorted(l2)
        assert sorted(s1) == sorted(s2)


class TestGetBondsString:
    def test_single_bond_sub(self):
        bs = get_bonds_string("[*]-C")
        assert bs == "-"

    def test_double_bond_sub(self):
        bs = get_bonds_string("[*]=C")
        assert bs == "="

    def test_two_single_bonds(self):
        bs = get_bonds_string("[*]-C-[*]")
        assert bs == "--"

    def test_mixed_bonds(self):
        bs = get_bonds_string("[*]-C=[*]")
        assert "-" in bs and "=" in bs


class TestGetCorpus:
    def test_returns_list(self):
        mols = [Chem.MolFromSmiles(smi) for smi in ["Cc1ccccc1", "c1ccc(Cc2ccccc2)cc1"]]
        corpus = get_corpus(mols)
        assert isinstance(corpus, list)

    def test_entries_have_5_fields(self):
        mols = [Chem.MolFromSmiles("Cc1ccccc1")]
        corpus = get_corpus(mols)
        assert all(len(e) == 5 for e in corpus)

    def test_sorted_by_occurrence_desc(self):
        mols = [Chem.MolFromSmiles(smi) for smi in ["Cc1ccccc1"] * 10]
        corpus = get_corpus(mols)
        counts = [e[1] for e in corpus]
        assert counts == sorted(counts, reverse=True)
