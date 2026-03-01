"""
tests/conftest.py — Shared pytest fixtures for LACAN tests
"""
import pytest
from rdkit import Chem
from lacan.lacan import load_profile


# ── Molecules ───────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def fluoxetine():
    return Chem.MolFromSmiles("CNCCC(c1ccccc1)Oc1ccc(C(F)(F)F)cc1")

@pytest.fixture(scope="session")
def aspirin():
    # Note: aspirin scores 0 against the chembl profile because it contains
    # a structural alert that was used to filter the ChEMBL training set.
    # It's kept as a fixture for tests that specifically want a low-scoring mol.
    return Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")

@pytest.fixture(scope="session")
def caffeine():
    return Chem.MolFromSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")

@pytest.fixture(scope="session")
def toluene():
    return Chem.MolFromSmiles("Cc1ccccc1")

@pytest.fixture(scope="session")
def diphenylmethane():
    return Chem.MolFromSmiles("c1ccc(Cc2ccccc2)cc1")

@pytest.fixture(scope="session")
def ibuprofen():
    return Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")

@pytest.fixture(scope="session")
def sildenafil():
    return Chem.MolFromSmiles(
        "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12"
    )


# ── Profile ─────────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def chembl_profile():
    return load_profile("chembl")


# ── Drug-like SMILES list ────────────────────────────────────────────────────

DRUG_SMILES = [
    "CNCCC(c1ccccc1)Oc1ccc(C(F)(F)F)cc1",  # fluoxetine
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",          # caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",          # ibuprofen
    "c1ccc2ccccc2c1",                       # naphthalene
]

@pytest.fixture
def drug_mols():
    return [Chem.MolFromSmiles(smi) for smi in DRUG_SMILES]
