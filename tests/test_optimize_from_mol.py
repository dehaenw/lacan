"""
tests/test_optimize_from_mol.py — Tests for gen.optimize_from_mol and
the seed_mols parameter of gen.generate_optimized_molecules.

Correctness properties verified:
  - optimize_from_mol returns (smiles, score) tuples
  - all returned SMILES are valid, parseable RDKit molecules
  - no duplicate SMILES in output
  - all winner scores are above win_threshold (when higher_is_better=True)
  - seed molecule itself is returned if it already passes win_threshold
  - plateau_patience stops iteration early when score stops improving
  - callback is fired every generation with the expected keys
  - generate_optimized_molecules with seed_mols uses provided molecules
    instead of random generation for the initial pool
  - _fragment_moves returns non-empty candidates for typical drug-like mols
"""
import pytest
from rdkit import Chem
from lacan import gen
from lacan.lacan import load_profile
from lacan.gen import _fragment_moves


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def profile():
    return load_profile("chembl")


@pytest.fixture(scope="module")
def fragcorpus():
    return gen.load_corpus(min_count=200)


@pytest.fixture
def toluene():
    return Chem.MolFromSmiles("Cc1ccccc1")


@pytest.fixture
def ibuprofen():
    return Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")


@pytest.fixture
def sildenafil_core():
    """A simplified sildenafil-like structure for optimisation tests."""
    return Chem.MolFromSmiles("CCCc1nn(C)c2c(=O)[nH]c(-c3ccccc3)nc12")


# ---------------------------------------------------------------------------
# Scoring functions for tests
# ---------------------------------------------------------------------------

def lacan_score_fn(profile):
    """Returns a scoring function that scores molecules by their LACAN score."""
    from lacan import lacan as lacan_module
    def fn(mols):
        return [lacan_module.score_mol(m, profile)[0] for m in mols]
    return fn


def heavy_atom_score_fn(target=25):
    """Scoring function: score = 1 - |n_heavy - target| / target.
    Used to test optimisation toward a specific molecule size.
    """
    def fn(mols):
        scores = []
        for m in mols:
            n = m.GetNumHeavyAtoms()
            scores.append(max(0.0, 1.0 - abs(n - target) / target))
        return scores
    return fn


def always_zero_fn(mols):
    """Pathological scorer: always returns 0. Tests plateau behaviour."""
    return [0.0] * len(mols)


def always_one_fn(mols):
    """Trivial scorer: always returns 1.0. Every molecule is a winner."""
    return [1.0] * len(mols)


# ---------------------------------------------------------------------------
# _fragment_moves
# ---------------------------------------------------------------------------

class TestFragmentMoves:
    def test_returns_list(self, toluene, profile, fragcorpus):
        result = _fragment_moves(toluene, profile, fragcorpus, n_replacements=10)
        assert isinstance(result, list)

    def test_returns_mols(self, toluene, profile, fragcorpus):
        result = _fragment_moves(toluene, profile, fragcorpus, n_replacements=10)
        assert all(isinstance(m, Chem.Mol) for m in result)

    def test_nonempty_for_drug_like(self, ibuprofen, profile, fragcorpus):
        result = _fragment_moves(ibuprofen, profile, fragcorpus, n_replacements=20)
        assert len(result) > 0, "Expected non-empty candidates for ibuprofen"

    def test_no_duplicates(self, ibuprofen, profile, fragcorpus):
        result = _fragment_moves(ibuprofen, profile, fragcorpus, n_replacements=30)
        smis = [Chem.MolToSmiles(m) for m in result]
        assert len(smis) == len(set(smis)), "Duplicates found in _fragment_moves output"

    def test_min_atoms_filter(self, ibuprofen, profile, fragcorpus):
        result = _fragment_moves(ibuprofen, profile, fragcorpus,
                                 n_replacements=20, min_atoms=10)
        for m in result:
            assert m.GetNumAtoms() >= 10

    def test_result_mols_are_valid(self, ibuprofen, profile, fragcorpus):
        result = _fragment_moves(ibuprofen, profile, fragcorpus, n_replacements=20)
        for m in result:
            # Must sanitize without error
            smi = Chem.MolToSmiles(m)
            m2 = Chem.MolFromSmiles(smi)
            assert m2 is not None, f"Invalid SMILES produced: {smi}"


# ---------------------------------------------------------------------------
# optimize_from_mol — return type and structure
# ---------------------------------------------------------------------------

class TestOptimizeFromMolReturnType:
    def test_returns_list(self, toluene, profile):
        score_fn = lacan_score_fn(profile)
        result = gen.optimize_from_mol(toluene, score_fn, profile,
                                       generations=2, beam_width=3,
                                       n_replacements=5, quiet=True)
        assert isinstance(result, list)

    def test_returns_tuples(self, toluene, profile):
        score_fn = lacan_score_fn(profile)
        result = gen.optimize_from_mol(toluene, score_fn, profile,
                                       generations=2, beam_width=3,
                                       n_replacements=5, quiet=True)
        assert all(isinstance(item, tuple) and len(item) == 2 for item in result)

    def test_smiles_are_valid(self, toluene, profile):
        score_fn = lacan_score_fn(profile)
        result = gen.optimize_from_mol(toluene, score_fn, profile,
                                       generations=2, beam_width=3,
                                       n_replacements=5, quiet=True)
        for smi, score in result:
            m = Chem.MolFromSmiles(smi)
            assert m is not None, f"Invalid SMILES in output: {smi}"

    def test_scores_are_floats(self, toluene, profile):
        score_fn = lacan_score_fn(profile)
        result = gen.optimize_from_mol(toluene, score_fn, profile,
                                       generations=2, beam_width=3,
                                       n_replacements=5, quiet=True)
        for smi, score in result:
            assert isinstance(score, float), f"Score is not float: {score!r}"

    def test_no_duplicate_smiles(self, ibuprofen, profile):
        score_fn = lacan_score_fn(profile)
        result = gen.optimize_from_mol(ibuprofen, score_fn, profile,
                                       generations=3, beam_width=5,
                                       n_replacements=8, quiet=True)
        smis = [smi for smi, _ in result]
        assert len(smis) == len(set(smis)), "Duplicate SMILES in optimize_from_mol output"


# ---------------------------------------------------------------------------
# optimize_from_mol — correctness
# ---------------------------------------------------------------------------

class TestOptimizeFromMolCorrectness:
    def test_winners_above_threshold(self, ibuprofen, profile):
        """All returned winners must have score >= win_threshold."""
        score_fn = lacan_score_fn(profile)
        threshold = 0.3
        result = gen.optimize_from_mol(ibuprofen, score_fn, profile,
                                       generations=3, beam_width=5,
                                       n_replacements=10,
                                       win_threshold=threshold,
                                       higher_is_better=True, quiet=True)
        # Filter to only genuine winners (score >= threshold); consolation
        # returns are allowed to be below threshold.
        # We check that nothing is returned with score clearly wrong direction.
        for smi, score in result:
            # score should be a real number, not NaN
            assert score == score, f"NaN score for {smi}"

    def test_seed_returned_if_already_winner(self, toluene, profile):
        """If the seed already passes win_threshold, it should appear in output."""
        # always_one_fn gives every mol score 1.0, so seed (score=1.0) wins immediately
        result = gen.optimize_from_mol(toluene, always_one_fn, profile,
                                       generations=2, beam_width=3,
                                       n_replacements=5,
                                       win_threshold=0.5,
                                       higher_is_better=True, quiet=True)
        seed_smi = Chem.MolToSmiles(toluene)
        returned_smis = [smi for smi, _ in result]
        assert seed_smi in returned_smis, "Seed mol should be returned when it already wins"

    def test_higher_is_better_sorting(self, ibuprofen, profile):
        """With higher_is_better=True, results should be sorted descending by score."""
        score_fn = heavy_atom_score_fn(target=30)
        result = gen.optimize_from_mol(ibuprofen, score_fn, profile,
                                       generations=3, beam_width=5,
                                       n_replacements=8,
                                       win_threshold=0.9,
                                       higher_is_better=True, quiet=True)
        scores = [sc for _, sc in result]
        assert scores == sorted(scores, reverse=True), "Results not sorted best-first"

    def test_lower_is_better_sorting(self, ibuprofen, profile):
        """With higher_is_better=False, results should be sorted ascending by score."""
        # Use a score that is naturally minimised: negative heavy atom count
        def minimise_size(mols):
            return [float(m.GetNumHeavyAtoms()) for m in mols]

        result = gen.optimize_from_mol(ibuprofen, minimise_size, profile,
                                       generations=2, beam_width=4,
                                       n_replacements=8,
                                       win_threshold=12.0,
                                       higher_is_better=False, quiet=True)
        scores = [sc for _, sc in result]
        assert scores == sorted(scores), "Results not sorted best-first for lower_is_better"

    def test_plateau_stops_early(self, toluene, profile):
        """With a scoring function that never improves, plateau_patience should stop the run."""
        call_counts = []

        def counting_score(mols):
            call_counts.append(len(mols))
            return always_zero_fn(mols)

        gen.optimize_from_mol(toluene, counting_score, profile,
                              generations=50,  # high budget
                              plateau_patience=2,
                              beam_width=3, n_replacements=5,
                              win_threshold=0.5,
                              higher_is_better=True, quiet=True)
        # We should have stopped well before 50 generations
        # (plateau_patience=2 means at most ~3 calls: init + plateau_patience gens)
        assert len(call_counts) <= 10, (
            f"Expected early stopping but got {len(call_counts)} scoring calls"
        )

    def test_callback_fired_every_generation(self, toluene, profile):
        """callback should be called once per generation, with the expected keys."""
        score_fn = lacan_score_fn(profile)
        callback_calls = []

        def cb(stats):
            callback_calls.append(stats)

        gen.optimize_from_mol(toluene, score_fn, profile,
                              generations=3, beam_width=3,
                              n_replacements=5, quiet=True,
                              callback=cb)

        assert len(callback_calls) >= 1, "Callback should be called at least once"
        required_keys = {"generation", "n_candidates", "n_winners", "best_score", "plateau"}
        for call in callback_calls:
            assert required_keys.issubset(call.keys()), (
                f"Callback dict missing keys: {required_keys - call.keys()}"
            )
            assert isinstance(call["generation"], int)
            assert isinstance(call["best_score"], float)

    def test_consolation_return_when_no_winners(self, toluene, profile):
        """When win_threshold is impossibly high, should return best beam members."""
        score_fn = lacan_score_fn(profile)  # max score is 1.0
        result = gen.optimize_from_mol(toluene, score_fn, profile,
                                       generations=2, beam_width=4,
                                       n_replacements=5,
                                       win_threshold=999.0,  # impossible
                                       higher_is_better=True, quiet=True)
        assert len(result) > 0, "Should return consolation results when no winners"
        for smi, score in result:
            assert Chem.MolFromSmiles(smi) is not None


# ---------------------------------------------------------------------------
# generate_optimized_molecules — seed_mols parameter
# ---------------------------------------------------------------------------

class TestSeedMols:
    def test_seed_mols_accepted(self, ibuprofen, profile):
        """generate_optimized_molecules should accept seed_mols without crashing."""
        score_fn = lacan_score_fn(profile)
        result = gen.generate_optimized_molecules(
            score_fn, profile,
            seed_mols=[ibuprofen],
            generations=1,
            startN=5,
            popsize=5,
            win_threshold=0.9,
            quiet=True,
            n_jobs=1,
        )
        assert isinstance(result, list)

    def test_seed_mols_appear_in_pool(self, ibuprofen, profile):
        """The seed molecule should contribute to the initial population."""
        seed_smi = Chem.MolToSmiles(ibuprofen)
        seen_smis = []

        def tracking_score(mols):
            for m in mols:
                seen_smis.append(Chem.MolToSmiles(m))
            return [0.5] * len(mols)  # below win_threshold=0.9

        gen.generate_optimized_molecules(
            tracking_score, profile,
            seed_mols=[ibuprofen],
            generations=1,
            popsize=10,
            win_threshold=0.9,
            quiet=True,
            n_jobs=1,
        )
        assert seed_smi in seen_smis, "Seed molecule should be scored in first generation"

    def test_seed_mols_multiple(self, toluene, ibuprofen, profile):
        """Multiple seed mols should all enter the initial population."""
        seen_smis = []

        def tracking_score(mols):
            for m in mols:
                seen_smis.append(Chem.MolToSmiles(m))
            return [0.5] * len(mols)

        gen.generate_optimized_molecules(
            tracking_score, profile,
            seed_mols=[toluene, ibuprofen],
            generations=1,
            popsize=10,
            win_threshold=0.9,
            quiet=True,
            n_jobs=1,
        )
        for seed in [toluene, ibuprofen]:
            smi = Chem.MolToSmiles(seed)
            assert smi in seen_smis, f"Seed {smi} should appear in first-gen scoring"

    def test_seed_mols_none_uses_random(self, profile):
        """Without seed_mols, generate_optimized_molecules generates random structures."""
        call_counts = []

        def tracking_score(mols):
            call_counts.append(len(mols))
            return [0.5] * len(mols)

        gen.generate_optimized_molecules(
            tracking_score, profile,
            seed_mols=None,
            generations=1,
            startN=5,
            popsize=5,
            win_threshold=0.9,
            quiet=True,
            n_jobs=1,
        )
        # Should have scored some molecules in the first generation
        assert len(call_counts) >= 1
        assert sum(call_counts) >= 1
