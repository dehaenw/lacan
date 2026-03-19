"""
lacan.py — Core scoring module for LACAN.

LACAN (Leveraging Adjacent Co-occurrence of Atomic Neighborhoods) scores molecules by asking how statistically likely each bond
is, given the chemical environments on both sides of it.

Scoring model
-------------
For each bond in a molecule, LACAN computes a pair of **ECFP2-like atom
environment identifiers** — one for each endpoint — and looks up their
co-occurrence statistics in a pre-built profile.

The atom environment identifier encodes five features:

* Atomic number
* Heavy-atom degree
* Total hydrogen count (explicit + implicit)
* Formal charge
* Ring type: 0 if acyclic, 1 if in a non-aromatic ring, 2 if in an aromatic ring

These five integers are hashed to a 32-bit integer via MD5
(:func:`hash_invariants`).  The two hashes for a bond are sorted and stored as
a tuple — the *bond pair* — so that the bond (A→B) and (B→A) map to the same
key.

A **profile** is a dict with three entries built from a molecule set:

``idx``
    ``{hash: count}`` — how often each atom environment appears.

``pairs``
    ``{(hash1, hash2): count}`` — how often each bond pair co-occurs.

``setsize``
    Total number of bonds seen, used to normalise counts to frequencies.

For a bond with environments *e1* and *e2* the **pointwise mutual information**
(PMI) score is::

    observed  = pairs[(e1, e2)] / setsize
    expected  = (idx[e1] / setsize / 2) * (idx[e2] / setsize / 2)
    PMI       = observed / expected

A PMI > 1 means the two environments co-occur more often than chance; a PMI
near 0 means the combination is chemically unusual.

The molecule-level score is derived from the **minimum** per-bond PMI::

    score = min_PMI / (1 + min_PMI)

This saturates toward 1.0 as the worst-bond PMI grows large, and approaches
0 when the worst bond is near zero.  Bonds below a threshold *t* (default
0.05) are reported as ``bad_bonds`` in :func:`score_mol`.

CLI usage
---------
::

    python -m lacan.lacan -i molecules.smi -p chembl -m score -t 0.05

Modes:

``score``
    Print per-molecule score and list of failing bond indices.

``profile``
    Build a new profile from the input SMILES file and save it to
    ``data/<profile_name>.pickle``.
"""

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from collections import Counter
import sys, os
import pickle
import argparse
import math
import hashlib
import multiprocessing
from itertools import chain
import copy

RDLogger.DisableLog('rdApp.*')
MFPGEN = rdFingerprintGenerator.GetMorganGenerator(1)
ao = rdFingerprintGenerator.AdditionalOutput()
ao.AllocateBitInfoMap()
ao.AllocateAtomToBits()
p = Chem.MolFromSmarts("[#0]")  # dummy atom
srflag = Chem.SanitizeFlags.SANITIZE_SYMMRINGS


def hash_invariants(invs):
    """Hash a list of integers to a reproducible, portable 32-bit signed int.

    Uses the first 4 bytes of an MD5 digest so that the hash is identical
    across Python versions and platforms — a requirement for profiles built on
    one machine to be usable on another.

    Parameters
    ----------
    invs : list of int
        Flattened atom-environment descriptor (atomic number, degree, H count,
        charge, ring type, bond type, …).

    Returns
    -------
    int
        32-bit signed integer hash.
    """
    h = hashlib.md5()
    h.update(str(invs).encode())
    return int.from_bytes(h.digest()[:4], "big", signed=True)


def mol_to_pairs(mol):
    """Compute the bond-pair identifiers for all bonds in a molecule.

    For each bond the function calculates two ECFP2-like atom environment
    hashes — one per endpoint — using the immediate neighbourhood (atom + all
    directly bonded neighbours, along with bond type to each).  The two hashes
    are sorted and stored as a tuple so that bond direction does not matter.

    This replaces RDKit's fingerprint generators for this task to avoid the
    overhead of bond fracturing and aromatic re-sanitisation that would be
    needed with the generator API.

    Parameters
    ----------
    mol : RDKit Mol

    Returns
    -------
    list of tuple(int, int)
        One sorted (h1, h2) bond-pair per bond, in bond-index order.
        Returns [] if *mol* is None or unparseable.
    """
    if mol:
        nb = get_neighbors(mol)
        a_invs = get_atom_invariants(mol)
        invs = []
        for b in mol.GetBonds():
            b1 = b.GetBeginAtomIdx()
            b2 = b.GetEndAtomIdx()
            nb1 = copy.copy(nb[b1])
            nb2 = copy.copy(nb[b2])
            nb1.remove(b2)
            nb2.remove(b1)
            bt = int(b.GetBondType())
            n_inv1 = [a_invs[b1] + [bt]] + sorted(
                [a_invs[n] + [int(mol.GetBondBetweenAtoms(b1, n).GetBondType())] for n in nb1])
            n_inv2 = [a_invs[b2] + [bt]] + sorted(
                [a_invs[n] + [int(mol.GetBondBetweenAtoms(b2, n).GetBondType())] for n in nb2])
            h1 = hash_invariants(list(chain.from_iterable(n_inv1)))
            h2 = hash_invariants(list(chain.from_iterable(n_inv2)))
            invs.append(tuple(sorted([h1, h2])))
        return invs
    else:
        print("there was a molecule that didn't parse.")
        return []


def get_neighbors(mol):
    """Return the neighbour atom indices for every atom in *mol*.

    Parameters
    ----------
    mol : RDKit Mol

    Returns
    -------
    list of list of int
        ``result[i]`` is the list of atom indices directly bonded to atom *i*.
    """
    return [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]


def get_atom_invariants(mol):
    """Compute ECFP2-like atom environment descriptors for every atom.

    Each descriptor is a list of five integers:

    1. Atomic number
    2. Heavy-atom degree (number of non-H neighbours)
    3. Total hydrogen count (explicit + implicit)
    4. Formal charge
    5. Ring type: 0 if acyclic, 1 if in a non-aromatic ring, 2 if in an aromatic ring

    These are used directly as the innermost layer of the bond-pair hash
    (see :func:`mol_to_pairs`).

    Parameters
    ----------
    mol : RDKit Mol

    Returns
    -------
    list of list of int
        One 5-element descriptor per atom, in atom-index order.
    """
    invs = []
    for idx, a in enumerate(mol.GetAtoms()):
        inv = [
            a.GetAtomicNum(),
            a.GetDegree(),
            a.GetNumExplicitHs() + a.GetNumImplicitHs(),
            a.GetFormalCharge(),
            int(a.IsInRing())+int(a.GetIsAromatic()),
        ]
        invs.append(inv)
    return invs


def get_profile_for_mols(suppl, profile_name, size=8192, n_jobs=1):
    """Build a LACAN profile from a molecule supplier and save it to disk.

    The profile is a dict ``{"idx": ..., "pairs": ..., "setsize": ...}``
    (see module docstring) pickled to ``data/<profile_name>.pickle`` inside
    the package directory.

    Parameters
    ----------
    suppl        : iterable of RDKit Mol objects (e.g. SmilesMolSupplier)
    profile_name : str — name under which to save the profile
    size         : int — keep only the *size*-1 most common atom environments
                   in ``idx`` (reduces profile file size; default 8192)
    n_jobs       : int — parallel workers; set < 1 for all CPU cores

    Returns
    -------
    dict
        The profile dict (same object that is written to disk).
    """
    if n_jobs == 1:
        all_pairs = [mol_to_pairs(m) for m in suppl if m]
    else:
        if n_jobs < 1:
            n_jobs = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=n_jobs)
        all_pairs = pool.imap_unordered(mol_to_pairs, suppl)
        pool.close()
        pool.join()
    all_pairs = [item for sublist in all_pairs for item in sublist]
    idx = [pair[0] for pair in all_pairs] + [pair[1] for pair in all_pairs]
    idx_occurences = dict(Counter(idx).most_common(size - 1))
    pair_occurences = dict(Counter(all_pairs))
    this_dir, this_filename = os.path.split(__file__)
    DATA_PATH = os.path.join(this_dir, f"data/{profile_name}.pickle")
    with open(DATA_PATH, 'wb') as file:
        pickle.dump({"idx": idx_occurences, "pairs": pair_occurences,
                     "setsize": len(all_pairs)}, file)
    return {"idx": idx_occurences, "pairs": pair_occurences, "setsize": len(all_pairs)}


def assess_per_bond(mol, profile=None):
    """Return the PMI score for every bond in *mol*.

    For each bond the PMI is::

        PMI = (pairs[(e1,e2)] / setsize) / ((idx[e1]/setsize/2) * (idx[e2]/setsize/2))

    Unknown environments (not seen in training) contribute an observed
    occurrence of 0, giving PMI = 0 for that bond.

    Parameters
    ----------
    mol     : RDKit Mol
    profile : LACAN profile dict; loads the default ChEMBL profile if None

    Returns
    -------
    list of float
        One PMI value per bond, in bond-index order.  Values near 0 indicate
        unusual bond environments; values > 1 indicate commonly co-occurring
        environments.
    """
    if profile is None:
        profile = get_default_profile()
    pairs = mol_to_pairs(mol)
    result = []
    for pair in pairs:
        try:
            o1 = profile["idx"][pair[0]] / profile["setsize"] / 2
        except Exception:
            o1 = 0
        try:
            o2 = profile["idx"][pair[1]] / profile["setsize"] / 2
        except Exception:
            o2 = 0
        expected = o1 * o2
        real = profile["pairs"].get(pair, 0) / profile["setsize"]
        result.append(0 if expected == 0 else real / expected)
    return result


def load_profile(profile_name):
    """Load a pickled LACAN profile from the package data directory.

    Parameters
    ----------
    profile_name : str
        Name of the profile (without ``.pickle`` extension).  The file must
        exist at ``<package_dir>/data/<profile_name>.pickle``.

    Returns
    -------
    dict
        Profile dict with keys ``"idx"``, ``"pairs"``, ``"setsize"``.

    Raises
    ------
    FileNotFoundError
        If the pickle file does not exist.
    """
    this_dir, _ = os.path.split(__file__)
    DATA_PATH = os.path.join(this_dir, f"data/{profile_name}.pickle")
    with open(DATA_PATH, 'rb') as file:
        profile = pickle.load(file)
    return profile


_DEFAULT_PROFILE = None


def get_default_profile():
    """Lazy-load the default ChEMBL profile the first time it is requested.

    Caches the result in the module-level ``_DEFAULT_PROFILE`` variable so
    subsequent calls return immediately without re-reading disk.

    Returns
    -------
    dict
        The ChEMBL LACAN profile.
    """
    global _DEFAULT_PROFILE
    if _DEFAULT_PROFILE is None:
        _DEFAULT_PROFILE = load_profile("chembl")
    return _DEFAULT_PROFILE


def score_mol(mol, profile=None, mode="score", t=0.05):
    """Score a molecule using the LACAN profile, respecting protected bonds.

    Protected bonds (marked with ``_lp`` via :mod:`lacan.protect`) are excluded
    from the score calculation.  This means a molecule with required-but-unusual
    motifs (e.g. a covalent warhead) can be scored fairly by protecting those
    bonds first with :func:`~lacan.protect.protect_rejected_bonds`.

    The score is derived from the minimum per-bond PMI over **unprotected**
    bonds only.  Two modes are available:

    ``"score"`` (default)
        Continuous score in [0, 1] computed as::

            score = min_PMI / (1 + min_PMI)

        where ``min_PMI`` is the lowest per-bond PMI over unprotected bonds.
        Saturates toward 1.0 for well-profiled bonds; approaches 0 for unusual
        ones.  The threshold *t* does not affect this value — it only determines
        which bond indices appear in ``info["bad_bonds"]``.

    ``"threshold"``
        Binary: 1 if all unprotected bonds pass (PMI ≥ *t*), else 0.

    Parameters
    ----------
    mol     : RDKit Mol (may have bonds with ``_lp`` protection property)
    profile : LACAN profile dict; loads the default ChEMBL profile if None
    mode    : ``"score"`` or ``"threshold"``
    t       : bond PMI threshold (default 0.05)

    Returns
    -------
    (score, info)
        ``score`` — float in [0, 1] (or 0/1 in threshold mode).
        ``info``  — dict with key ``"bad_bonds"``: list of unprotected bond
        indices whose PMI is below *t*.
    """
    apb = assess_per_bond(mol, profile)
    if not apb:
        apb = [0.0]

    # Identify protected bonds — skip them from scoring
    protected = set()
    for b in mol.GetBonds():
        if b.HasProp("_lp") and b.GetBoolProp("_lp"):
            protected.add(b.GetIdx())

    apb_active = [score for i, score in enumerate(apb) if i not in protected]
    if not apb_active:
        return 1.0, {"bad_bonds": []}  # all bonds protected — trivially passes

    info = {"bad_bonds": [i for i, sc in enumerate(apb)
                          if sc < t and i not in protected]}

    if mode == "threshold":
        score = 0 if min(apb_active) < t else 1
    elif mode == "score":
        score = min(apb_active)/(1+min(apb_active))
    else:
        print("mode not supported yet, sorry.")
        score = 0
        info["bad_bonds"] = []
    return score, info


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lacan CLI")
    parser.add_argument("-i", "--input", type=str,
                        help="input file. should be one smiles per line", required=True)
    parser.add_argument("-m", "--mode", type=str, default="score",
                        help="mode: score or profile", required=False)
    parser.add_argument("-p", "--profile", type=str, default="chembl",
                        help="name of profile to run/generate", required=False)
    parser.add_argument("-s", "--size", type=int, default=8192,
                        help="top N fragments to keep", required=False)
    parser.add_argument("-t", "--threshold", type=float, default=0.05,
                        help="rejection threshold for bond occurrence", required=False)
    parser.add_argument("-c", "--cpus", type=int, default=1,
                        help="number of CPU threads (0 = all)", required=False)
    args = vars(parser.parse_args())

    if args["cpus"] == 1:
        suppl = Chem.SmilesMolSupplier(args["input"], titleLine=False)
    else:
        suppl = Chem.MultithreadedSmilesMolSupplier(
            args["input"], titleLine=False, numWriterThreads=args["cpus"])

    if args["mode"] == "score":
        mols = [m for m in suppl]
        PROFILE = load_profile(args["profile"])
        scores = [score_mol(m, PROFILE, t=args["threshold"], mode="threshold") for m in mols]
        print("overview of failed compounds:")
        print("score\tbad bonds idx\tSMILES")
        for i, s in enumerate(scores):
            print(str(s[0]) + "\t" + str(s[1]["bad_bonds"]) + "\t" + Chem.MolToSmiles(mols[i]))
        print(100 * sum([s[0] for s in scores]) / len(scores), "percent of molecules passed")
    elif args["mode"] == "profile":
        get_profile_for_mols(suppl, args["profile"], size=args["size"], n_jobs=args["cpus"])
