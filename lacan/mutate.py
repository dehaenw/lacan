"""
mutate.py — Atom-level mutation operations for LACAN.

This module provides fine-grained single-step mutations: adding, removing, or
changing individual atoms, bonds, or ring sizes.  These are the operations can
be considered as the exploitation operation in the adaptive GA.

Each reaction in ``mutate_smarts`` targets a specific chemotype change.
The dictionary is compiled into ``mutate_ops`` (RDKit Reaction objects) at
import time.

Protection
----------
Before running any reaction, :func:`apply_mutations` calls
:func:`~lacan.protect.reaction_touches_protected` to check whether the
reaction's SMARTS template would match a protected atom.  If so, the entire
reaction is skipped for that molecule.

Multiprocessing
---------------
Molecule atom/bond properties (including protection marks stored under the
``_lp`` property) would normally be lost when RDKit Mol objects are pickled
for ``multiprocessing.Pool``.  :func:`apply_mutations_mols` works around this
by serialising mols to SDF block strings (which preserve properties) before
dispatching to worker processes.
"""

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions, inchi
from lacan import lacan
from lacan.protect import reaction_touches_protected, score_mol_ignoring_protected_bonds
import random, argparse
import multiprocessing

RDLogger.DisableLog('rdApp.*')

mutate_smarts = {
    "addC":        "[H1,H2,H3:0]>>[*:0][CH3]",
    "addO":        "[H1,H2,H3:0]>>[*:0][OH]",
    "addCO":       "[CH2:0]>>[*:0]=[O]",
    "addN":        "[H1,H2,H3:0]>>[*:0][NH2]",
    "addF":        "[H1,H2,H3:0]>>[*:0][F]",
    "addCl":       "[H1,H2,H3:0]>>[*:0][Cl]",
    "addBr":       "[H1,H2,H3:0]>>[*:0][Br]",
    "contractAroO":  "[ar6:0]:[cH,nH0;r6:1]:[cH,nH0;r6:2]:[ar6:3]>>[*:0]:[o]:[*:3].[*:1][*:2]",
    "contractAroS":  "[ar6:0]:[cH,nH0;r6:1]:[cH,nH0;r6:2]:[ar6:3]>>[*:0]:[s]:[*:3].[*:1][*:2]",
    "contractAroNH": "[ar6:0]:[cH,nH0;r6:1]:[cH,nH0;r6:2]:[ar6:3]>>[*:0]:[nH]:[*:3].[*:1][*:2]",
    "expandAroCC": "[ar5:0]:[nH,o,s;r5:1]:[ar5:2]>>[*:0]:[cH]:[cH]:[*:2].[*:1]",
    "expandAroCN": "[ar5:0]:[nH,o,s;r5:1]:[ar5:2]>>[*:0]:[cH]:[n]:[*:2].[*:1]",
    "insertC":     "[*:0]-[*:1]>>[*:0]-[CH2]-[*:1]",
    "insertN":     "[*:0]-[*:1]>>[*:0]-[NH]-[*:1]",
    "insertO":     "[*:0]-[*:1]>>[*:0]-[O]-[*:1]",
    "insertS":     "[*:0]-[*:1]>>[*:0]-[S]-[*:1]",
    "replaceC":    "[!C;A;d4,d3,d2,d1;v4,v3,v2,v1:0]>>[C:0]",
    "replaceN":    "[!N;!$([CH0]);d3,d2,d1;A:0]>>[N:0]",
    "replaceO":    "[!O;!$([CH0]);$([*](-[*])(-[*]));d2,d1;A:0]>>[O:0]",
    "replaceS":    "[!S;!$([CH0]);$([*](-[*])(-[*]));d2,d1;A:0]>>[S:0]",
    "aroCtoN":     "[cH:0]>>[nH0:0]",
    "aroNtoC":     "[nH0X2:0]>>[cH:0]",
    "openring":    "[R:0]@!:[R:1]>>([*:0].[*:1])",
    "close3ring":  "[!R;H1,H2,H3:0][*:1][!R;H1,H2,H3:2]>>[*:0]1~[*:1]~[*:2]1",
    "close4ring":  "[!R;H1,H2,H3:0][*:1][*:2][!R;H1,H2,H3:3]>>[*:0]1~[*:1]~[*:2]~[*:3]1",
    "close5ring":  "[!R;H1,H2,H3:0][*:1][*:2][!a:3][!R;H1,H2,H3:4]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]1",
    "close6ring1": "[!R;H1,H2,H3:0][*:1][*:2][!a:3][!a:4][!R;H1,H2,H3:5]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]~[*:5]1",
    "close6ring2": "[!R;H1,H2,H3:0][!a:1][*:2][*:3][!a:4][!R;H1,H2,H3:5]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]~[*:5]1",
    "arofuse5":    "[aH1:0]:[a:1]-[*:2][*:3]-[!R;H1,H2,H3:4]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]1",
    "arofuse6na1": "[aH1:0]:[a:1]-[*:2][*:3][!a:4]-[!R;H1,H2,H3:5]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]~[*:5]1",
    "arofuse6na2": "[aH1:0]:[a:1]-[!a:2][*:3][*:4]-[!R;H1,H2,H3:5]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]~[*:5]1",
    "bond2to3":    "[!R;H1,H2:0]=[!R;H1,H2:1]>>[*:0]#[*:1]",
    "bond1to2":    "[H1,H2,H3:0]-[H1,H2,H3:1]>>[*:0]=[*:1]",
    "bond2to1":    "[A:0]=[A:1]>>[*:0]-[*:1]",
    "bond3to2":    "[*:0]#[*:1]>>[*:0]=[*:1]",
    "deleteD1":    "[!$([nX3]):0][d1:1]>>[*:0].[*:1]",
    "deleteD2":    "[*:0][d2;A:1][*:2]>>[*:0][*:2].[*:1]",
}
"""Dictionary of named reaction SMARTS for single-step atom-level mutations.

Covers: atom addition (C/O/N/F/Cl/Br), ring contraction/expansion (5↔6,
heteroatom swaps), chain insertion (C/N/O/S), element replacement, aromatic
C↔N swaps, ring opening/closure (3–6-membered), ring fusion, bond order
changes (single/double/triple), and degree-1/2 atom deletion.
"""

mutate_ops = {name: rdChemReactions.ReactionFromSmarts(mutate_smarts[name])
              for name in mutate_smarts}
"""Pre-compiled RDKit Reaction objects for all mutations in ``mutate_smarts``."""


def _apply_mutations_worker(args):
    """Multiprocessing worker: deserialise mol, mutate, re-serialise products.

    Mols are passed as SDF block strings (not pickled Mol objects) so that
    bond properties — in particular the ``_lp`` bond protection mark — survive
    the round-trip through ``multiprocessing.Pool``.  Note that
    ``protect_smarts`` is re-derived on each call, so it does not need
    to survive serialisation.

    Parameters
    ----------
    args : (sdf_block: str, p: dict, score_threshold: float)

    Returns a list of SDF block strings, one per accepted mutation product.
    """
    sdf_block, p, score_threshold = args
    mol = Chem.MolFromMolBlock(sdf_block, removeHs=False)
    results = apply_mutations(mol, p, score_threshold)
    return [Chem.MolToMolBlock(m) for m in results]


def apply_mutations_mols(mols, p, score_threshold, n_jobs=-1):
    """Apply all mutations to a list of molecules, optionally in parallel.

    This is the batch version of :func:`apply_mutations`.  When ``n_jobs != 1``
    it uses ``multiprocessing.Pool`` and serialises mols via SDF blocks to
    preserve atom/bond properties.

    Parameters
    ----------
    mols            : list of RDKit Mol objects
    p               : LACAN profile dict
    score_threshold : minimum LACAN score for output molecules
    n_jobs          : number of parallel workers; -1 uses all CPU cores

    Returns a flat list of RDKit Mol objects (all products from all input mols).
    """
    if n_jobs == 1:
        result = []
        for m in mols:
            result += apply_mutations(m, p, score_threshold)
        return result
    else:
        if n_jobs == -1:
            n_jobs = multiprocessing.cpu_count()
        args = [(Chem.MolToMolBlock(m), p, score_threshold) for m in mols]
        pool = multiprocessing.Pool(processes=n_jobs)
        muts = pool.map(_apply_mutations_worker, args)
        pool.close()
        pool.join()
        result = []
        for blocks in muts:
            result += [Chem.MolFromMolBlock(b, removeHs=False) for b in blocks]
        return result


def apply_mutations(mol, p, score_threshold, mode="all", protect_smarts=None):
    """Apply single-step atom-level mutations to a molecule and return passing variants.

    For each reaction in ``mutate_ops`` (or a randomly chosen one if
    ``mode="random"``):

    1. Skip the reaction if it would touch an atom matching ``protect_smarts``
       (checked via :func:`~lacan.protect.reaction_touches_protected`).
    2. Run the reaction on the molecule; sanitize each product.
    3. Score each product with :func:`~lacan.protect.score_mol_ignoring_protected_bonds`.
    4. Keep products whose score exceeds ``score_threshold``.

    Deduplication is by InChIKey.

    Parameters
    ----------
    mol             : RDKit Mol (may have protected bonds)
    p               : LACAN profile dict
    score_threshold : minimum LACAN score; set 0.0 to accept all LACAN-passing
                      products, 0.8 for stricter drug-likeness filtering
    mode            : ``"all"`` (try every reaction) or ``"random"`` (one random
                      reaction)
    protect_smarts  : SMARTS string; reactions touching any matching atom are
                      skipped entirely.  ``None`` = no exclusion (default).

    Returns a deduplicated list of RDKit Mol objects.

    Note
    ----
    Products of reactions do **not** inherit the parent's bond protection marks.
    This is intentional for normal use (the score filter selects good products),
    but :func:`~lacan.protect.mol_cleaner` bypasses this function in favour of
    :func:`~lacan.protect._raw_mutations` to avoid discarding partial fixes.
    """
    mutate_prods = []
    if mode == "all":
        ops = mutate_ops
    elif mode == "random":
        random_key = random.sample(list(mutate_ops.keys()), 1)[0]
        ops = {random_key: mutate_ops[random_key]}
    else:
        print("this mode doesnt exist")
        ops = {}

    for op in ops:
        rxn = ops[op]
        if reaction_touches_protected(mol, rxn, protect_smarts):
            continue
        prods = rxn.RunReactants((mol,))
        if len(prods) > 0:
            for prod in prods:
                try:
                    Chem.SanitizeMol(prod[0])
                    if prod[0]:
                        mutate_prods.append(prod[0])
                except Exception:
                    pass

    filtered_prods = []
    iks = []
    for m in mutate_prods:
        score, info = score_mol_ignoring_protected_bonds(m, p)
        if score > score_threshold:
            ik = inchi.MolToInchiKey(m)
            if ik not in iks:
                filtered_prods.append(m)
                iks.append(ik)
    return filtered_prods


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lacan Mutate CLI")
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-p", "--profile", type=str, default="chembl")
    parser.add_argument("-m", "--mode", type=str, default="all")
    parser.add_argument("-t", "--threshold", type=float, default=0.8)
    args = vars(parser.parse_args())
    p = lacan.load_profile(args["profile"])
    mol = Chem.MolFromSmiles(args["input"])
    for m in apply_mutations(mol, p, args["threshold"], args["mode"]):
        print(Chem.MolToSmiles(m))
