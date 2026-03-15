"""
breed.py — Molecular crossover via fragment-based recombination.

This module implements the crossover operation used by the adaptive GA
in gen.py.  Two molecules are cut at their non-ring bonds and their fragments
are recombined: substituents from one molecule are attached to the core of the
other, producing offspring that inherit structural features from both parents.

Algorithm
---------
1. **Fragmentation** (:func:`fragment_molecule`) — cut the molecule at *n*
   non-ring bonds in all combinations.  Keep only cuts that yield exactly *n*
   single-dummy substituents and one *n*-dummy core.

2. **Recombination** (:func:`crossover_fragments`) — randomly pair a core from
   one molecule with substituents from the other, fill all attachment points,
   and accept offspring that:

   * pass the LACAN score threshold (``score_mol > 0.5``;
     uses the current ``min_PMI / (1 + min_PMI)`` formula),
   * have heavy-atom count in ``[hacmin, hacmax]``,
   * have a "balance ratio" (fraction of atoms from parent 2) in
     ``[min_ratio, max_ratio]`` — this ensures the offspring genuinely
     interpolates between the two parents rather than being almost identical
     to one of them.

3. **Deduplication** by InChIKey — no duplicate offspring are returned.

Fallback
--------
If three-cut fragmentation yields no substituents for a molecule (e.g. because
it has fewer than three non-ring bonds), :func:`breed` falls back to two-cut
fragmentation.  This behaviour is silent by default; set ``debug=True`` to
print a message when the fallback fires.
"""

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions, inchi
from lacan import lacan
import argparse
import random
from itertools import combinations
import multiprocessing

RDLogger.DisableLog('rdApp.*')

combine_frags = rdChemReactions.ReactionFromSmarts(
    '[*:0][#0:1].[*:2][#0:3]>>[*:0][*:2].[*:1][*:3]'
)
"""Join two dummy attachment points (any bond type)."""

def fragment_molecule(mol,n=3):
    """Cut a molecule at *n* non-ring bonds and return substituents and cores.

    All combinations of *n* non-ring bond indices are tried.
    A cut is retained only if it produces exactly *n* single-dummy fragments
    (substituents) and one *n*-dummy fragment (the core).

    Parameters
    ----------
    mol : RDKit Mol
    n   : int — number of cuts to make (default 3)

    Returns
    -------
    (substituents, cores) : (list of str, list of str)
        Canonical SMILES strings with ``*`` dummies.
        Either list may be empty if no valid cuts exist at the requested depth.
    """
    substituents = set([])
    cores = set([])
    nonringbonds = [b.GetIdx() for b in mol.GetBonds() if b.IsInRing()==False]
    neededlengths = [1]*n+[n]
    for bonds in combinations(nonringbonds,n):
        newmol = Chem.FragmentOnBonds(mol,bonds)
        newsmi = Chem.MolToSmiles(newmol)
        frags = newsmi.split(".")
        lengths = sorted([x.count("*") for x in frags])
        if lengths == neededlengths:
            for frag in frags:
                if frag.count("*")==1:
                    substituents.add(frag)
                else:
                    cores.add(frag)
    return list(substituents),list(cores)
    
def crossover_fragments(s1, s2, c1, c2, profile, nmols=10, randomseed=123, max_steps=1000,
                        hacmin=0, hacmax=30, min_ratio=0.25, max_ratio=0.75):
    """Recombine substituents and cores from two molecules to produce offspring.

    At each step a coin-flip chooses which parent provides the core and which
    provides the substituents.  All attachment points of the core are filled
    with randomly sampled substituents, and the offspring is accepted if it
    passes the LACAN filter and size/ratio constraints.

    The ``ratio`` is ``atoms_from_parent2 / total_atoms``; restricting it to
    ``[min_ratio, max_ratio]`` ensures the offspring genuinely blends both
    parents rather than being almost identical to one.

    Parameters
    ----------
    s1, s2      : list of str — substituent SMILES from parents 1 and 2
    c1, c2      : list of str — core SMILES from parents 1 and 2
    profile     : LACAN profile dict
    nmols       : int — target number of offspring to produce (default 10)
    randomseed  : int — random seed for reproducibility (default 123)
    max_steps   : int — maximum sampling attempts before giving up (default 1000)
    hacmin      : int — minimum heavy-atom count (default 0)
    hacmax      : int — maximum heavy-atom count (default 30)
    min_ratio   : float — minimum parent-2 atom fraction (default 0.25)
    max_ratio   : float — maximum parent-2 atom fraction (default 0.75)

    Returns
    -------
    list of RDKit Mol — deduplicated offspring (may be fewer than *nmols* if
    max_steps is reached first)
    """
    mols = []
    inchis = []
    steps = 0
    random.seed(randomseed)
    while len(mols)<nmols and steps<max_steps:
        atoms1 = 0
        atoms2 = 0
        coinflip = 0
        if random.random() > 0.5:
            s=s1
            c=random.sample(c2,1)[0]
            coinflip = 0
        else:
            s=s2
            c=random.sample(c1,1)[0]
            coinflip = 1
        cmol = Chem.MolFromSmiles(c)
        if coinflip == 1:
            atoms1 = len(cmol.GetAtoms())
        else:
            atoms2 = len(cmol.GetAtoms())
        try:
            for r in range(c.count("*")):
                smol = Chem.MolFromSmiles(random.sample(s,1)[0])
                cmol = combine_frags.RunReactants((cmol,smol))[0][0]
                Chem.SanitizeMol(cmol)
                if coinflip == 1:
                    atoms2 += len(smol.GetAtoms())
                else:
                    atoms1 += len(smol.GetAtoms())
            ratio = atoms2/(atoms1+atoms2)
            if lacan.score_mol(cmol,profile)[0]>0.5:
                if hacmin < len(cmol.GetAtoms()) < hacmax:
                    if min_ratio < ratio < max_ratio:
                        ik = inchi.MolToInchiKey(cmol)
                        if  ik not in inchis:
                            mols.append(cmol)
                            inchis.append(ik)
        except Exception as e:
            pass
        steps += 1
    return mols
    
def breed(m1, m2, profile, nmols=10, cuts=3, hacrange=(0.8, 1.2), interprange=(0.3, 0.7), debug=False):
    """Cross two molecules and return a list of offspring molecules.

    Fragments both parents at ``cuts`` non-ring bonds, then calls
    :func:`crossover_fragments` to combine their pieces.  The heavy-atom count
    range for offspring is derived from the parents' sizes via ``hacrange``:
    ``[hacrange[0] * min(n1,n2), hacrange[1] * max(n1,n2)]``.

    If three-cut fragmentation fails for either parent (too few non-ring bonds),
    the function silently falls back to two-cut fragmentation.  Set
    ``debug=True`` to see a message when this happens.

    Parameters
    ----------
    m1, m2      : RDKit Mol — parent molecules
    profile     : LACAN profile dict
    nmols       : int — number of offspring to request (default 10)
    cuts        : int — number of cuts to make in each parent (default 3)
    hacrange    : (float, float) — (min_frac, max_frac) multiplied by parent
                  sizes to set the offspring heavy-atom count window
    interprange : (float, float) — (min_ratio, max_ratio) fraction of atoms
                  from parent 2 (controls interpolation balance)
    debug       : bool — if True, print a message when falling back to 2-cut
                  fragmentation (default False)

    Returns
    -------
    list of RDKit Mol
    """
    s1,c1 = fragment_molecule(m1,cuts)
    s2,c2 = fragment_molecule(m2,cuts)
    n1 = len(m1.GetAtoms())
    n2 = len(m2.GetAtoms())
    if len(s1) == 0 or len(s2) == 0:
        if debug:
            print("didnt find a way to cut >2 times, reverting to twice")
        s1,c1 = fragment_molecule(m1,2)
        s2,c2 = fragment_molecule(m2,2)
    mols = crossover_fragments(s1,s2,c1,c2,profile,nmols,hacmin=int(hacrange[0]*min(n1,n2)),hacmax=int(hacrange[1]*max(n1,n2)),min_ratio=interprange[0],max_ratio=interprange[1])
    return mols
    
def cross_breed_mols(mols, p, score_threshold, nmols=1, n_jobs=-1, debug=False):
    """Apply :func:`breed` across an entire population in parallel.

    Each molecule in *mols* is paired with a randomly chosen partner from the
    same list, and :func:`breed` is called for each pair via a
    ``multiprocessing.Pool``.  Workers always run with ``debug=False``
    regardless of the caller's *debug* flag, because worker processes cannot
    print to a Jupyter notebook.

    Parameters
    ----------
    mols           : list of RDKit Mol — population to breed
    p              : LACAN profile dict
    score_threshold: passed through to :func:`breed` (currently unused there,
                     kept for API consistency with other operations)
    nmols          : int — offspring requested per pair (default 1)
    n_jobs         : int — parallel workers; -1 uses all CPU cores
    debug          : bool — not forwarded to workers (see note above)

    Returns
    -------
    list of RDKit Mol — flat list of all offspring from all pairs
    """
    inputs = []
    if n_jobs==-1:
        n_jobs = multiprocessing.cpu_count()
    for m in mols:
        # Pass debug=False always in multiprocessing — workers can't print to notebook
        inputs.append((m, random.sample(mols,1)[0], p, nmols, 3, (0.8,1.2), (0.3,0.7), False))
    pool = multiprocessing.Pool(processes=n_jobs)
    muts = pool.starmap(breed, inputs)
    pool.close()
    pool.join()
    result = []
    for mut in muts:
        result += mut
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lacan Breed CLI")
    parser.add_argument("-i1", "--input1", type=str, help="input file. should be a smiles", required=True)
    parser.add_argument("-i2", "--input2", type=str, help="input file 2. should be a smiles", required=True)
    parser.add_argument("-p", "--profile", type=str, default="chembl", required=False)
    parser.add_argument("-n", "--nmols", type=int, default=10, required=False)
    parser.add_argument("--debug", action="store_true", default=False)
    args = vars(parser.parse_args())
    m1 = Chem.MolFromSmiles(args["input1"])
    m2 = Chem.MolFromSmiles(args["input2"])
    p = lacan.load_profile(args["profile"])
    for m in breed(m1, m2, p, args["nmols"], debug=args["debug"]):
        print(Chem.MolToSmiles(m))
