"""
decompose.py — Molecule fragmentation into rings, linkers, and substituents.

This module breaks drug-like molecules into three canonical fragment types and
provides utilities for building and saving fragment corpora.

Fragment taxonomy
-----------------
After exhaustive recursive cutting at all non-ring exo-bonds, each resulting
fragment is classified by the number and position of its dummy atoms (``*``):

**Ring**
    Fragment contains at least one ring atom adjacent to a dummy (matched by
    ``ringdummy``).  Examples: ``[*]c1ccccc1``, ``[*]c1ccncc1[*]``.

**Linker**
    Fragment has no ring atoms adjacent to any dummy, but has **two or more**
    non-ring dummy neighbours (matched by ``nonringdummy`` twice or more).
    Examples: ``[*]CC[*]``, ``[*]C(=O)[*]``.

**Substituent**
    Fragment has exactly **one** non-ring dummy neighbour — a terminal group.
    Examples: ``[*]CH3``, ``[*]CF3``, ``[*]OH``.

Single-atom molecules, atoms that failed sanitisation, or molecules that do not
decompose further (no non-ring exo-bonds) are silently ignored.

Reaction SMARTS
---------------
``decompose1``  cuts single-bond ring-exo bonds: ``[!#0]!@;-[R] >> [!R]-[*].[R]-[*]``
``decompose2``  cuts double-bond ring-exo bonds: ``[!#0]!@;=[R] >> [!R]=[*].[R]=[*]``

Note the ``!@`` (not-in-ring) flag ensures only exo-bonds are cut, not ring
bonds themselves.  The ``!#0`` guard prevents re-cutting already-cut dummy
atoms.

Corpus format
-------------
Each corpus entry is a list::

    [smiles, count, degree, ftype, bonds]

where ``smiles`` is the canonical fragment SMILES (with ``*`` dummies),
``count`` is its occurrence frequency, ``degree`` is the number of attachment
points, ``ftype`` is ``"Ring"``, ``"Linker"``, or ``"Sub"``, and ``bonds`` is
a string like ``"-"``, ``"--"``, or ``"=-"`` encoding the bond types at each
dummy in order of dummy atom index.
"""

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions, inchi
import random, argparse, csv
from collections import Counter
import multiprocessing

RDLogger.DisableLog('rdApp.*')

decompose1 = rdChemReactions.ReactionFromSmarts(
    "[!#0:0]!@;-[R:1]>>[!R:0]-[#0].[R:1]-[#0]"
)
"""Cut a single-bond ring-exo bond.

Yields two fragments: the non-ring part with a single-bond dummy, and the ring
part with a single-bond dummy.  The ``!#0`` guard prevents re-cutting dummies
from a previous iteration.  Modified from the original to also decompose
biphenyls and other directly-fused ring systems.
"""

decompose2 = rdChemReactions.ReactionFromSmarts(
    "[!#0:0]!@;=[R:1]>>[!R:0]=[#0].[R:1]=[#0]"
)
"""Cut a double-bond ring-exo bond (e.g. a carbonyl attached to a ring)."""

ringdummy = Chem.MolFromSmarts("[#0]~[R]")
"""Matches a dummy atom adjacent to a ring atom — identifies Ring fragments."""

nonringdummy = Chem.MolFromSmarts("[#0]~[!R]")
"""Matches a dummy atom adjacent to a non-ring atom — identifies Linker/Sub fragments."""

singledummy = Chem.MolFromSmarts("[#0;$([*]-[*])]")
"""Matches a single-bond dummy attachment point."""

doubledummy = Chem.MolFromSmarts("[#0;$([*]=[*])]")
"""Matches a double-bond dummy attachment point."""


def get_bonds_string(smi):
    """Encode the bond types at all dummy attachment points as a string.

    The string contains one character per dummy atom, sorted by dummy atom
    index: ``"-"`` for single bonds and ``"="`` for double bonds.  This is
    stored in the corpus so that fragment matching can filter by compatible
    bond types.

    Examples: ``"[*]CC"`` → ``"-"``, ``"[*]CC[*]"`` → ``"--"``,
    ``"[*]C(=O)[*]"`` → ``"-="``.

    Parameters
    ----------
    smi : str — fragment SMILES containing at least one ``*`` dummy atom

    Returns
    -------
    str
        Bond-type string, one character per dummy, sorted by dummy atom index.
    """
    mol = Chem.MolFromSmiles(smi)
    bonds = []
    for db in mol.GetSubstructMatches(singledummy):
        bonds.append(("-", db[0]))
    for db in mol.GetSubstructMatches(doubledummy):
        bonds.append(("=", db[0]))
    bonds = sorted(bonds, key=lambda x: x[1])
    return "".join([b[0] for b in bonds])


def decompose_molecule(smi, asSmi=True):
    """Recursively decompose a molecule into rings, linkers, and substituents.

    The algorithm repeatedly applies ``decompose1`` and ``decompose2`` to cut
    all non-ring exo-bonds, one at a time, until no further cuts are possible.
    The resulting terminal fragments are then classified into rings, linkers,
    and substituents based on their dummy-atom topology (see module docstring).

    Stereochemistry is removed before decomposition so that stereoisomers map
    to the same fragment SMARTS.

    Single-fragment decompositions (molecules with no exo-bonds, e.g. bare
    rings) are discarded — they yield no useful corpus entries.

    Parameters
    ----------
    smi   : str or RDKit Mol — the molecule to decompose
    asSmi : bool — if True (default), treat *smi* as a SMILES string and parse
            it; if False, treat it directly as an RDKit Mol object

    Returns
    -------
    (rings, linkers, subs) : three lists of canonical SMILES strings
    """
    rings = []
    linkers = []
    subs = []
    mol = Chem.MolFromSmiles(smi) if asSmi else smi
    if mol:
        Chem.RemoveStereochemistry(mol)
        frags = [mol]
        final_frags = []
        while len(frags) > 0:
            prods = decompose1.RunReactants((frags[0],))
            prods += decompose2.RunReactants((frags[0],))
            if len(prods) > 0:
                try:
                    frags.pop(0)
                    Chem.SanitizeMol(prods[0][0])
                    Chem.SanitizeMol(prods[0][1])
                    frags += [prods[0][0], prods[0][1]]
                except Exception as e:
                    frags = []
                    print(f"decompose_molecule: sanitization failed ({e}), skipping molecule")
            else:
                final_frags.append(frags[0])
                frags.pop(0)
        if len(final_frags) > 1:
            for frag in final_frags:
                rec = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(frag)))
                if frag.HasSubstructMatch(ringdummy):
                    rings.append(rec)
                else:
                    sm = len(frag.GetSubstructMatches(nonringdummy))
                    if sm > 1:
                        linkers.append(rec)
                    elif sm == 1:
                        subs.append(rec)
    return rings, linkers, subs


def get_all_decompositions(smis, n_jobs=-1):
    """Decompose a large list of SMILES strings in parallel.

    Distributes :func:`decompose_molecule` across a multiprocessing pool using
    ``imap_unordered`` for memory-efficient streaming (important for millions
    of ChEMBL molecules).

    Parameters
    ----------
    smis   : iterable of SMILES strings
    n_jobs : int — number of worker processes; -1 (default) uses all CPU cores

    Returns
    -------
    (all_rings, all_linkers, all_subs) : three flat lists of SMILES strings
        (unsorted, may contain duplicates — use :func:`save_frags` to count
        and filter)
    """
    all_rings = []
    all_linkers = []
    all_subs = []
    if n_jobs == -1:
        n_jobs = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=n_jobs)
    rlss = pool.imap_unordered(decompose_molecule, smis, chunksize=10000)
    pool.close()
    pool.join()
    for r, l, s in rlss:
        all_rings += r
        all_linkers += l
        all_subs += s
    return all_rings, all_linkers, all_subs



def save_frags(r, l, s, output_path, minN=10):
    """Count fragment occurrences and save a filtered corpus CSV.

    Only fragments appearing more than ``minN`` times are written.  The output
    CSV has columns ``smiles, occurrence, degree, type, bonds`` and is sorted
    by occurrence descending.  The resulting file can be loaded with
    ``replace.load_corpus(path="my_corpus.csv")`` as a drop-in replacement
    for the built-in ChEMBL corpus.

    Parameters
    ----------
    r           : list of ring SMILES (from :func:`get_all_decompositions`)
    l           : list of linker SMILES
    s           : list of substituent SMILES
    output_path : str — path to write the CSV
    minN        : int — minimum occurrence count to include (default 10)
    """
    cr = dict(Counter(r))
    cl = dict(Counter(l))
    cs = dict(Counter(s))
    entries = []
    for k, cnt in cr.items():
        if cnt > minN:
            entries.append([k, cnt, k.count("*"), "Ring", get_bonds_string(k)])
    for k, cnt in cl.items():
        if cnt > minN:
            entries.append([k, cnt, k.count("*"), "Linker", get_bonds_string(k)])
    for k, cnt in cs.items():
        if cnt > minN:
            entries.append([k, cnt, k.count("*"), "Sub", get_bonds_string(k)])
    entries = sorted(entries, key=lambda x: -x[1])
    with open(output_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile, delimiter=",", quotechar='"',
                            quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["smiles", "occurrence", "degree", "type", "bonds"])
        for entry in entries:
            writer.writerow(entry)
    print(f"Wrote {len(entries)} fragments to {output_path}")

def get_corpus(mols):
    """Build a fragment corpus from a list of RDKit Mol objects.

    Decomposes each molecule and counts all fragment occurrences.  Unlike
    :func:`save_frags`, this function returns the corpus as a list in memory
    (no minimum count filter, no file I/O) so it can be passed directly to
    :func:`~lacan.gen.get_corpus_from_mols` for corpus biasing.

    Parameters
    ----------
    mols : iterable of RDKit Mol objects (None entries are silently skipped)

    Returns
    -------
    list of [smiles, count, degree, ftype, bonds]
        Sorted by count descending.  Contains all fragment types with count ≥ 1.
    """
    r, l, s = [], [], []
    for mol in mols:
        if mol:
            r0, l0, s0 = decompose_molecule(mol, asSmi=False)
            r += r0
            l += l0
            s += s0
    cr = dict(Counter(r))
    cl = dict(Counter(l))
    cs = dict(Counter(s))
    entries = []
    for k in cr:
        entries.append([k, cr[k], k.count("*"), "Ring", get_bonds_string(k)])
    for k in cl:
        entries.append([k, cl[k], k.count("*"), "Linker", get_bonds_string(k)])
    for k in cs:
        entries.append([k, cs[k], k.count("*"), "Sub", get_bonds_string(k)])
    entries = sorted(entries, key=lambda x: -x[1])
    return entries



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Build a LACAN fragment corpus from a SMILES/CSV file.\n\n"
            "Examples:\n"
            "  # ChEMBL (space-separated, SMILES col 0, has header):\n"
            "  python -m lacan.decompose -i chembl.smi -o chembl_rls.csv --delimiter ' ' --skip-header\n\n"
            "  # COCONUT (comma-separated, SMILES col 0):\n"
            "  python -m lacan.decompose -i coconut.csv -o coconut_rls.csv\n\n"
            "Load the result with:\n"
            "  corpus = replace.load_corpus(path='coconut_rls.csv')"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-i", "--input",  required=True, help="Input SMILES or CSV file.")
    parser.add_argument("-o", "--output", required=True, help="Output CSV path.")
    parser.add_argument("--delimiter", default=",", help="Column delimiter (default: comma).")
    parser.add_argument("--smiles-col", type=int, default=0, help="Zero-based SMILES column index (default: 0).")
    parser.add_argument("--skip-header", action="store_true", help="Skip the first (header) line.")
    parser.add_argument("--min-count", type=int, default=10, help="Minimum fragment occurrence (default: 10).")
    parser.add_argument("--n-jobs", type=int, default=-1, help="Parallel workers; -1 = all CPUs.")
    args = parser.parse_args()

    smis = []
    with open(args.input, newline="") as f:
        reader = csv.reader(f, delimiter=args.delimiter)
        for i, row in enumerate(reader):
            if i == 0 and args.skip_header:
                continue
            if row and len(row) > args.smiles_col:
                smi = row[args.smiles_col].strip()
                if smi:
                    smis.append(smi)

    print(f"Read {len(smis):,} SMILES from {args.input}")
    r, l, s = get_all_decompositions(smis, n_jobs=args.n_jobs)
    print(f"Decomposed: {len(r):,} ring  {len(l):,} linker  {len(s):,} substituent")
    save_frags(r, l, s, output_path=args.output, minN=args.min_count)
