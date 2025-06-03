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
p = Chem.MolFromSmarts("[#0]") # dummy atom
srflag = Chem.SanitizeFlags.SANITIZE_SYMMRINGS

def hash_invariants(invs):
    """
    md5 hash folded to 32bit int to have
    reproducible and portable hashes for environments
    """
    h = hashlib.md5()
    h.update(str(invs).encode())
    return int.from_bytes(h.digest()[:4],signed=True) #32 bit prefix

def mol_to_pairs(mol):
    """
    function that fractures every bond and reports the two ECFP2like
    identifiers at the fracture point.
    New code calculates them explicitly instead of using fpgenetators
    from rdkit. This is to avoid the time sink of bond fracturing and
    aromatic sanitization.
    """
    if mol:
        nb = get_neighbors(mol)
        a_invs = get_atom_invariants(mol)
        btd = get_bond_type_dict(mol)
        invs = []
        for b in mol.GetBonds():
            b1 = b.GetBeginAtomIdx()
            b2 = b.GetEndAtomIdx()
            nb1 = copy.copy(nb[b1])
            nb2 = copy.copy(nb[b2])
            nb1.remove(b2)
            nb2.remove(b1)
            bt = int(b.GetBondType())
            n_inv1 = [a_invs[b1]+[bt]]+ sorted([a_invs[n] + [int(mol.GetBondBetweenAtoms(b1,n).GetBondType())] for n in nb1])
            n_inv2 = [a_invs[b2]+[bt]]+ sorted([a_invs[n] + [int(mol.GetBondBetweenAtoms(b2,n).GetBondType())] for n in nb2])
            h1 = hash_invariants(list(chain.from_iterable(n_inv1)))
            h2 = hash_invariants(list(chain.from_iterable(n_inv2)))
            invs.append(tuple(sorted([h1,h2])))
        return invs
    else:
        print("there was a molecule that didn't parse.")
        return []
    
def get_neighbors(mol):
    """
    return the idx of each atoms directly bounds neighbots
    """
    return [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]

def get_bond_type_dict(mol):
    """
    return dict where btd[i][j] give the bond type of the bond between
    atom with idx i and j
    """
    na = len(mol.GetAtoms())
    btd = {i:{j:0 for j in range(na)} for i in range(na)}
    for b in mol.GetBonds():
        b1 = b.GetBeginAtomIdx()
        b2 = b.GetEndAtomIdx()
    return btd
    
def get_atom_invariants(mol):
    """
    get ECFP like atom identifiers. these are
    - atom number
    - degree
    - h count
    - formal charge
    - smallest ring atom is in. set to 0 if not in ring
    """
    invs = []
    sssr = Chem.GetSSSR(mol)
    min_ring = {}
    for ring in sssr:
        for a in ring:
            if a not in min_ring:
                min_ring[a] = len(ring)
            else:
                if min_ring[a] > len(ring):
                    min_ring[a] = len(ring)
    for idx,a in enumerate(mol.GetAtoms()):
        inv = []
        inv.append(a.GetAtomicNum())
        inv.append(a.GetDegree())
        inv.append(a.GetNumExplicitHs() + a.GetNumImplicitHs())
        inv.append(a.GetFormalCharge())
        try:
            inv.append(min_ring[idx])
        except:
            inv.append(0)
        invs.append(inv)
    return invs


def get_profile_for_mols(suppl,profile_name,size=2048,n_jobs=1):
    if n_jobs==1:
        all_pairs = [mol_to_pairs(m) for m in suppl if m]
    else:
        if n_jobs<1:
            n_jobs = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=n_jobs)
        all_pairs = pool.imap_unordered(mol_to_pairs, suppl) #much better mem usage than pool.map
        pool.close()
        pool.join()
    all_pairs = [item for sublist in all_pairs for item in sublist] #flatten
    idx = [pair[0] for pair in all_pairs] + [pair[1] for pair in all_pairs]
    idx_occurences = dict(Counter(idx).most_common(size-1))
    pair_occurences = dict(Counter(all_pairs))
    this_dir, this_filename = os.path.split(__file__)
    DATA_PATH = os.path.join(this_dir, f"data/{profile_name}.pickle")
    with open(DATA_PATH, 'wb') as file:
        pickle.dump({"idx":idx_occurences,"pairs":pair_occurences,"setsize":len(all_pairs)},file)
    return {"idx":idx_occurences,"pairs":pair_occurences,"setsize":len(all_pairs)}

def assess_per_bond(mol,profile=None):
    if profile==None:
        profile = PROFILE
    pairs = mol_to_pairs(mol)
    assess_per_bond = []
    for pair in pairs:
        try:
            o1 = profile["idx"][pair[0]]/profile["setsize"]/2
        except Exception as e:
            o1 = 0
        try: 
            o2 = profile["idx"][pair[1]]/profile["setsize"]/2
        except Exception as e:
            o2 = 0
        expected_occurence = o1*o2
        if pair in profile["pairs"]:
            real_occurence = profile["pairs"][pair]/profile["setsize"]
        else:
            real_occurence = 0
        if expected_occurence == 0:
            assess_per_bond.append(0)
        else:
            assess_per_bond.append(real_occurence/expected_occurence)
    return assess_per_bond
    
def load_profile(profile_name):
    this_dir, this_filename = os.path.split(__file__)
    DATA_PATH = os.path.join(this_dir, f"data/{profile_name}.pickle")
    with open(DATA_PATH, 'rb') as file:
        profile = pickle.load(file)
    return profile
        
def score_mol(mol,profile=None,mode="score",t=0.05):
    apb = assess_per_bond(mol,profile)
    info = {}
    if len(apb) == 0:
        apb=[0]
    if mode == "threshold":
        info["bad_bonds"] = [i for i,b in enumerate(apb) if b < t]
        if min(apb) < t:
            score = 0
        else:
            score = 1
    elif mode == "score":
        info["bad_bonds"] = [i for i,b in enumerate(apb) if b < t]
        #score set so when threshold is reached, score is 0.5. 
        score = min(0.5*(min(apb)/t)**0.5,1.0)
    else:
        print("mode not supported yet, sorry.")
    return score, info

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lacan CLI")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="input file. should be one smiles per line and nothing else",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        default="score",
        help="mode to run. you can choose between score, profile",
        required=False,
    )
    parser.add_argument(
        "-p",
        "--profile",
        type=str,
        default="chembl",
        help="name of profile to run/generate (depending on mode)",
        required=False,
    )
    parser.add_argument(
        "-s",
        "--size",
        type=int,
        default=2048,
        help="top N fragments to keep",
        required=False,
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.05,
        help="rejection threshold for lower than expected occurence for bonds",
        required=False,
    )
    
    parser.add_argument(
        "-c",
        "--cpus",
        type=int,
        default=1,
        help="number of cpu threads to use (set to 0 for all)",
        required=False,
    )

    
    args = vars(parser.parse_args())
    if args["cpus"] == 1:
        suppl = Chem.SmilesMolSupplier(args["input"],titleLine=False,)
    else:
        suppl = Chem.MultithreadedSmilesMolSupplier(args["input"],titleLine=False,numWriterThreads=args["cpus"])
    if args["mode"] == "score":
        mols = [m for m in suppl]
        PROFILE = load_profile(args["profile"])
        scores = [score_mol(m,PROFILE,t=args["threshold"],mode="threshold") for m in mols]
        print("overview of failed compounds:")
        print("score\tbad bonds idx\tSMILES")
        for i,s in enumerate(scores):
            if s[0] == 0:
                print(str(s[0])+"\t"+str(s[1]["bad_bonds"])+"\t"+Chem.MolToSmiles(mols[i]))
            else:
                print(str(s[0])+"\t"+str(s[1]["bad_bonds"])+"\t"+Chem.MolToSmiles(mols[i]))
        print(100*sum([s[0] for s in scores])/len(scores),"percent of molecules passed")
    elif args["mode"] == "profile":
        get_profile_for_mols(suppl,args["profile"],size=args["size"],n_jobs=args["cpus"])
