from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from collections import Counter
import sys, os
import pickle
import argparse
import math

RDLogger.DisableLog('rdApp.*')
MFPGEN = rdFingerprintGenerator.GetMorganGenerator(1)
ao = rdFingerprintGenerator.AdditionalOutput()
ao.AllocateBitInfoMap()
ao.AllocateAtomToBits()
p = Chem.MolFromSmarts("[#0]") # dummy atom
srflag = Chem.SanitizeFlags.SANITIZE_SYMMRINGS



def mol_to_pairs(m):
    """
    function that fractures every bond and reports the two ECFP2
    (including dummy) at the fracture point.
    """
    id_pairs = []
    ri_full = m.GetRingInfo()
    ar = ri_full.AtomRings()
    for b in m.GetBonds():
        bidx = [b.GetBeginAtomIdx(),b.GetEndAtomIdx()]
        newmol=Chem.FragmentOnBonds(m,[b.GetIdx()])
        try:
            if b.IsInRing():
                Chem.SanitizeMol(newmol,sanitizeOps=srflag)
                ri = newmol.GetRingInfo() #reset ringinfo trick
                for ring in ar:
                    for idx in ring:
                        ri.AddRing((idx,), (0,))
            else:
                Chem.SanitizeMol(newmol)
            MFPGEN.GetSparseFingerprint(newmol,fromAtoms=bidx,additionalOutput=ao)
            id_pairs.append(tuple(sorted([ao.GetAtomToBits()[idx][1] for idx in bidx])))
        except Exception as e:
            pass #silent for now
    return id_pairs


def get_profile_for_mols(suppl,profile_name,size=1024):
    all_pairs = [mol_to_pairs(m) for m in suppl if m]
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
        except:
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

    
    args = vars(parser.parse_args())
    suppl = Chem.SmilesMolSupplier(args["input"],titleLine=False)
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
        get_profile_for_mols(suppl,args["profile"],size=args["size"])
