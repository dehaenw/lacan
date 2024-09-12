from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions, inchi
from lacan import lacan
import argparse
import random
from itertools import combinations

RDLogger.DisableLog('rdApp.*')
combine_frags = rdChemReactions.ReactionFromSmarts('[*:0][#0:1].[*:2][#0:3]>>[*:0][*:2].[*:1][*:3]')

def fragment_molecule(mol,n=3):
    """
    this function takes a molecule and cuts it n times. the cuts that are
    retained are those that have n fragments with one dummy and 1 fragment
    with n dummies ("substituents and cores")
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
    
def crossover_fragments(s1,s2,c1,c2,profile,nmols=10,randomseed=123,max_steps=50000,hacmin = 0,hacmax= 30, min_ratio=0.25, max_ratio=0.75):
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
            if lacan.score_mol(cmol,profile)[0]>0:
                if hacmin < len(cmol.GetAtoms()) < hacmax:
                    if min_ratio < ratio < max_ratio:
                        ik = inchi.MolToInchiKey(cmol)
                        if  ik not in inchis:
                            mols.append(cmol)
                            inchis.append(ik)
        except Exception as e:
            print(e)
            pass
        steps += 1
    print(len(mols)/steps)
    return mols
    
def breed(m1,m2,profile,nmols=10,cuts=3,hacrange=(0.8,1.2),interprange=(0.3,0.7)):
    s1,c1 = fragment_molecule(m1,cuts)
    s2,c2 = fragment_molecule(m2,cuts)
    n1 = len(m1.GetAtoms())
    n2 = len(m2.GetAtoms())
    if len(s1) == 0:
        print("didnt find a way to cut >2 times, reverting to twice")
        s1,c1 = fragment_molecule(m1,2)
        s2,c2 = fragment_molecule(m2,2)
    mols = crossover_fragments(s1,s2,c1,c2,profile,nmols,hacmin = int(hacrange[0]*min(n1,n2)),hacmax = int(hacrange[1]*max(n1,n2)),min_ratio=interprange[0],max_ratio=interprange[1])
    return mols


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lacan Breed CLI")
    parser.add_argument(
        "-i1",
        "--input1",
        type=str,
        help="input file. should be a smiles",
        required=True,
    )
    parser.add_argument(
        "-i2",
        "--input2",
        type=str,
        help="input file 2. should be a smiles",
        required=True,
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
        "-n",
        "--nmols",
        type=int,
        default=10,
        help="amount of mols to output",
        required=False,
    )

    
    args = vars(parser.parse_args())
    m1 = Chem.MolFromSmiles(args["input1"])
    m2 = Chem.MolFromSmiles(args["input2"])
    
    
    
    p = lacan.load_profile(args["profile"])
    for m in breed(m1,m2,p,args["nmols"]):
        print(Chem.MolToSmiles(m))
