from rdkit import Chem
from lacan import lacan
from lacan import mutate, breed, decompose
import random
import csv
import os
from rdkit import DataStructs
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdFingerprintGenerator
import multiprocessing


MFPGEN = rdFingerprintGenerator.GetMorganGenerator(2,fpSize=4096)
rxn1 = rdChemReactions.ReactionFromSmarts("[*:0]-[#0:1].[*:2]-[#0:3]>>[*:0]-[*:2].[*:1]-[*:3]")
rxn2 = rdChemReactions.ReactionFromSmarts("[*:0]=[#0:1].[*:2]=[#0:3]>>[*:0]=[*:2].[*:1]=[*:3]")

entries = []
this_dir, this_filename = os.path.split(__file__)
with open( os.path.join(this_dir, f"data/rls.csv"), newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    for e in reader:
        if e[0] != "smiles":
            if int(e[1]) > 200: #make this customizable
                entries.append([e[0],int(e[1]),int(e[2]),e[3],e[4]])
FRAGCOUNT = {"Sub":0,"Linker":0,"Ring":0}
for e in entries:
    FRAGCOUNT[e[3]] += e[1]


def generate_random_molecule(fragments=[],fragcorpus=entries,needed_sample=[],bondstring="",mol=None):
    if len(fragments)>0:
        new_needed_sample = []
        new_bondstring = ""
        es = []
        for i,cns in enumerate(needed_sample):
            fragtype = random.choices(cns,[FRAGCOUNT[ft] for ft in cns])[0]
            fe = [e for e in fragcorpus if e[3] == fragtype and bondstring[i]==e[4][0]]
            es.append(random.choices([k for k in fe],[k[1]**0.666 for k in fe])[0]) #change or make customizable
        for i,e in enumerate(es):
            fragments.append(e[0])
            if bondstring[i]=="-":
                mol = rxn1.RunReactants((mol,Chem.MolFromSmiles(e[0])))[0][0]
            else:
                mol = rxn2.RunReactants((mol,Chem.MolFromSmiles(e[0])))[0][0]
            Chem.SanitizeMol(mol)
            if e[3] == "Ring" and e[2]>1:
                new_needed_sample += (e[2]-1)*[["Sub","Linker"]]
                new_bondstring += e[4][1:]
            if e[3] == "Linker":
                new_needed_sample += (e[2]-1)*[["Ring"]]
                new_bondstring += e[4][1:]
        if len(new_needed_sample)>0:
            mol,fragments = generate_random_molecule(fragments,fragcorpus,new_needed_sample,new_bondstring,mol)
    else:
        e = random.choices([k for k in fragcorpus],[k[1]**0.666 for k in fragcorpus])[0]
        mol = Chem.MolFromSmiles(e[0])
        if e[3] == "Linker":
            mol,fragments = generate_random_molecule([e[0]],fragcorpus,e[2]*[["Ring"]],e[4],mol)
        elif e[3] == "Sub":
            mol,fragments = generate_random_molecule([e[0]],fragcorpus,[["Ring"]],e[4],mol)
        else:
            mol,fragments = generate_random_molecule([e[0]],fragcorpus,e[2]*[["Sub","Linker","Ring"]],e[4],mol)
    return mol,fragments

def generate_filtered_molecule(profile,fragcorpus=entries,threshold=0.5,seed=None,min_atoms=14,max_atoms=35):
    score = 0
    if seed:
        random.seed(seed)
    while score < threshold:
        mol, _ = generate_random_molecule(fragcorpus=fragcorpus)
        score, _ = lacan.score_mol(mol, profile)
        if min_atoms<len(mol.GetAtoms())<max_atoms:
            pass
        else:
            score=0
    return mol
    
def generate_filtered_molecules(profile,fragcorpus=entries,threshold=0.5,seed=None,min_atoms=14,max_atoms=35,n_jobs=1,n_molecules=10):
    if n_jobs==1:
        mols = [generate_filtered_molecule(**kwargs) for i in range(n_molecules)]
    else:
        if n_jobs==-1:
            n_jobs = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=n_jobs)
        mols = pool.starmap(generate_filtered_molecule, [(profile,fragcorpus,threshold,seed+823848*i,min_atoms,max_atoms) for i in range(n_molecules)])
        pool.close()
        pool.join()
    return mols
    

def get_corpus_from_mols(mols,fragcorpus=entries,ratio=1):
    custom_entries = decompose.get_corpus(mols) # this can take a while ...
    custom_fragcount = sum([e[1] for e in custom_entries])
    fragcount = sum([e[1] for e in custom_entries])
    mult = fragcount/custom_fragcount*ratio
    fragcounts = {e[0]:e[1] for e in fragcorpus}
    custom_fragcounts = {e[0]:e[1] for e in custom_entries}
    merged_entries = []
    for e in custom_entries:
        if e[0] in fragcounts:
            merged_entries.append([e[0]]+[int(e[1]*mult)+fragcounts[e[0]]]+e[2:])
        else:
            merged_entries.append([e[0]]+[int(e[1]*mult)]+e[2:])
    for e in fragcorpus:
        if e[0] not in custom_fragcounts:
            merged_entries.append(e)
    return merged_entries


def generate_optimized_molecules(scoring_function,profile,
              seed=123,
              startN=50,
              generations=10,
              popsize=20,
              maxmutations=99999,
              win_threshold=0.8,
              sim_threshold=0.45,
              higher_is_better=True,
              quiet=False,
              fragcorpus=entries):
    random.seed(seed)
    winners = []
    winnerfps = []
    scores = []
    starting_mols = generate_filtered_molecules(profile,fragcorpus=fragcorpus,min_atoms=15,n_molecules=startN,n_jobs=-1,seed=seed)
    rs = scoring_function(starting_mols)
    if higher_is_better:
        win_threshold *= -1
        rs = [-x for x in rs]
    for i,starting_mol in enumerate(starting_mols):
        if rs[i]<win_threshold:
            winners.append((Chem.MolToSmiles(starting_mol),rs[i]))
            winnerfps.append(MFPGEN.GetFingerprint(starting_mol))
        else:
            scores.append((Chem.MolToSmiles(starting_mol),rs[i])) 
    #mutate best N
    for j in range(generations):
        if seed:
            seed *= 6543 # to ensure different mols
        scores = sorted(scores,key=lambda x:x[1])[:len(scores)] 
        if len(winnerfps)>0:
            sims = [max(DataStructs.BulkTanimotoSimilarity(
                MFPGEN.GetFingerprint(Chem.MolFromSmiles(s[0])),winnerfps)) for s in scores]
            scores = [score for i,score in enumerate(scores) if sims[i]<sim_threshold]
        if not quiet:
            print(len(winners),"winners")
        if len(scores)==0:
            bestNmols = ["CC","CCC"]
        else:
            if not quiet:
                print("best score so far",scores[0][1])
            bestNmols = [scores[i][0] for i in range(min(popsize,len(scores)))]
        if len(scores)>len(bestNmols):
            scores = scores[len(bestNmols):]
        else:
            if not quiet:
                print("generating new mols")
            starting_mols = generate_filtered_molecules(profile,fragcorpus=fragcorpus,min_atoms=15,n_molecules=startN,n_jobs=-1,seed=seed)
            rs = scoring_function(starting_mols)
            if higher_is_better:
                rs = [-x for x in rs]
            for i,starting_mol in enumerate(starting_mols):
                if rs[i]<win_threshold:
                    winners.append((Chem.MolToSmiles(starting_mol),rs[i]))
                    winnerfps.append(MFPGEN.GetFingerprint(starting_mol))
                else:
                    scores.append((Chem.MolToSmiles(starting_mol),rs[i]))
            scores = sorted(scores,key=lambda x:x[1])[:len(scores)]
        mols = [Chem.MolFromSmiles(smi) for smi in bestNmols]
        all_mutated_mols = mutate.apply_mutations_mols(mols,profile,0.8)
        all_mutated_mols += generate_filtered_molecules(profile,fragcorpus=fragcorpus,min_atoms=15,n_molecules=startN,n_jobs=-1,seed=seed)
        rs = scoring_function(all_mutated_mols)
        if higher_is_better:
            rs = [-x for x in rs]       
        for i,mmol in enumerate(all_mutated_mols):
            msmi = Chem.MolToSmiles(mmol)
            if msmi not in [k[0] for k in scores]:
                if rs[i]<win_threshold:
                    mfp = MFPGEN.GetFingerprint(mmol)
                    if len(winnerfps)>0:
                        if max(DataStructs.BulkTanimotoSimilarity(mfp,winnerfps))<sim_threshold:
                            winners.append((msmi,rs[i]))
                            winnerfps.append(MFPGEN.GetFingerprint(mmol))
                    else:
                        winners.append((msmi,rs[i]))
                        winnerfps.append(MFPGEN.GetFingerprint(mmol))
                else:
                    scores.append((msmi,rs[i]))
    winners = sorted(winners,key=lambda x:x[1])
    if higher_is_better:
        winners = [(c[0],-1*c[1]) for c in winners]
    return winners
    
    
def next_population(population):
    mutated_mols = mutate.apply_mutations_mols(mols,profile,0.8)
    replaced_mols = random_replacement_mols(mols,profile,0.8)# 
    breed_mols = breed.cross_breed_mols(mols,profile,0.8)#
        
