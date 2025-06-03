from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions, inchi
#from lacan import lacan
import random,argparse,csv
from collections import Counter

RDLogger.DisableLog('rdApp.*')

decompose1 = rdChemReactions.ReactionFromSmarts("[!R;!#0:0]-[R:1]>>[!R:0]-[#0].[R:1]-[#0]")
decompose2 = rdChemReactions.ReactionFromSmarts("[!R;!#0:0]=[R:1]>>[!R:0]=[#0].[R:1]=[#0]")
ringdummy = Chem.MolFromSmarts("[#0]~[R]")
nonringdummy = Chem.MolFromSmarts("[#0]~[!R]")
singledummy = Chem.MolFromSmarts("[#0;$([*]-[*])]")
doubledummy = Chem.MolFromSmarts("[#0;$([*]=[*])]")

def get_bonds_string(smi):
    mol = Chem.MolFromSmiles(smi)
    bonds = []
    for db in mol.GetSubstructMatches(singledummy):
        bonds.append(("-",db[0]))
    for db in mol.GetSubstructMatches(doubledummy):
        bonds.append(("=",db[0]))
    bonds = sorted(bonds,key=lambda x:x[1])
    return "".join([b[0] for b in bonds])

def decompose_molecule(mol):
    frags = [mol]
    final_frags = []
    while len(frags)>0:
        prods = decompose1.RunReactants((frags[0],))
        prods += decompose2.RunReactants((frags[0],))
        if len(prods) > 0:
            try:
                frags.pop(0)
                Chem.SanitizeMol(prods[0][0])
                Chem.SanitizeMol(prods[0][1])
                frags += [prods[0][0],prods[0][1]]
            except:
                frags = []
                print("problem. this should never happen")
        else:
            final_frags.append(frags[0])
            frags.pop(0)
    rings = []
    linkers = []
    subs = []
    if len(final_frags)>1:
        for frag in final_frags:
            rec = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(frag)))
            if frag.HasSubstructMatch(ringdummy):
                rings.append(rec)
            else:
                sm = len(frag.GetSubstructMatches(nonringdummy))
                if sm>1:
                    linkers.append(rec)
                elif sm == 1:
                    subs.append(rec)
    return rings,linkers,subs
    
def get_all_decompositions(suppl):
    all_rings = []
    all_linkers = []
    all_subs = []
    for mol in suppl:
        if mol:
            Chem.RemoveStereochemistry(mol)
            r, l, s = decompose_molecule(mol)
            all_rings += r
            all_linkers += l
            all_subs += s
    return all_rings,all_linkers,all_subs
    
def save_frags(r,l,s,minN=10):
    cr = dict(Counter(r))
    cl = dict(Counter(l))
    cs = dict(Counter(s))
    entries = []
    for k in cr:
        if cr[k] > minN:
            entries.append([k,cr[k],k.count("*"),"Ring",get_bonds_string(k)])
    for k in cl:
        if cl[k] > minN:
            entries.append([k,cl[k],k.count("*"),"Linker",get_bonds_string(k)])
    for k in cs:
        if cs[k] > minN:
            entries.append([k,cs[k],k.count("*"),"Sub",get_bonds_string(k)])

    entries=sorted(entries,key=lambda x:-x[1])
    with open('rls_test.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["smiles","occurrence","degree","type","bonds"])
        for entry in entries:
            writer.writerow(entry)
    return
    
def get_corpus(mols):
    r=[]
    l=[]
    s=[]
    for mol in mols:
        if mol:
            r0,l0,s0 = decompose_molecule(mol)
            r+=r0
            l+=l0
            s+=s0
    cr = dict(Counter(r))
    cl = dict(Counter(l))
    cs = dict(Counter(s))
    entries = []
    for k in cr:
        entries.append([k,cr[k],k.count("*"),"Ring",get_bonds_string(k)])
    for k in cl:
        entries.append([k,cl[k],k.count("*"),"Linker",get_bonds_string(k)])
    for k in cs:
        entries.append([k,cs[k],k.count("*"),"Sub",get_bonds_string(k)])
    entries=sorted(entries,key=lambda x:-x[1])
    return entries
       
    
if __name__ == "__main__":
    suppl = Chem.SmilesMolSupplier("/home/wim/Downloads/chembl_35_cleaned.csv")
    r,l,s = get_all_decompositions(suppl)
    save_frags(r,l,s)
