from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions, inchi
from lacan import lacan
import random,argparse

RDLogger.DisableLog('rdApp.*')

mutate_smarts = {"addC":"[H1,H2,H3:0]>>[*:0][CH3]",
                 "addO":"[H1,H2,H3:0]>>[*:0][OH]",
                 "addCO":"[CH2:0]>>[*:0]=[O]",
                 "addN":"[H1,H2,H3:0]>>[*:0][NH2]",
                 "addF":"[H1,H2,H3:0]>>[*:0][F]",
                 "addCl":"[H1,H2,H3:0]>>[*:0][Cl]",
                 "addBr":"[H1,H2,H3:0]>>[*:0][Br]",
                 "contractAroO":"[ar6:0]:[cH,nH0;r6:1]:[cH,nH0;r6:2]:[ar6:3]>>[*:0]:[o]:[*:3].[*:1][*:2]",
                 "contractAroS":"[ar6:0]:[cH,nH0;r6:1]:[cH,nH0;r6:2]:[ar6:3]>>[*:0]:[s]:[*:3].[*:1][*:2]",
                 "contractAroNH":"[ar6:0]:[cH,nH0;r6:1]:[cH,nH0;r6:2]:[ar6:3]>>[*:0]:[nH]:[*:3].[*:1][*:2]",
                 "expandAroCC":"[ar5:0]:[nH,o,s;r5:1]:[ar5:2]>>[*:0]:[cH]:[cH]:[*:2].[*:1]",
                 "expandAroCN":"[ar5:0]:[nH,o,s;r5:1]:[ar5:2]>>[*:0]:[cH]:[n]:[*:2].[*:1]",
                 "insertC":"[*:0]-[*:1]>>[*:0]-[CH2]-[*:1]",
                 "insertN":"[*:0]-[*:1]>>[*:0]-[NH]-[*:1]",
                 "insertO":"[*:0]-[*:1]>>[*:0]-[O]-[*:1]",
                 "insertS":"[*:0]-[*:1]>>[*:0]-[S]-[*:1]",
                 "replaceC":"[!C;A;d4,d3,d2,d1:0]>>[C:0]",
                 "replaceN":"[!N;!$([CH0]);d3,d2,d1;A:0]>>[N:0]",
                 "replaceO":"[!O;!$([CH0]);$([*](-[*])(-[*]));d2,d1;A:0]>>[O:0]",
                 "replaceS":"[!S;!$([CH0]);$([*](-[*])(-[*]));d2,d1;A:0]>>[S:0]",
                 "aroCtoN":"[cH:0]>>[nH0:0]",
                 "aroNtoC":"[nH0X2:0]>>[cH:0]",
                 "openring":"[R]@!:[R]>>([*:0].[*:1])",
                 "close3ring":"[!R;H1,H2,H3:0][*:1][!R;H1,H2,H3:2]>>[*:0]1~[*:1]~[*:2]1",
                 "close4ring":"[!R;H1,H2,H3:0][*:1][*:2][!R;H1,H2,H3:3]>>[*:0]1~[*:1]~[*:2]~[*:3]1",
                 "close5ring":"[!R;H1,H2,H3:0][*:1][*:2][!a:3][!R;H1,H2,H3:4]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]1",
                 "close6ring1":"[!R;H1,H2,H3:0][*:1][*:2][!a:3][!a:4][!R;H1,H2,H3:5]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]~[*:5]1",
                 "close6ring2":"[!R;H1,H2,H3:0][!a:1][*:2][*:3][!a:4][!R;H1,H2,H3:5]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]~[*:5]1",
                 "arofuse5":"[aH1:0]:[a:1]-[*:2][*:3]-[!R;H1,H2,H3:4]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]1",
                 "arofuse6na1":"[aH1:0]:[a:1]-[*:2][*:3][!a:4]-[!R;H1,H2,H3:5]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]~[*:5]1",
                 "arofuse6na2":"[aH1:0]:[a:1]-[!a:2][*:3][*:4]-[!R;H1,H2,H3:5]>>[*:0]1~[*:1]~[*:2]~[*:3]~[*:4]~[*:5]1",
                 "bond2to3":"[!R;H1,H2:0]=[!R;H1,H2:1]>>[*:0]#[*:1]",
                 "bond1to2":"[H1,H2,H3:0]-[H1,H2,H3:1]>>[*:0]=[*:1]",
                 "bond2to1":"[A:0]=[A:1]>>[*:0]-[*:1]",
                 "bond3to2":"[*:0]#[*:1]>>[*:0]=[*:1]",
                 "deleteD1":"[*:0][d1:1]>>[*:0].[*:1]",
                 "deleteD2":"[*:0][d2;A:1][*:2]>>[*:0][*:2].[*:1]"}
                 
mutate_ops = {name:rdChemReactions.ReactionFromSmarts(mutate_smarts[name]) for name in mutate_smarts}


def apply_mutations(mol,p,score_threshold,mode="all"):
    mutate_prods = []
    if mode=="all":
        ops = mutate_ops
    elif mode=="random":
        ops = random.sample(list(mutate_ops.keys()),1)[0]
    else:
        print("this mode doesnt exist")
    for op in mutate_ops:
        rxn = mutate_ops[op]
        prods = rxn.RunReactants((mol,))
        if len(prods)>0:
            for prod in prods:
                try:
                    Chem.SanitizeMol(prod[0])
                    if prod[0]:
                        mutate_prods.append(prod[0])
                except Exception as e:
                    print(op,"got some exception",e)

    filtered_prods = []
    iks = []
    for m in mutate_prods:
        score,info = lacan.score_mol(m,p)   
        if score>score_threshold: 
            ik = inchi.MolToInchiKey(m)
            if ik not in iks:
                filtered_prods.append(m)
                iks.append(ik)
    return filtered_prods
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lacan Mutate CLI")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="input structure. should be a smiles",
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
        "-m",
        "--mode",
        type=str,
        default="all",
        help="name of mode to run",
        required=False,
    )
    
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.8,
        help="threshold for lacan filter score",
        required=False,
    )
    

    
    args = vars(parser.parse_args())
    p = lacan.load_profile(args["profile"])
    mol = Chem.MolFromSmiles(args["input"])
    for m in apply_mutations(mol,p,args["threshold"],args["mode"]):
        print(Chem.MolToSmiles(m))
