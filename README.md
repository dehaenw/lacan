# LACAN
LACAN filter: Leveraging adjacent co-ocurrence of atomic neighborhoods for molecular filtering

> "All sorts of things in the world behave like mirrors"
> -Jacques Lacan

Some molecular fragments are common, but they have the tendency not to occur together. For example, alkyloxy radicals are frequent motifs in medicinal chemistry datasets, whereas the linkage of both radicals into a peroxide is rather uncommon. Likewise, halides and amines are some of the most commonly occurring atomic neighborhoods, and yet their pairing results in the unstable and toxic haloamine motif. We apply this concept using co-occurences of ECFP2 like atomic neighborhoods at the bond interface, and leverage co-occurence patterns to construct a molecular filter that highlights uncommon linkages.

## Installation

clone this repo, activate your environment, navigate to root dir and run:

```
pip install .
```

## Basic usage: Localizing problem bonds

import lacan and inspect a molecule by running the following commands:

```python
from lacan import lacan
from rdkit import Chem
p = lacan.load_profile("chembl")
m = Chem.MolFromSmiles("c1ccccc1CCN(OCCc1occc1)")
score,info = lacan.score_mol(m,p)
print(info["bad_bonds"])
```

which will output a dictionary with an entry for every bond in the molecule. Currently the filter is binary, so the score is 1 if the molecule passes the filter and 0 if it doesn't. The problem bonds output
follow rdkit bond numbering which means we can visualize problem bonds in our
molecules easily as follows:

```python
from rdkit.Chem import Draw
d = Draw.MolToImage(m,highlightBonds=info["bad_bonds"])
display(d)
```

giving the following result:

![image](https://github.com/user-attachments/assets/5758aace-c6aa-4aaf-a04a-a58a31fe48af)


This correctly identified the N-O linkage as problematic.

## Breeding molecules

This filter enables us to recombine fragments and filter out linkages that are rare in the reference set. Lacan has a "breeding" or crossover functionality where two molecules get fragmented and recombined. By subjecting the recombinations to LACAN filter we can retain only decent looking "median molecules".

example:
```python
from lacan import breed
from rdkit import Chem
from rdkit.Chem import Draw

m1 = Chem.MolFromSmiles("c1cc(ccc1[C@@H]2CCNC[C@H]2COc3ccc4c(c3)OCO4)F")
m2 = Chem.MolFromSmiles("CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F")
median_molecules = breed.breed(m1,m2,p,nmols=9)
```

this outputs the following molecules that are "in between" its parents fluoxetine and sertraline:
```python
d = Draw.MolsToGridImage(median_molecules)
display(d)
```
![image](https://github.com/user-attachments/assets/e6609b81-21ca-4f9a-9c3c-3fe15cbc38d8)


## Building a profile

If you want to build a custom profile using your own reference data set, this can be done through the LACAN cli as follows

`python lacan.py -i your_dataset_here.smi -m profile -p my_new_profile`

This will create a pickled profile in the data folder which you can then invoke using:

`p = lacan.load_profile("my_new_profile")`
