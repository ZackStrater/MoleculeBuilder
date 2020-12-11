
import pandas as pd
from Smiles_to_Structure import convert_to_structure, MoleculeStructure
from scipy import stats
import os
from math import log

path = "."

df1 = pd.read_csv("LogP_fragments_data.csv")

n = []
o = []
c = []
nc = []
oc = []
nmw = []
omw = []

for smil in df1["smiles"]:
    mol = convert_to_structure(MoleculeStructure(), smil)
    carbon_num = 0
    nitrogen_num = 0
    oxygen_num = 0
    mw = mol.mw()
    for atom in mol.atom_list:
        if atom.symbol == "C" or atom.symbol == "c":
            carbon_num += 1
        if atom.symbol == "N" or atom.symbol == "n":
            nitrogen_num += 1
        if atom.symbol == "O" or atom.symbol == "o":
            oxygen_num += 1
    n.append(nitrogen_num)
    o.append(oxygen_num)
    c.append(carbon_num)
    try:
        nc.append(nitrogen_num/carbon_num)
    except(TypeError, ZeroDivisionError):
        nc.append("NA")

    try:
        oc.append(oxygen_num/carbon_num)
    except(TypeError, ZeroDivisionError):
        oc.append("NA")

    nmw.append(nitrogen_num/log(mw))
    omw.append(oxygen_num/log(mw))

# z score = (value - mean)/stdev
n = stats.zscore(n)
o = stats.zscore(o)
c = stats.zscore(c)

print(f"n{n}")
print(f"o{o}")
print(f"c{c}")
print(f"nc{nc}")
print(f"oc{oc}")
print(f"nmw{nmw}")
print(f"omw{omw}")

df2 = pd.DataFrame({"Nitrogen(Zscore)": n,
                    "Oxygen(Zscore)": o,
                    "Carbon(Zscore": c,
                    "Nitrogen/Carbon": nc,
                    "Oxygen/Carbon": oc,
                    "Nitrogen/MW": nmw,
                    "Oxygen/MW": omw,
                    })

df3 = pd.concat([df1, df2], axis=1)
updated_LogP_fragments_data = os.path.join(path, "updated_LogP_fragments_data.csv")
df3.to_csv(updated_LogP_fragments_data, index=False)

