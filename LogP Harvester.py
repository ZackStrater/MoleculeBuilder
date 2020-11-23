import re
from cas_name_logp_data import Kow_data
from Smiles_CAS_name import Smiles_data
import pandas as pd
from find_fragments import fragmentize
from fragments_library import biomolecules, peptide_amino_acids, heterocycles, arenes, functional_groups, hydrocarbons

# list of tuples with (CAS, Chemical Name, LogP value)
CAS_LogP = re.findall(r"(\d{6}-\d{2}-\d)(.+?)(\D\d{1}\.\d{2}\D)", Kow_data)

# dict with CAS as keys and properly formatted LogP as values (also gets rid of duplicates)
CAS_LogP_dict = {t[0]: t[2][1:-1] for t in CAS_LogP}
print(len(CAS_LogP_dict))

# list of tuples with (CAS, Chemical Name, smiles string)
CAS_smiles = re.findall(r"(\d{6}-\d{2}-\d)(.+?)(\s\S+\s\S\s)(?=\d{6}-\d{2}-\d)", Smiles_data)

# dict with CAS as keys and properly formatted smiles string as values
CAS_smiles_dict = {t[0]: t[2][1:-3] for t in CAS_smiles}
print(len(CAS_smiles_dict))

# add the CAS_LogP data to the combined list, excluding any CAS numbers that
# don't show up in the CAS_smiles data.  the values are put inside a list
combined_dict = {k: [v] for k, v in CAS_LogP_dict.items() if k in CAS_smiles_dict}

print(len(combined_dict))

# append the smiles strings to the value (which is a list containing the logP)
# of each CAS key
for k, v in CAS_smiles_dict.items():
    if k in combined_dict:
        combined_dict[k].append(v)

combined_dict = {k: v for k, v in combined_dict.items() if "[" not in v[1]}
print(len(combined_dict))

CAS_data = [key for key in combined_dict]
LogP_data = [float(v[0]) for k, v in combined_dict.items()]  # convert values to numerical values
smiles_data = [v[1] for k, v in combined_dict.items()]

df = pd.DataFrame({"CAS": CAS_data, "LogP": LogP_data, "smiles": smiles_data})
#df["fragments"] = df["smiles"].map(lambda x: fragmentize(x, biomolecules, peptide_amino_acids, heterocycles, arenes, functional_groups, hydrocarbons))



import time

from timeit import default_timer as timer


starter = timer()
start = time.process_time()

print("time:")
print(time.process_time() - start)

for e in smiles_data[0:100]:
    print(e)
    fragmentize(e, biomolecules, peptide_amino_acids, heterocycles, arenes, functional_groups, hydrocarbons)

ender = timer()
print("timeit:")
print(ender - starter)