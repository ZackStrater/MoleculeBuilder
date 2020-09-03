import re
from cas_name_logp_data import Kow_data
from Smiles_CAS_name import Smiles_data

Chem_info1 = re.findall(r"(\d{6}-\d{2}-\d)(.+?)(\D\d{1}\.\d{2}\D)", Kow_data)
# tuples with (CAS, Chemical Name, LogP value) - 12534 matches

CAS = [ele[0] for ele in Chem_info1]
CAS = list(dict.fromkeys(CAS))
# removing duplicates
# all CAS from cheminfo1 - 12534 matches
print("CAS")
print(len(CAS))

Chemical_Names = [ele[1] for ele in Chem_info1]

LogP = [ele[2] for ele in Chem_info1]

Chem_info2 = re.findall(r"(\d{6}-\d{2}-\d)(.+?)(\s\S+\s\S\s)", Smiles_data)
# tuples with (CAS, Chemical Name, Smiles Formula - 107108 matches

print(len(Chem_info2))

Corrected_Chem_info2 = [ele for ele in Chem_info2 if ele[0] in CAS]
# Pairing down Cheminfo2 to filter out nonmatches wrt data from CAS1/cheminfo1- 12257 matches
print("Corrected_Chem_info2")
print(len(Corrected_Chem_info2))

Corrected_CAS = [ele[0] for ele in Corrected_Chem_info2]
print("corrected CAs")
print(len(Corrected_CAS))

CAS = [ele for ele in CAS if ele in Corrected_CAS]
print("CAS")
print(len(CAS))


Smiles_formulas = [ele[2] for ele in Corrected_Chem_info2]
# Smiles formulas; filtered out non-matches wrt data from CAS1/Chem_info1 - 12257 matches

Smiles_formulas = [ele[1:-3] for ele in Smiles_formulas]
print(Smiles_formulas)
print(len(Smiles_formulas))
# reformatting

bracket_count = 0
for ele in Smiles_formulas:
    if "[" in ele:
        bracket_count += 1
print(bracket_count)
#  497 strings with brackets



# make dictionary - > CAS - > name
# make dictionary - > CAS _-> SMILES
# make dictionary -> CAS -> Log P


