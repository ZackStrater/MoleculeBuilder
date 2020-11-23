
import re
smiles_string = "C(H)CC[H]CC([H])"




corrected_smiles_string = re.sub(r"\(\[H\]\)|\(H\)", "", smiles_string)


print(corrected_smiles_string)