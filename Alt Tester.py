import re


smiles_string = "[H][C@@]12N(C3=C([C@@]14CCN5CC=C[C@@]([C@H]([C@@]2(C(OC)=O)O)OC(C)=O)([C@@]45[H])CC)C=C([C@@]6(C(OC)=O)C[C@@H]7CN(CCC8=C6NC9=CC=CC=C89)C[C@](CC)(C7)O)C(OC)=C3)C"

smiles_string1 = re.sub(r"[\[\]H@]", "", smiles_string)

print(smiles_string1)