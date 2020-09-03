
import re
from molecule_builder import Atom


smiles_string = "[R]CC(=O)NCCc1c[nH]c2ccc(OC)cc12"
atom_map = []

for match in re.finditer(r"H|B|C|N|O|P|S|Si|F|Cl|Br|I|\[R\]|b|c|n|o|p|s", smiles_string):
        print(match.group())
        atom_map.append((match.start(), Atom(match.group())))


bond_map = re.findall(r"(?:[A-Za-z])([^A-Za-z]*)", smiles_string)

print(atom_map)
print(bond_map)