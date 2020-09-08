from Smiles_to_Structure import convert_to_structure, MoleculeStructure
from fragments import heterocycles
from collections import Counter

def fragmentize(smiles_string, fragment_string):  #specify fragment list

    total_molecule = convert_to_structure(MoleculeStructure(), smiles_string)
    fragment = convert_to_structure(MoleculeStructure(), fragment_string)

    def find_atom():
        for ele in ["S", "P", "B", "Si", "I", "Br", "Cl", "F", "O", "N", "C"]:
            for atom in fragment.atom_list:
                if atom.symbol == ele:
                    return atom
    base_atom = find_atom()

    def abbr_bond(bond):
        return bond.bond_code, bond.atom.symbol

    base_atom_bonds = [abbr_bond(bond) for bond in base_atom.bonded_to]
    count_base_atom_bonds = Counter(base_atom_bonds)
    starting_atoms = []

    def find_starter_atoms(atom, fragment_counter, atom_list):
        if atom.symbol == base_atom.symbol:
            bond_list = [abbr_bond(bond) for bond in atom.bonded_to]
            count_bond_list = Counter(bond_list)
            for key in fragment_counter:
                if key not in count_bond_list or fragment_counter[key] > count_bond_list[key]:
                    return
            print(atom.symbol)
            print(count_bond_list)
            atom_list.append(atom)

    for atom in total_molecule.atom_list:
        find_starter_atoms(atom, count_base_atom_bonds, starting_atoms)




    #start substructure search
    # start based atom see if bonds match
    # go outward and see if bonding pattern matches
    # then go back and do internal bonding

fragmentize("O=C1OC(C2=C(C(N3C)=CC(N4C)=NC3=O)C(C5=NC6=C4C=C(C7=CC=C8C(OC=C8)=C7)C=C6N=C5)=C9N=COC9=C2)CC1C%10=CC(C%11=NC(C%12=C%13C(C(C%14=C%15C=CN=NC%15=CC=C%14)=CC=N%13)=CC(C%16=NC=CC=N%16)=C%12)=NN%11C)=NO%10", "C1=NCCC1")
