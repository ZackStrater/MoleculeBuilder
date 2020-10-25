from fragments_building_library import heterocycles, functionalized_arenes, functional_groups, hydrocarbons, amines, amino_acids
import random
from Smiles_to_Structure import MoleculeStructure, convert_to_structure
from Strucutre_to_Smiles import convert_to_smiles, Bond

def build_molecule(molecule, frag_num, *fragment_libraries):   #TODO remove molecule from arg and just include in func
    combined_frag_library = {}
    for library in fragment_libraries:
        combined_frag_library.update(library)
    combined_frag_list = list(combined_frag_library)
    frag_list = [MoleculeStructure() for _ in range(frag_num)]
    # list of empty molecular structures
    # structures will get filled with random fragments and then fragments will get pieced together
    frag_name = random.choice(combined_frag_list)
    # choose random fragment key from fragment list
    print(frag_name)
    frag = combined_frag_library[frag_name]
    # frag = string of that fragment
    convert_to_structure(frag_list[0], frag)
    # convert the seed fragment to structure using the first empty MolecularStructre(), where frag = fragment smiles
    for atom in frag_list[0].atom_list:
        molecule.atom_list.append(atom)
    for c in range(1, frag_num):
        print_name = random.choice(combined_frag_list)
        print(print_name)
        frag = combined_frag_library[print_name]
        convert_to_structure(frag_list[c], frag)
        frag_atom = random.choice([atom for atom in frag_list[c].atom_list if atom.can_bond])
        if frag_atom.heteroatom:
            molecule_atom = random.choice([atom for atom in molecule.atom_list if atom.can_bond and not atom.heteroatom])
        else:
            molecule_atom = random.choice([atom for atom in molecule.atom_list if atom.can_bond])
        frag_atom.bonded_to.append(Bond(molecule_atom, 1))
        molecule_atom.bonded_to.append(Bond(frag_atom, 1))
        if len(frag_atom.bonded_to) > 2:
            frag_atom.can_bond = False
        if len(molecule_atom.bonded_to) > 2:
            molecule_atom.can_bond = False
        for atom in frag_list[c].atom_list:
            molecule.atom_list.append(atom)

    print(convert_to_smiles(molecule, molecule.atom_list[0]))


# TODO add option for two point frag add, which would add a new ring


for i in range(100):
    build_molecule(MoleculeStructure(), 5, heterocycles, functionalized_arenes, functional_groups, hydrocarbons)


