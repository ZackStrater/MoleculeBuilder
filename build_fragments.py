from fragments import fragments
import random
from Smiles_to_Structure import MoleculeStructure, convert_to_structure
from molecule_builder import convert_to_smiles, Bond


def build_molecule(molecule, frag_num):
    frag_list = [MoleculeStructure() for _ in range(frag_num)]
    frag = fragments[random.choice(list(fragments))]
    convert_to_structure(frag_list[0], frag)
    for atom in frag_list[0].atom_list:
        molecule.atom_list.append(atom)
    for i in range(1, frag_num):
        frag = fragments[random.choice(list(fragments))]
        convert_to_structure(frag_list[i], frag)
        frag_atom = random.choice([atom for atom in frag_list[i].atom_list if atom.can_bond])
        # TODO need to add system to detect if atom is maxed out on bonds (i.e. has 3 or more bonds)

        molecule_atom = random.choice([atom for atom in molecule.atom_list if atom.can_bond])
        frag_atom.bonded_to.append(Bond(molecule_atom, 1))
        molecule_atom.bonded_to.append(Bond(frag_atom, 1))
        for atom in frag_list[i].atom_list:
            molecule.atom_list.append(atom)

    convert_to_smiles(molecule, molecule.atom_list[0])


for i in range(100000):
    build_molecule(MoleculeStructure(), 10)


