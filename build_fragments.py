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
        for atom in frag_list[i].atom_list:
            molecule.atom_list.append(atom)

    convert_to_smiles(molecule, molecule.atom_list[0])


for i in range(100000):
    build_molecule(MoleculeStructure(), 10)


