from fragments import fragments
import random
from Smiles_to_Structure import MoleculeStructure, convert_to_structure
from molecule_builder import convert_to_smiles, Bond


frag1 = MoleculeStructure()
frag2 = MoleculeStructure()
frag3 = MoleculeStructure()
frag4 = MoleculeStructure()
frag_list = [frag1, frag2, frag3, frag4]


def build_molecule(molecule):

    for i in range(4):
        if i == 0:
            get_frag = random.choice(list(fragments))
            frag = fragments[get_frag]
            convert_to_structure(frag1, frag)  # maybe make this an object somehow???
            for atom in frag1.atom_list:
                molecule.atom_list.append(atom)
        if i > 0:
            get_frag = random.choice(list(fragments))
            frag = fragments[get_frag]
            convert_to_structure(frag_list[i], frag)
            frag_atom = random.choice(frag_list[i].atom_list)
            molecule_atom = random.choice(molecule.atom_list)
            frag_atom.bonded_to.append(Bond(molecule_atom, 1))
            molecule_atom.bonded_to.append(Bond(frag_atom, 1))
            for atom in frag_list[i].atom_list:
                molecule.atom_list.append(atom)

    convert_to_smiles(molecule, molecule.atom_list[0])
    frag1.atom_list.clear()
    frag2.atom_list.clear()
    frag3.atom_list.clear()
    frag4.atom_list.clear()


for i in range(100000):
    build_molecule(MoleculeStructure())


