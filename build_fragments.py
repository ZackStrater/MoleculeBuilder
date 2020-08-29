from fragments import fragments
import random
from Smiles_to_Structure import MoleculeStructure, convert_to_structure
from molecule_builder import convert_to_smiles

new_molecule = MoleculeStructure()
frag1 = MoleculeStructure
frag2 = MoleculeStructure
frag3 = MoleculeStructure
frag4 = MoleculeStructure
frag_list = [frag1, frag2, frag3, frag4]


def build_molecule(molecule):

    for i in range(4):
        get_frag = random.choice(list(fragments))
        frag = fragments[get_frag]
        convert_to_structure(frag_list[i], frag)
        # choose random atom in frag and random atom in molecule, bond them together, then transfer atoms from frag# to molecule

    convert_to_smiles(molecule, molecule.atom_list[0])


build_molecule(new_molecule)