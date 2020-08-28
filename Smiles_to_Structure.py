from molecule_builder import Atom, Bond, convert_to_smiles
import re


class MoleculeStructure:

    def __init__(self):
        self.atom_list = []


bond_encoder = {
    "=": 2,  # double bond
    "#": 3,  # triple bond
    # aromatic bonds (4) denoted by lowercase element symbol
    "$": 5,  # ionic bond
    "/": 6,  # vinyl up
    "\\": 7,  # vinyl down
    ".": 8,  # non-bond
}


def encode_bond(bonding_info):
    code = None
    for key in bond_encoder:
        if key in bonding_info:
            code = (bond_encoder[key])
    if code is None:
        return 1
    else:
        return code


new_molecule = MoleculeStructure()

smiles_string_input = "C1%11(=C(C(=C(C(=C1)C2(=CC=C(C=C2)CCCOC4(CCC3(=CC=C5(C=C7(C(=CC5=C34)C6(CC67)))))N))C(C)CC(C)C)C8(=C9(C(=CC(=C8)C=O)CCCN9)))C%10(CC%10%11))"
atom_map = []
for match in re.finditer(r"[A-Z][a-z]?", smiles_string_input):
    atom_map.append((match.start(), Atom(match.group())))
    # ordered list of tuples for atoms -> (index, Atom)
    # atom object created in this step

bond_map = re.findall(r"(?:[A-Z][a-z]?)([^A-Za-z]*)", smiles_string_input)
# ordered list containing all bonding information between atoms

for ele in atom_map:
    new_molecule.atom_list.append(ele[1])
    # adding all atoms to new_molecule

left_parens_list = []
# ordered list of atoms that are followed by a "("
ring_closure_dict = {}
# keeps track of ring closures


for i in range(0, len(new_molecule.atom_list)):
    # for loop through all Atoms, but allows access to Atom (i + 1)
    all_bond_info = (bond_map[i])
    # string with all bonding info for current atom
    ring_closure_info = re.findall(r"(\D*[%]\d{2}|\D*\d)", all_bond_info)
    # pulls out all ring #'s and the bond info that proceeds them

    for ele in ring_closure_info:
        ring_closure_split = re.split(r"(\D*(?=\d))", ele, 1)
        # separates bonding [1] info from ring # [2]
        if ring_closure_split[2] in ring_closure_dict:
            # if this ring # has already been encountered -> make bond between those two Atoms
            closure_partner = ring_closure_dict[ring_closure_split[2]]
            combined_bonding_info = closure_partner[1] + ring_closure_split[1]
            # often bonding info only included for one of the ring closure partners, so combine and use for both Bonds
            new_molecule.atom_list[i].bonded_to.append(Bond(closure_partner[0], encode_bond(combined_bonding_info)))
            closure_partner[0].bonded_to.append(Bond(new_molecule.atom_list[i], encode_bond(combined_bonding_info)))
            # adding bonds to both closure partners
        else:
            ring_closure_dict[ring_closure_split[2]] = (new_molecule.atom_list[i], ring_closure_split[1])
            # enter atom in ring_closure_dict -> ring #: (Atom, bonding info)

    bond_info = re.split(r"(\D*$)", all_bond_info, 1)[1]
    # captures the rest of bonding info after ring closure #'s
    if ")" not in bond_info:
        # make Bond between current Atom and Atom (i + 1)
        if i != (len(new_molecule.atom_list) - 1):
            # prevents index error on last Atom
            new_molecule.atom_list[i].bonded_to.append(Bond(new_molecule.atom_list[i + 1], encode_bond(bond_info)))
            new_molecule.atom_list[i + 1].bonded_to.append(Bond(new_molecule.atom_list[i], encode_bond(bond_info)))
            if "(" in bond_info:
                left_parens_list.append(new_molecule.atom_list[i])

    else:
        # if ")" in bond_info, count the # of ")"'s and go back that many "("'s in the string
        # and make Bond to preceding Atom
        if i != (len(new_molecule.atom_list) - 1):
            parentheses_count = bond_info.count(")")
            # counts # of ")" in bond_info and bonds current Atom to
            bond_atom1 = left_parens_list[-parentheses_count]
            # atom that opened the original bracket
            bond_atom2 = new_molecule.atom_list[i + 1]
            bond_atom2.bonded_to.append(Bond(bond_atom1, encode_bond(bond_info)))
            bond_atom1.bonded_to.append(Bond(bond_atom2, encode_bond(bond_info)))
            if "(" in bond_info:  # bond info has "C)(C"
                for c in range(parentheses_count - 1):
                    left_parens_list.pop()
                    # this keeps track of last add "C(" instance
                    # i.e. A(B(C))(DE)F , for bond info between C and D -> C))(D
                    # we get two ")" so we should delete two atoms off left_parens_list
                    # but we also have "(" so we need to add an atom back on the list
                    # the correct atom is not the previous one, but actually the second one we deleted from left_parens_list
                    # so instead we just delete one less from the list and now it will correctly bond A to F

            else:
                for c in range(parentheses_count):
                    left_parens_list.pop()


convert_to_smiles(new_molecule, new_molecule.atom_list[0])


