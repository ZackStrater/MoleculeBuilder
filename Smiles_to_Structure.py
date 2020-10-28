from Strucutre_to_Smiles import Atom, Bond, convert_to_smiles
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


smiles_string_input = "c1ccccc1"

# TODO remove molecule from this func, make it just part of the func
# TODO add removal of H, [], and @ from string (for now)
def convert_to_structure(molecule, smiles_string):

    corrected_smiles_string = re.sub(r"[\[\]H@]", "", smiles_string)  # TODO need to change this at some point

    for match in re.findall(r"N|O|P|Si|S|F|Cl|Br|I|C|B|\[R\]|b|c|n|o|p|s|W|X|Y|Z", corrected_smiles_string):
        molecule.atom_list.append(Atom(match))
        # create Atom object for each elemental symbol found in the smiles_string, keeps order of Atom in the string

    bond_map = re.findall(r"(?:N|O|P|Si|S|F|Cl|Br|I|C|B|R|b|c|n|o|p|s|W|X|Y|Z)([^A-Za-z]*)", corrected_smiles_string)
    # ordered list containing all bonding symbol denoting bonding information following each element symbol
    # TODO right now this regex finds R instead of [R], either should remove "[" and "]" from fragment library
    # TODO at some point might need to add comprehension for charged parts i.e [nH4+]
    # TODO figure out what to do about H, (i.e. C[@@H] etc.... 

    left_parens_list = []
    # ordered list of atoms that are followed by a "("
    ring_closure_dict = {}
    # keeps track of ring closures

    for i in range(0, len(molecule.atom_list)):
        # for loop through all Atoms, but allows access to Atom (i + 1)
        # notice this range includes the last item in the list cause ring closure can be included in last atom
        all_bond_info = (bond_map[i])
        # string with all bonding info for current atom
        ring_closure_info = re.findall(r"(\D*[%]\d{2}|\D*\d)", all_bond_info)
        # pulls out all ring #'s and the bond info that proceeds them

        for ele in ring_closure_info:
            ring_closure_split = re.split(r"(\D*(?=\d))", ele, 1)
            # separates bonding info [1] from ring # [2]
            if ring_closure_split[2] in ring_closure_dict:
                # if this ring # has already been encountered -> make bond between those two Atoms
                closure_partner = ring_closure_dict[ring_closure_split[2]]
                combined_bonding_info = closure_partner[1] + ring_closure_split[1]
                # often bonding info only included for one of the ring closure partners, so combine and use for both Bonds
                molecule.atom_list[i].bonded_to.append(Bond(closure_partner[0], encode_bond(combined_bonding_info)))
                closure_partner[0].bonded_to.append(Bond(molecule.atom_list[i], encode_bond(combined_bonding_info)))
                # adding bonds to both closure partners
            else:
                ring_closure_dict[ring_closure_split[2]] = (molecule.atom_list[i], ring_closure_split[1])
                # enter atom in ring_closure_dict -> ring #: (Atom, bonding info)

        bond_info = re.split(r"(\D*$)", all_bond_info, 1)[1]
        # captures the rest of bonding info after ring closure #'s
        if ")" not in bond_info:
            # make Bond between current Atom and Atom (i + 1)
            if i != (len(molecule.atom_list) - 1):
                # prevents index error on last Atom TODO see if there is a better way to do this line
                # # TODO maybe split into range 0 - len(molecule_list) for above stuff and a seperate 0-range_moleculelist -1 for this part
                molecule.atom_list[i].bonded_to.append(Bond(molecule.atom_list[i + 1], encode_bond(bond_info)))
                molecule.atom_list[i + 1].bonded_to.append(Bond(molecule.atom_list[i], encode_bond(bond_info)))
                if "(" in bond_info:
                    left_parens_list.append(molecule.atom_list[i])

        else:
            # if ")" in bond_info, count the # of ")"'s and go back that many "("'s in the string
            # and make Bond from atom before last "(" and the next atom (i + 1)
            if i != (len(molecule.atom_list) - 1):  # TODO see if there is a better way to do this line
                # important not to include last item here, ")" at the end of a string don't actually give structural information
                parentheses_count = bond_info.count(")")
                # counts # of ")" in bond_info and bonds current Atom to
                bond_atom1 = left_parens_list[-parentheses_count]
                # atom that opened the original bracket
                bond_atom2 = molecule.atom_list[i + 1]
                # atom that comes after all the ")"'s
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

    for atom in molecule.atom_list:
        if atom.symbol.islower():
            for bond in atom.bonded_to:
                if bond.atom.symbol.islower():
                    bond.bond_code = 4
    # smiles strings sometimes use lowercase to denote an aromatic ring
    # without explicit bond notation, encode_bond will return a 1 (single bond) for bond_code
    # this loop records the aromatic bonding information for strings with such notation,
    # only counting bonds between two lower case atom symbols

        for bond in atom.bonded_to:  # TODO add [Z] element to refer to priority bonding
            if bond.atom.symbol == "[R]":
                atom.can_bond = True
        # here the [R] atoms represent sites that the fragment can bond
        # first atom is bonded to element [R], then can_bond for atom changed to True, then delete [R] atoms from
        # molecule.atom_list and any bonds to [R] atoms

        atom.bonded_to = [bond for bond in atom.bonded_to if bond.atom.symbol != "[R]"]
    molecule.atom_list = [atom for atom in molecule.atom_list if atom.symbol != "[R]"]
    # removing [R] Atom objects and the Bonds to them

    return molecule




