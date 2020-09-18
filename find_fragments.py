from Smiles_to_Structure import convert_to_structure, MoleculeStructure
from fragments import heterocycles
from collections import Counter


class AtomData:
    def __init__(self, symbol):
        self.symbol = symbol
        self.bond = None
        self.daughter_branches = []


class Branch:
    def __init__(self, bond):
        self.bond = bond
        self.sequence = []


def abbr_bond(bond):
    return bond.bond_code, bond.atom.symbol


def find_fragment(molecule_string, fragment_string):  #specify fragment list

    molecule_structure = convert_to_structure(MoleculeStructure(), molecule_string)
    fragment_structure = convert_to_structure(MoleculeStructure(), fragment_string)

    def find_anchor_atom(fragment):
        for ele in ["Si", "P", "I", "Br", "Cl", "F", "B", "S", "O", "N", "C"]:
            for atom in fragment.atom_list:
                if atom.symbol == ele:
                    return atom

    fragment_anchor_atom = find_anchor_atom(fragment_structure)
    # the actual atom object of highest priority in the fragment structure

    def map_fragment(fragment, anchor_atom):

        visited = {}
        for atom in fragment.atom_list:
            visited[atom] = False
        # keeps track of which atoms have been visited

        def dfs(current_atom, previous_atom, current_branch):
            visited[current_atom] = True

            current_atom_data = AtomData(current_atom.symbol)
            # data object for current atom

            if current_branch:
                current_branch.sequence.append(current_atom_data)
                # append atom info to branch sequence
                # if current_branch b/c first atom does not have a branch(maybe it should???)

            previous_bond = None
            for bond in current_atom.bonded_to:
                if bond.atom == previous_atom:
                    previous_bond = bond

            unchecked_bonds = [bond for bond in current_atom.bonded_to if bond != previous_bond]
            num_unchecked_bonds = len(unchecked_bonds)

            # if more than 1 unchecked bonds (i.e. a branch point), create new branch for each unchecked bond
            if num_unchecked_bonds > 1:
                print("new branches")
                for bond in unchecked_bonds:
                    if not visited[bond.atom]:
                        new_branch(bond.atom, current_atom, current_atom_data, bond.bond_code)

            # if a contiguous section of branch, add bond info
            elif num_unchecked_bonds == 1:
                print("contin branch")
                for bond in unchecked_bonds:
                    current_atom_data.bond = abbr_bond(bond)
                for bond in unchecked_bonds:
                    if not visited[bond.atom]:
                        dfs(bond.atom, current_atom, current_branch)

            else:
                print("end point")

            if not current_branch:
                return current_atom_data
                # this returns anchor atom to the map_fragment function

        def new_branch(current_atom, previous_atom, previous_atom_data, bond_code):
            current_branch = Branch((bond_code, current_atom.symbol))
            # create new branch with bonding info to first atom in branch
            previous_atom_data.daughter_branches.append(current_branch)
            # add new branch to the atom which spawned it
            dfs(current_atom, previous_atom, current_branch)
            # start dfs on first atom in branch
            # need to pass previous_atom in order to not travel backwards

        # TODO add ring detection eventually

        return dfs(anchor_atom, None, None)
        # starts process of mapping fragment, but also returns the anchor atom

    anchored_fragment_map = map_fragment(fragment_structure, fragment_anchor_atom)
    # the map base is the atom_data representation of the anchor atom
    # the rest of the map is stored in the daughter branches

    print(anchored_fragment_map.symbol)
    for branch in anchored_fragment_map.daughter_branches:
        print("branch:")
        print(branch.bond)
        for atom in branch.sequence:
            print(atom.symbol)
            print(atom.bond)

    def is_anchor(atom, fragment_anchor_atom, atom_list):
        # searches through all atoms in molecules in total_molecule to see if they match the fragment base atom
        # atom -> current atom its checking

        fragment_anchor_atom_bonds = Counter([abbr_bond(bond) for bond in fragment_anchor_atom.bonded_to])
        # count bonds from anchor atom

        if atom.symbol == fragment_anchor_atom.symbol:
            # check to see if atom is the same element

            atom_bond_list = Counter([abbr_bond(bond) for bond in atom.bonded_to])
            for key in fragment_anchor_atom_bonds:
                if key not in atom_bond_list or fragment_anchor_atom_bonds[key] > atom_bond_list[key]:
                    # check 1: are there bonds types in fragment base atom that current atom doesn't have
                    # check 2: does current atom have >= the amount of each bond type compared to fragment base atom
                    # i.e. are the bonds in fragment anchor atom a subset of the bonds of current atom
                    return
            atom_list.append(atom)
            # if all checks passed, atom is a potential base atom and is  stored in a list

    potential_anchor_atoms = []
    # keeping track of atoms that match fragment base atom
    for atom in molecule_structure.atom_list:
        is_anchor(atom, fragment_anchor_atom, potential_anchor_atoms)

    for atom in potential_anchor_atoms:
        print(atom.symbol)


    def check_anchor_atom(potential_anchor_atom, fragment_map, molecule):


    for atom in potential_anchor_atoms:
        check_anchor_atom(atom, anchored_fragment_map, molecule_structure)
find_fragment("BrCCSi(CCCl)(CCS)CCSiC", "BrCCSiCCCl")