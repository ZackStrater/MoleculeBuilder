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

            unchecked_bonds = [bond for bond in current_atom.bonded_to if bond.atom != previous_atom]
            num_unchecked_bonds = len(unchecked_bonds)

            # if more than 1 unchecked bonds (i.e. a branch point), create new branch for each unchecked bond
            if num_unchecked_bonds > 1:
                print("new branches")
                for bond in unchecked_bonds:
                    if not visited[bond.atom]:
                        new_branch(bond.atom, current_atom, current_atom_data, bond.bond_code)

            # if a contiguous section of branch, add bond info
            elif num_unchecked_bonds == 1:
                if current_branch:
                    print("contin branch")
                    current_atom_data.bond = abbr_bond(unchecked_bonds[0])
                    dfs(unchecked_bonds[0].atom, current_atom, current_branch)
                else:
                    print("new branch")
                    for bond in unchecked_bonds:
                        new_branch(bond.atom, current_atom, current_atom_data, bond.bond_code)
                    # if the anchor atom only has 1 bond, need to start a branch

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
        print("branch length: ", end="")
        print(len(branch.sequence))
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
        print("potential anchor: ")
        print(atom.symbol)

    print("\n")

    def check_anchor_atom(potential_anchor_atom, fragment_map, molecule):

        def check_branch_point(map_atom_info, current_molecule_atom, previous_molecule_atom):
            print("I'm trying a branch point")
            map_atom_info.daughter_branches.sort(key=lambda x: len(x.sequence), reverse=True)
            # this makes longer branches go first -> have to search the longest branch first
            # otherwise a shorter branch might be identified in what is actually the long branch
            # i.e. if atom has ethyl and propyl group, you could find the ethyl group where the propyl group is
            for branch in map_atom_info.daughter_branches:
                for bond in current_molecule_atom.bonded_to:
                    if branch.bond == abbr_bond(bond):
                        try_branch(branch.sequence, 0, bond.atom, current_molecule_atom)

        def check_atom(current_molecule_atom, previous_molecule_atom, branch_sequence, index):  # TODO need a way to pass atom info to branches
            if len(branch_sequence[index].daughter_branches) > 0:
                # atom is branch point and need to check branches
                check_branches(current_molecule_atom, previous_molecule_atom, branch_sequence[index])
            else:
                # atom is either an endpoint or contiguous segment:
                if not branch_sequence[index].bond:
                    # if no bond data, means we have matched the entire branch, return True
                    return True
                else:
                    # else: this is a contiguous segment look for appropriate bonds
                    unchecked_bonds = [bond for bond in current_molecule_atom.bonded_to if bond.atom != previous_molecule_atom]
                    if len(unchecked_bonds) == 0:
                        # actual branch has ended, but map says there should be another atom bonded here, therefore return False
                        return False
                    elif len(unchecked_bonds) == 1:
                        # actual molecule only has a contiguous segment here
                        if branch_sequence[index].bond == abbr_bond(unchecked_bonds[0]):
                            check_atom(unchecked_bonds[0].atom, current_molecule_atom, branch_sequence, index + 1)  # check next atom
                            # uncheck_bonds[0].atom becomes new current_molecule_atom, current_molecule_atom becomes previous_molecule_atom
                            # also pass the branch sequence and the index of the next atom_info in the branch
                    else:
                        # check all ways
                        for bond in unchecked_bonds:
                            # check bond for each possible direction

        def check_branches(current_molecule_atom, previous_molecule_atom, atom_data):
            print("hello")

        def try_branch(branch_sequence, index, molecule_atom, previous_molecule_atom):
            print("I'm trying a branch!")
            branch_atoms = []
            unchecked_bonds = [bond for bond in molecule_atom if bond.atom != previous_molecule_atom]
            if len(unchecked_bonds) == 1:
                for bond in unchecked_bonds:
                    if abbr_bond(bond) == branch_sequence[index]:
                        print("hello")








    for atom in potential_anchor_atoms:
        check_anchor_atom(atom, anchored_fragment_map, molecule_structure)
find_fragment("BrCCC(CCCl)(CCS)CCSiC", "BrCCCCCCl")