from Smiles_to_Structure import convert_to_structure, MoleculeStructure
from fragments import heterocycles
from collections import Counter
from termcolor import cprint

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


def calc_branch_length(branch):  # TODO not sure if this is the right place for this
    branch_length = 0

    def add_daughter_branch_length(daughter_branch):
        nonlocal branch_length
        branch_length += len(daughter_branch.sequence)
        if len(daughter_branch.sequence[-1].daughter_branches) > 0:
            for b in daughter_branch.sequence[-1].daughter_branches:
                add_daughter_branch_length(b)

    add_daughter_branch_length(branch)
    print(branch_length)
    return branch_length


def find_fragment(molecule_string, fragment_string):

    molecule_structure = convert_to_structure(MoleculeStructure(), molecule_string)
    fragment_structure = convert_to_structure(MoleculeStructure(), fragment_string)

    def find_anchor_atom(fragment):
        for ele in ["Si", "P", "S", "I", "Br", "Cl", "F", "B", "O", "N", "C"]:
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

            # if more than 1 unchecked bonds (i.e. a branch point), create new branch for each unchecked bond
            if len(unchecked_bonds) > 1:
                for bond in unchecked_bonds:
                    if not visited[bond.atom]:
                        print("new branch")
                        new_branch(bond.atom, current_atom, current_atom_data, bond.bond_code)

            # if a contiguous section of branch, add bond info
            elif len(unchecked_bonds) == 1:
                if current_branch:
                    if not visited[unchecked_bonds[0].atom]:
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

    def expand_map(atom_map):
        for branch in atom_map.daughter_branches:
            print("branch:")
            print("branch length: ", end="")
            print(len(branch.sequence))
            print("total branch length: ", end="")
            calc_branch_length(branch)
            print("bond to branch:", end="")
            print(branch.bond)
            for atom_info in branch.sequence:
                print(atom_info.symbol)
                if atom_info.bond:
                    print(atom_info.bond)
                if len(atom_info.daughter_branches) > 0:
                    cprint("branch point", "yellow")
                    expand_map(atom_info)

    print("\n")
    cprint("expanded map:", "yellow")
    expand_map(anchored_fragment_map)
    print("\n")

    def is_potential_anchor(atom, fragment_anchor_atom, atom_list):
        # searches through all atoms in molecules in total_molecule to see if they match the fragment base atom
        # atom -> current atom its checking
        # atom_list is list where potential anchor atoms are stored
        # fragment_anchor_atom is the actual atom object from the fragment structure

        fragment_anchor_atom_bonds = Counter([abbr_bond(bond) for bond in fragment_anchor_atom.bonded_to])
        # count bonds from anchor atom

        if atom.symbol == fragment_anchor_atom.symbol:
            # check to see if atom is the same element

            atom_bonds = Counter([abbr_bond(bond) for bond in atom.bonded_to])
            for key in fragment_anchor_atom_bonds:
                if key not in atom_bonds or fragment_anchor_atom_bonds[key] > atom_bonds[key]:
                    # check 1: are there bonds types in fragment base atom that current atom doesn't have
                    # check 2: does current atom have >= the amount of each bond type compared to fragment base atom
                    # i.e. are the bonds in fragment anchor atom a subset of the bonds of current atom
                    return
            atom_list.append(atom)
            # if all checks passed, atom is a potential base atom and is  stored in a list

    potential_anchor_atoms = []
    # keeping track of atoms that match fragment base atom
    for atom in molecule_structure.atom_list:
        is_potential_anchor(atom, fragment_anchor_atom, potential_anchor_atoms)

    for atom in potential_anchor_atoms:
        print("potential anchor: ")
        print(atom.symbol)

    print("\n")

    def check_anchor_atom(potential_anchor_atom, fragment_map, molecule):
        molecule_atoms = set()
        molecule_atoms.add(potential_anchor_atom)
        # list to keep track of which atoms in the molecule constitute a matched fragment

        def check_branch_point(current_molecule_atom, previous_molecule_atom, map_atom_info, branch_atoms):
            # TODO maybe have it so if branch points add their atoms to the branch that spawned them
            # TODO and then only a branch point with no branch can deliver a successful fragment to the actual molecule list
            branch_point_atoms = set()
            cprint("I'm trying a branch point", "yellow")
            map_atom_info.daughter_branches.sort(key=calc_branch_length, reverse=True) # TODO need to calculate total length of branch and daughter branches
            # this makes longer branches go first -> have to search the longest branch first
            # otherwise a shorter branch might be identified in what is actually the long branch
            # i.e. if atom has ethyl and propyl group, you could find the ethyl group where the propyl group is and then be unable to find propyl group
            # also important - need to caculate the total branch length (including length of all its daughter branches)

            cprint("branch point bonds check", "yellow")
            unchecked_bonds = Counter([abbr_bond(bond) for bond in current_molecule_atom.bonded_to if bond.atom != previous_molecule_atom])
            fragment_branch_point_bonds = Counter([branch.bond for branch in map_atom_info.daughter_branches])
            print(unchecked_bonds)
            print(fragment_branch_point_bonds)
            # subset check on branch point, just to make sure current atom has all the bonds the fragment branchpoint has
            for key in fragment_branch_point_bonds:
                if key not in unchecked_bonds or fragment_branch_point_bonds[key] > unchecked_bonds[key]:
                    cprint("branch point doesn't contain necessary bonds", "red")
                    return False

            branch_check = {}
            for branch in map_atom_info.daughter_branches:
                branch_check[branch] = False
            # set all branches to unconfirmed

            trial_paths = [bond for bond in current_molecule_atom.bonded_to if bond.atom != previous_molecule_atom]
            # available routes to see if branch matches
            checked_paths = []
            # keeps track of bonds that have been used to successfully identify branches
            for branch in map_atom_info.daughter_branches:
                # take each branch
                for bond in trial_paths:
                    if branch.bond == abbr_bond(bond) and bond not in checked_paths:
                        # if the bond to the branch matches the current bond (and the current bond hasn't already been used to identify a branch):
                        if try_branch(branch.sequence, bond.atom, current_molecule_atom, branch_point_atoms):
                            # test to see if the current branch works on this bond path
                            cprint("branch successful", "blue")
                            branch_check[branch] = True
                            checked_paths.append(bond)
                            # if true, the branch was successfully found, turn branch to True in branch check
                            # add bond used to checked_paths so it isn't used for further branches
                            break
                            # need to stop the for loop so it doesn't search the matched branch in further trial_paths
                        else:
                            cprint("branch not successful", "red")

            if all(value is True for value in branch_check.values()):
                cprint("branch point match", "blue")
                if branch_atoms:
                    branch_atoms.update(branch_point_atoms)  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@22
                else:
                    molecule_atoms.update(branch_point_atoms)
                    # first branch point does not have a branch that spawned it
                return True
                # if all branches have been found, they will be True in branch_check, branch point is a match, return True
            else:
                cprint("branch point not matched", "red")
                return False
                # one or more branches were found, branch point wasn't a match, return False

        def check_atom_bonds(current_molecule_atom, previous_molecule_atom, branch_sequence, index, branch_atoms):
            cprint("checking atom", "yellow")
            print(current_molecule_atom.symbol)
            branch_atoms.add(current_molecule_atom)
            if len(branch_sequence[index].daughter_branches) > 0:
                # atom is branch point and need to check branches
                return check_branch_point(current_molecule_atom, previous_molecule_atom, branch_sequence[index], branch_atoms)
            else:
                # atom is either an endpoint or contiguous segment:
                if not branch_sequence[index].bond:
                    cprint("reached branch end", "blue")
                    # if no bond data, means we have matched the entire branch, return True
                    return True
                else:
                    # else: this is a contiguous segment look for appropriate bonds
                    unchecked_bonds = [bond for bond in current_molecule_atom.bonded_to if bond.atom != previous_molecule_atom]
                    if len(unchecked_bonds) == 0:
                        cprint("branch ended too early", "red")
                        # actual branch has ended, but map says there should be another atom bonded here, therefore return False
                        return False
                    elif len(unchecked_bonds) == 1:
                        # actual molecule only has a contiguous segment here
                        print(branch_sequence[index].bond)
                        print(abbr_bond(unchecked_bonds[0]))
                        if branch_sequence[index].bond == abbr_bond(unchecked_bonds[0]):
                            return check_atom_bonds(unchecked_bonds[0].atom, current_molecule_atom, branch_sequence, index + 1, branch_atoms)  # check next atom
                            # all branches should either return a function, True, or False.  All child functions should do the same
                            # uncheck_bonds[0].atom becomes new current_molecule_atom, current_molecule_atom becomes previous_molecule_atom
                            # also pass the branch sequence and the index of the next atom_info in the branch
                        else:
                            cprint("wrong bond", "red")
                            return False
                            # the next atom in the branch either has the wrong bond or atom symbol
                    else:
                        # there are multiple possible paths branch could go
                        cprint("checking multiple paths for branch", "yellow")
                        # check all ways
                        for bond in unchecked_bonds:
                            if branch_sequence[index].bond != abbr_bond(bond):  # this is purely for seeing what's happening
                                print(abbr_bond(bond))
                                print(branch_sequence[index].bond)
                            if branch_sequence[index].bond == abbr_bond(bond):
                                print(abbr_bond(bond))
                                print(branch_sequence[index].bond)
                                # looks at all possible ways that match the correct bond
                                midway_fork = set()
                                if check_atom_bonds(bond.atom, current_molecule_atom, branch_sequence, index + 1, midway_fork):
                                    branch_atoms.update(midway_fork)
                                    return True
                                    # return True if any of the paths work (also returns first found)
                        return False
                        # if this else clause didn't return True (i.e. none of the paths succeeded)
                        # then none of the paths are productive, return false

        def try_branch(branch_sequence, current_molecule_atom, previous_molecule_atom, branch_point_atoms):
            branch_atoms = set()
            cprint("I'm trying a branch!", "yellow")
            if check_atom_bonds(current_molecule_atom, previous_molecule_atom, branch_sequence, 0, branch_atoms):
                branch_point_atoms.update(branch_atoms)  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
                # add branch_atoms to branch point_atoms
                return True

        if check_branch_point(potential_anchor_atom, None, fragment_map, None):
            cprint("number of atoms in fragment: ", "yellow", end="")
            print(len(molecule_atoms))
            for atom in molecule_atoms:
                print(atom.symbol)
            cprint("matched fragment to anchor atom", "magenta")
            return True
        else:
            cprint("anchor atom not matched to fragment", "red")
            return False
        # start from check_branch point on the potential anchor atom
        # the anchor atom in map is treated as a branch point, even if it only has 1 branch

    fragment_counter = 0
    for atom in potential_anchor_atoms:
        print("checking anchor atom")
        if check_anchor_atom(atom, anchored_fragment_map, molecule_structure):
            fragment_counter += 1

    print("\n")
    cprint("number of fragments found:", "yellow", end="")
    cprint(fragment_counter, "magenta")


find_fragment("O=C2(C=CC1(=CC=C(C(=C1N2)C3(=NC(=NC=N3)C4(=C(N=C(N=C4)N5(C(=O)CCC5))C6(C(=O)CCCC6))))C%11(=C7(C(C=C(O7)N8(NCCC8C9(=CC=C%10(OC=CCC%10=C9))))=CC=C%11C%12(=NC=C(S%12)N%14(C%13(=CC=CC=C%13C=C%14)))))))", "O=C1C=CC2=CC=C(C(C3=NC=NC=N3)=C2N1)C4=C5C(C=CO5)=CC=C4")