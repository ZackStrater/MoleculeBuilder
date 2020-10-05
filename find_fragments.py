from Smiles_to_Structure import convert_to_structure, MoleculeStructure
from collections import Counter
from termcolor import cprint


class AtomData:
    def __init__(self, symbol):
        self.symbol = symbol
        self.bond = None
        self.daughter_branches = []
        self.ring_closures = set()


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

        atom_info_dict = {}
        # links the molecule_atom and the atom_info representing that atom, used to pass ring_closure info to map

        ring_closure_counter = 1

        def dfs(current_atom, previous_atom, current_branch):
            visited[current_atom] = True

            current_atom_data = AtomData(current_atom.symbol)
            # data object for current atom

            atom_info_dict[current_atom] = current_atom_data

            if current_branch:
                current_branch.sequence.append(current_atom_data)
                # append atom info to branch sequence
                # if current_branch b/c first atom does not have a branch(maybe it should???)

            unchecked_bonds = [bond for bond in current_atom.bonded_to if bond.atom != previous_atom]

            nonlocal ring_closure_counter

            # if more than 1 unchecked bonds (i.e. a branch point), create new branch for each unchecked bond
            if len(unchecked_bonds) > 1:
                for bond in unchecked_bonds:
                    if not visited[bond.atom]:
                        print("new branch")
                        new_branch(bond.atom, current_atom, current_atom_data, bond.bond_code)
                    elif not bool(current_atom_data.ring_closures & atom_info_dict[bond.atom].ring_closures):
                        # if visited[bond.atom], we are at a ring closure
                        # this bool sees if the atom_info of these two atoms (current atom and the atom its bonded to) share any values
                        # if they do, this ring closure has already been documented and we don't want to double count it
                        current_atom_data.ring_closures.add(ring_closure_counter)
                        atom_info_dict[bond.atom].ring_closures.add(ring_closure_counter)
                        # add matching values to each atom_info.ring_closure
                        ring_closure_counter += 1

            # if a contiguous section of branch, add bond info
            elif len(unchecked_bonds) == 1:
                if current_branch:
                    if not visited[unchecked_bonds[0].atom]:
                        print("contin branch")
                        current_atom_data.bond = abbr_bond(unchecked_bonds[0])
                        dfs(unchecked_bonds[0].atom, current_atom, current_branch)
                    elif not bool(current_atom_data.ring_closures & atom_info_dict[unchecked_bonds[0].atom].ring_closures):
                        current_atom_data.ring_closures.add(ring_closure_counter)
                        atom_info_dict[unchecked_bonds[0].atom].ring_closures.add(ring_closure_counter)
                        ring_closure_counter += 1
                        # same as above
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
                if len(atom_info.ring_closures) > 0:
                    cprint("ring closures:", "yellow")
                    for num in atom_info.ring_closures:
                        print(num)
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
            branch_point_atoms = set()
            cprint("I'm trying a branch point", "yellow")
            map_atom_info.daughter_branches.sort(key=calc_branch_length, reverse=True)
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
                    branch_atoms.update(branch_point_atoms)
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
                                # need to separate the branch_atoms here since we don't know if any of the paths will work
                                if check_atom_bonds(bond.atom, current_molecule_atom, branch_sequence, index + 1, midway_fork):
                                    branch_atoms.update(midway_fork)
                                    # if one of the paths works, add all the atoms from the midway_fork "branch"
                                    return True
                                    # return True if any of the paths work (also returns first found)
                        return False
                        # if this else clause didn't return True (i.e. none of the paths succeeded)
                        # then none of the paths are productive, return false

        def try_branch(branch_sequence, current_molecule_atom, previous_molecule_atom, branch_point_atoms):
            branch_atoms = set()
            cprint("I'm trying a branch!", "yellow")
            if check_atom_bonds(current_molecule_atom, previous_molecule_atom, branch_sequence, 0, branch_atoms):
                branch_point_atoms.update(branch_atoms)
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

    # TODO for cycle detection, give atom_info a ring closing attribute, starts as NONe
    # TODO when identify ring closing event, need to give both atom_info a ring closing number (i.e. 1 for ring 1, 2 for ring 2)
    # TODO maybe need to make a dictionary while making map with atoms as keys and atom info as values (otherwise how do you modify previous atom_info???)
    # TODo when traversing anchor atom, and hit ring closing atom -> mark it somehow with the appropriate number
    # TODO then when you get to the ring closing atom of the corresponding number, need to make sure its bonded to the other atom with that number


find_fragment("O=C(OC1=CC=C(C2=C(C3=NC=CN3)NC=C2)C=C41)C=C4C5=C(C6COC(C7=NC8=CC=CC(C9=CC=CO9)=C8O7)CO6)C%10=C(N%11NCCC%11)C(C%12NCCC%12)=CC(C%13=NC=NC%14=NC=CN=C%13%14)=C%10N=C5", "O=C1C=C2C3=C(C4C5C2)C(C6C(C4=C7CC5)=C(C=C7)CCC6CC8)=C8C=C3O1")