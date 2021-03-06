from Smiles_to_Structure import convert_to_structure, MoleculeStructure
from collections import Counter
from termcolor import cprint
from fragments_library import biomolecules, peptide_amino_acids, heterocycles, arenes, functional_groups, hydrocarbons


class AtomData:
    def __init__(self, symbol):
        self.symbol = symbol
        self.bond = None
        self.daughter_branches = []
        self.ring_closures = set()
        self.phantom_bonds = None
        self.phantom_atom = False


class Branch:
    def __init__(self, bond):
        self.bond = bond
        self.sequence = []


def abbr_bond(bond):
    return bond.bond_code, bond.atom.symbol


def check_bond(bond, map_bond, bond_atom_info):  # bond_atom_info is the atom_info of the atom for the bond being checked (to see if it's a phantom bond/discovered)

    # if the bond.atom is discovered already, it should give back false (can't retread over discovered atoms)
    # unless the atom_info for that atom shows that the bond.atom should be a phantom atom (which can be searched for in discovered atoms)
    if bond.atom.discovered:
        if not bond_atom_info.phantom_atom:
            return False

    if abbr_bond(bond) == map_bond:
        return True
    # need to cover (correct, R) (correct, Q) (9, R) (9, Q) (9, correct)
    elif bond.bond_code == map_bond[0] or map_bond[0] == 9: # covers (correct, R) (correct, Q) (9, R) (9, Q)
        if map_bond[1] == "R":
            return True
        elif map_bond[1] == "Q" and bond.atom.heteroatom:
            return True
    elif bond.atom.symbol == map_bond[1]:
        if map_bond[0] == 9:
            return True
    else:
        return False


def calc_branch_length(branch):
    branch_length = 0

    def add_daughter_branch_length(daughter_branch):
        nonlocal branch_length
        branch_length += len(daughter_branch.sequence)
        if len(daughter_branch.sequence[-1].daughter_branches) > 0:
            for b in daughter_branch.sequence[-1].daughter_branches:
                add_daughter_branch_length(b)

    add_daughter_branch_length(branch)
    return branch_length


def find_fragment(fragment_string, molecule_string, structure=None):

    if structure:
        molecule_structure = structure
    else:
        molecule_structure = convert_to_structure(MoleculeStructure(), molecule_string)
    fragment_structure = convert_to_structure(MoleculeStructure(), fragment_string)

    def find_anchor_atom(fragment):
        for ele in ["Si", "P", "p", "S", "s", "I", "Br", "Cl", "F", "B", "b", "O", "o", "N", "n", "C", "c"]:
            for atom in fragment.atom_list:
                if atom.symbol == ele:
                    return atom

    fragment_anchor_atom = find_anchor_atom(fragment_structure)
    # the actual atom object of highest priority in the fragment structure

    def is_potential_anchor(atom, fragment_anchor_atom, atom_list):
        # searches through all atoms in molecules in total_molecule to see if they match the fragment base atom
        # atom -> current atom its checking
        # atom_list is list where potential anchor atoms are stored
        # fragment_anchor_atom is the actual atom object from the fragment structure

        if atom.discovered and not fragment_anchor_atom.phantom_atom:  # if fragment_anchor atom is a phantom atom, it can use discovered atoms as potential anchors
            return
        # atom has already been used to find a fragment

        if atom.symbol != fragment_anchor_atom.symbol:
            return
            # check to see if atom is the same element

        fragment_anchor_atom_bonds = Counter([abbr_bond(bond) for bond in fragment_anchor_atom.bonded_to])
        # count bonds from anchor atom

        atom_bonds = Counter([abbr_bond(bond) for bond in atom.bonded_to])

        # count bonds in potential anchor atom where the bond's atom haven't been discovered yet (as we won't be able to use those bonds)
        for key in fragment_anchor_atom_bonds:
            if key[1] != "R" and key[1] != "Q" and key[0] != 9:  # TODO better way to do this???
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

    if potential_anchor_atoms == []:
        #cprint("no anchor atoms found", "yellow")
        return 0
    #else:
    #    print("potential anchor atoms: ")
    #    for atom in potential_anchor_atoms:
    #        print(atom.symbol)
    #        for bond in atom.bonded_to:
    #            print(abbr_bond(bond))

    def map_fragment(fragment, anchor_atom):

        visited = {}
        for atom in fragment.atom_list:
            visited[atom] = False
        # keeps track of which atoms have been visited

        atom_info_dict = {}
        # links the molecule_atom and the atom_info representing that atom, used to pass ring_closure info to map

        ring_closure_counter = 1

        def traverse(current_atom, previous_atom, current_branch):
            visited[current_atom] = True

            current_atom_data = AtomData(current_atom.symbol)
            # data object for current atom

            # atom_data will reflect that the atom is a phantom_atom
            if current_atom.phantom_atom:
                current_atom_data.phantom_atom = True

            atom_info_dict[current_atom] = current_atom_data

            if current_branch:
                current_branch.sequence.append(current_atom_data)
                # append atom info to branch sequence
                # if current_branch b/c first atom does not have a branch

            current_atom_data.phantom_bonds = current_atom.phantom_bonds

            unchecked_bonds = [bond for bond in current_atom.bonded_to if bond.atom != previous_atom]

            nonlocal ring_closure_counter

            # if more than 1 unchecked bonds (i.e. a branch point), create new branch for each unchecked bond
            if len(unchecked_bonds) > 1:
                for bond in unchecked_bonds:
                    if not visited[bond.atom]:
                        #print("new branch")
                        new_branch(bond.atom, current_atom, current_atom_data, bond.bond_code)
                    elif not bool(current_atom_data.ring_closures & atom_info_dict[bond.atom].ring_closures):
                        # if visited[bond.atom], we are at a ring closure
                        # this bool sees if the atom_info of these two atoms (current atom and the atom its bonded to) share any values (& operator)
                        # if they do, this ring closure has already been documented and we don't want to double count it
                        #print("ring closure")
                        current_atom_data.ring_closures.add((ring_closure_counter, bond.bond_code))
                        atom_info_dict[bond.atom].ring_closures.add((ring_closure_counter, bond.bond_code))
                        # add matching values to each atom_info.ring_closure
                        # ring closure data in format (ring closure #, bond_code)
                        ring_closure_counter += 1

            # if a contiguous section of branch, add bond info
            elif len(unchecked_bonds) == 1:
                if current_branch:
                    if not visited[unchecked_bonds[0].atom]:
                        #print("contin branch")
                        current_atom_data.bond = abbr_bond(unchecked_bonds[0])
                        traverse(unchecked_bonds[0].atom, current_atom, current_branch)
                    elif not bool(current_atom_data.ring_closures & atom_info_dict[unchecked_bonds[0].atom].ring_closures):
                        #print("ring closure")
                        current_atom_data.ring_closures.add((ring_closure_counter, unchecked_bonds[0].bond_code))
                        atom_info_dict[unchecked_bonds[0].atom].ring_closures.add((ring_closure_counter, unchecked_bonds[0].bond_code))
                        ring_closure_counter += 1
                        # same as above
                else:
                    #print("new branch")
                    for bond in unchecked_bonds:
                        new_branch(bond.atom, current_atom, current_atom_data, bond.bond_code)
                    # if the anchor atom only has 1 bond, need to start a branch

            #else:
                #print("end point")

            if not current_branch:
                return current_atom_data
                # this returns anchor atom to the map_fragment function as the anchor atom is not spawned from a branch @

        def new_branch(current_atom, previous_atom, previous_atom_data, bond_code):
            current_branch = Branch((bond_code, current_atom.symbol))
            # create new branch with bonding info to first atom in branch
            previous_atom_data.daughter_branches.append(current_branch)
            # add new branch to the atom which spawned it
            traverse(current_atom, previous_atom, current_branch)
            # start traverse on first atom in branch
            # need to pass previous_atom in order to not travel backwards

        return traverse(anchor_atom, None, None)
        # starts process of mapping fragment, but also returns the anchor atom

    anchored_fragment_map = map_fragment(fragment_structure, fragment_anchor_atom)
    # the map base is the atom_data representation of the anchor atom
    # the rest of the map is stored in the daughter branches

    def expand_map(anchor_atom):
        print("anchor atom")
        print(anchor_atom.symbol)
        if len(anchor_atom.ring_closures) > 0:
            cprint("ring closures:", "yellow")
            for num in anchor_atom.ring_closures:
                cprint(num, "blue")
        if anchor_atom.phantom_bonds:
            print("phantom bonds = ", end="")
            print(anchor_atom.phantom_bonds)

        def expand_branch_point(atom_map):
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
                            cprint(num, "blue")
                    if atom_info.phantom_bonds:
                        print("phantom bonds = ", end="")
                        print(atom_info.phantom_bonds)
                    if atom_info.bond:
                        print(atom_info.bond)
                    if len(atom_info.daughter_branches) > 0:
                        cprint("branch point", "yellow")
                        expand_branch_point(atom_info)
        expand_branch_point(anchor_atom)

    #print("\n")
    #cprint("expanded map:", "yellow")
    #expand_map(anchored_fragment_map)
    #print("\n")

    def check_anchor_atom(potential_anchor_atom, fragment_map):
        molecule_atoms = {potential_anchor_atom}  # TODO is this not needed? just use currently visited???
        # list to keep track of which atoms in the molecule constitute a matched fragment

        currently_visited = {potential_anchor_atom: fragment_map}
        # dictionary that keeps track of which atoms have been used to find the fragment at any given step

        def check_branch_point(current_molecule_atom, previous_molecule_atom, map_atom_info, branch_atoms):
            if map_atom_info.phantom_bonds:
                bond_num = len(current_molecule_atom.bonded_to)
                if bond_num != map_atom_info.phantom_bonds:
                    #cprint("wrong amount of phantom bonds", "red")
                    return False
            # phantom_bonds is a way to ensure the current atom is bonded to the specified number of atoms
            # note that phantom bonds includes any bonds for the current molecule_atom, including those to atoms that are "discovered"

            branch_point_atoms = set()
            nonlocal currently_visited
            #cprint("I'm trying a branch point", "yellow")
            map_atom_info.daughter_branches.sort(key=calc_branch_length, reverse=True)
            # this makes longer branches go first -> have to search the longest branch first
            # otherwise a shorter branch might be identified in what is actually the long branch
            # i.e. if atom has ethyl and propyl group, you could find the ethyl group where the propyl group is and then be unable to find propyl group
            # also important - need to caculate the total branch length (including length of all its daughter branches)

            #cprint("branch point bonds check", "yellow")
            unchecked_bonds = Counter([abbr_bond(bond) for bond in current_molecule_atom.bonded_to if bond.atom != previous_molecule_atom])
            fragment_branch_point_bonds = Counter([branch.bond for branch in map_atom_info.daughter_branches])
            #print(unchecked_bonds)
            #print(fragment_branch_point_bonds)
            # subset check on branch point, just to make sure current atom has all the bonds the fragment branchpoint has
            for key in fragment_branch_point_bonds:
                if key[1] != "R" and key[1] != "Q" and key[0] != 9:  # TODO better way to do this?
                    if key not in unchecked_bonds or fragment_branch_point_bonds[key] > unchecked_bonds[key]:
                        #cprint("branch point doesn't contain necessary bonds", "red")
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
                    if check_bond(bond, branch.bond, branch.sequence[0]) and bond not in checked_paths and bond.atom not in currently_visited:
                        # if the bond to the branch matches the current bond (and the current bond hasn't already been used to identify a branch):
                        if try_branch(branch.sequence, bond.atom, current_molecule_atom, branch_point_atoms):
                            # test to see if the current branch works on this bond path
                            #cprint("branch successful", "blue")
                            branch_check[branch] = True
                            checked_paths.append(bond)
                            # if true, the branch was successfully found, turn branch to True in branch check
                            # add bond used to checked_paths so it isn't used for further branches
                            break
                            # need to stop the for loop so it doesn't search the matched branch in further trial_paths
                        #else:
                            #cprint("branch not successful", "red")

            if all(value is True for value in branch_check.values()):
                #cprint("branch point match", "blue")
                if branch_atoms:
                    branch_atoms.update(branch_point_atoms)
                else:
                    molecule_atoms.update(branch_point_atoms)
                    # first branch point does not have a branch that spawned it @
                return True
                # if all branches have been found, they will be True in branch_check, branch point is a match, return True
            else:
                #cprint("branch point not matched", "red")
                return False
                # one or more branches were found, branch point wasn't a match, return False

        def try_branch(branch_sequence, current_molecule_atom, previous_molecule_atom, branch_point_atoms):
            branch_atoms = set()
            #cprint("I'm trying a branch!", "yellow")
            if check_atom_bonds(current_molecule_atom, previous_molecule_atom, branch_sequence, 0, branch_atoms):
                branch_point_atoms.update(branch_atoms)
                # add branch_atoms to branch point_atoms
                return True
            else:
                nonlocal currently_visited
                for a in branch_atoms:
                    currently_visited.pop(a)

        def check_ring_closure(current_atom, atom_info):  # already check if atom_info has ring closure
            ring_closures = set()  # all the ring closure numbers in currently_visited
            for value in currently_visited.values():
                ring_closures.update(value.ring_closures)

            for closure_info in atom_info.ring_closures:  # checking each ring closure atom has
                #print("ring closure")
                #print(closure_info)
                if closure_info in ring_closures:  # we've already hit the other half of the ring closure
                    for key in currently_visited:  # looking for matching ring closure
                        if closure_info in currently_visited[key].ring_closures:  # matched ring closure, key = atom it should be bonded to
                            ring_closure_partner = key

                            if ring_closure_partner in [bond.atom for bond in current_atom.bonded_to]:
                                ring_closure_bond = None
                                for bond in current_atom.bonded_to:
                                    if bond.atom == ring_closure_partner:
                                        ring_closure_bond = bond
                                if ring_closure_bond.bond_code != closure_info[1]:
                                    #cprint("closure bond incorrect", "red")
                                    return False
                            else:
                                #cprint("atom not bonded to correct closure partner", "red")
                                return False
                else:
                    return True
                # first time encountering that ring closure number, don't need to do any further checks
            #cprint("all ring closures acounted for", "blue")
            return True

        def check_atom_bonds(current_molecule_atom, previous_molecule_atom, branch_sequence, index, branch_atoms):
            #cprint("checking atom", "yellow")
            #print(current_molecule_atom.symbol)
            nonlocal currently_visited
            current_atom_info = branch_sequence[index]

            if current_atom_info.phantom_bonds:
                bond_num = len(current_molecule_atom.bonded_to)
                if bond_num != current_atom_info.phantom_bonds:
                    #cprint("wrong amount of phantom bonds", "red")
                    return False

            if current_atom_info.ring_closures:
                if not check_ring_closure(current_molecule_atom, current_atom_info):
                    return False

            currently_visited[current_molecule_atom] = current_atom_info
            #for a in currently_visited:
                #print(a.symbol, end="")
                #if currently_visited[a].ring_closures:
                    #print(currently_visited[a].ring_closures, end="")
            #print("\n")

            branch_atoms.add(current_molecule_atom)
            if len(current_atom_info.daughter_branches) > 0:
                # atom is branch point and need to check branches
                return check_branch_point(current_molecule_atom, previous_molecule_atom, current_atom_info, branch_atoms)
            else:
                # atom is either an endpoint or contiguous segment:
                if not current_atom_info.bond:
                    #cprint("reached branch end", "blue")

                    # if no bond data, means we have matched the entire branch, return True
                    return True
                else:
                    # else: this is a contiguous segment look for appropriate bonds
                    unchecked_bonds = [bond for bond in current_molecule_atom.bonded_to if bond.atom != previous_molecule_atom]
                    if len(unchecked_bonds) == 0:
                        #cprint("branch ended too early", "red")
                        # actual branch has ended, but map says there should be another atom bonded here, therefore return False
                        return False
                    elif len(unchecked_bonds) == 1:
                        # actual molecule only has a contiguous segment here
                        #print(current_atom_info.bond)
                        #print(abbr_bond(unchecked_bonds[0]))
                        if check_bond(unchecked_bonds[0], current_atom_info.bond, branch_sequence[index + 1]) and unchecked_bonds[0].atom not in currently_visited:
                            return check_atom_bonds(unchecked_bonds[0].atom, current_molecule_atom, branch_sequence, index + 1, branch_atoms)  # check next atom
                            # all branches should either return a function, True, or False.  All child functions should do the same
                            # uncheck_bonds[0].atom becomes new current_molecule_atom, current_molecule_atom becomes previous_molecule_atom
                            # also pass the branch sequence and the index of the next atom_info in the branch
                        else:
                            #cprint("wrong bond or already visited", "red")
                            return False
                            # the next atom in the branch either has the wrong bond or atom symbol
                    else:
                        # there are multiple possible paths branch could go
                        #cprint("checking multiple paths for branch", "yellow")
                        # check all ways
                        for bond in unchecked_bonds:
                            #if current_atom_info.bond != abbr_bond(bond):  # this is purely for seeing what's happening
                                #print(abbr_bond(bond))
                                #print(current_atom_info.bond)
                            if check_bond(bond, current_atom_info.bond, branch_sequence[index + 1]) and bond.atom not in currently_visited:
                                #print(abbr_bond(bond))
                                #print(current_atom_info.bond)
                                # looks at all possible ways that match the correct bond
                                midway_fork = set()
                                # need to separate the branch_atoms here since we don't know if any of the paths will work
                                if check_atom_bonds(bond.atom, current_molecule_atom, branch_sequence, index + 1, midway_fork):
                                    branch_atoms.update(midway_fork)
                                    # if one of the paths works, add all the atoms from the midway_fork "branch"
                                    return True
                                    # return True if any of the paths work (also returns first found)
                                else:
                                    for a in midway_fork:
                                        currently_visited.pop(a)
                                    # if midway_fork doesn't work, need to remove those atoms from currently_visited
                        return False
                        # if this else clause didn't return True (i.e. none of the paths succeeded)
                        # then none of the paths are productive, return false

        if check_branch_point(potential_anchor_atom, None, fragment_map, None):
            #print("phantom atom check")
            for atom in currently_visited:
                if not currently_visited[atom].phantom_atom:
                    # using currently visited to see if the atom_data that was used to find that atom was marked as a phantom_atom
                    # if it is a phantom atom, don't mark as discovered
                    atom.discovered = True
                #else:
                    #print("this atom should not be counted")

            # running checks
            #cprint("number of atoms in fragment: ", "yellow", end="")
            #print(len(molecule_atoms))
            #cprint("number of atoms in currently_visited", "yellow", end="")
            #print(len(currently_visited))
            #if (len(molecule_atoms)) != len(currently_visited):
                #cprint("error in number of atoms found", "red")
            #for atom in molecule_atoms:
                #if atom not in currently_visited:
                    #cprint("error: descrepancy between currently_visited and molecule_atoms", "red")

            #for atom in molecule_atoms:
                #print(atom.symbol)
            #cprint("matched fragment to anchor atom", "magenta")
            return True
        else:
            #cprint("anchor atom not matched to fragment", "red")
            return False
        # start from check_branch point on the potential anchor atom
        # the anchor atom in map is treated as a branch point, even if it only has 1 branch

    fragment_counter = 0
    for atom in potential_anchor_atoms:
        #print("checking anchor atom")
        #for bond in atom.bonded_to:
            #print(abbr_bond(bond))
        if check_anchor_atom(atom, anchored_fragment_map):
            fragment_counter += 1

    #print("\n")
    #cprint("number of fragments found:", "yellow", end="")
    #cprint(fragment_counter, "magenta")

    return fragment_counter  # TODO may want to change the return to something else eventually (the particular atoms???)


# TODO exact match on NC1(CCC2=CC=C3C=C4C(C5C4C5)=CC3=C21)OCCCC(C=C6)=CC=C6C7=CC8=C(C9CC98)C(C%10=C%11C(CCCN%11)=CC(C=O)=C%10)=C7C(CC(C)C)C


def fragmentize(molecule_string, *fragment_libraries, numeric=False):

    molecular_structure = convert_to_structure(MoleculeStructure(), molecule_string)
    fragments = []
    fragments_counter = []
    for lib in fragment_libraries:
        for frag in lib:
            # print("\n\n")
            # cprint(frag, "magenta")
            if numeric:
                frag_num = 0
                for frag_res_structure in lib[frag]:
                    frag_num += (find_fragment(frag_res_structure, None, molecular_structure))
                fragments_counter.append(frag_num)
            else:
                for frag_res_structure in lib[frag]:
                    for f in range(find_fragment(frag_res_structure, None, molecular_structure)):
                        fragments.append(frag)
    atoms_not_discovered = 0
    for atom in molecular_structure.atom_list:
        if not atom.discovered:
            atoms_not_discovered += 1
    cprint(f"atoms not found: {atoms_not_discovered}", "red")
    if atoms_not_discovered > 0:
        total_frags = 0
        for lib in fragment_libraries:
            total_frags += len(lib)
        return ["NA" for _ in range(total_frags)]  # TODO make this total length of libraries
    if numeric:
        return fragments_counter
    else:
        return fragments


print(fragmentize(r"O=C1C=C(CN)ON1", biomolecules, peptide_amino_acids, heterocycles, arenes, functional_groups, hydrocarbons))




