from Smiles_to_Structure import convert_to_structure, MoleculeStructure
from fragments import heterocycles
from collections import Counter


class FragmentAtomData:
    def __init__(self, atom, symbol):  # TODO is adding object here necessary????
        self.atom = atom
        self.symbol = symbol
        self.bonds = []
        self.daughter_branches = []  #TODO add child branch to atoms


class Branch:
    def __init__(self, branch_num):
        self.branch_num = branch_num
        self.sequence = []


class FragmentBranches:
    def __init__(self):
        self.map = {}


def abbr_bond(bond):
    return bond.bond_code, bond.atom.symbol
# tuple representation of a bond

def fragmentize(smiles_string, fragment_string):  #specify fragment list

    total_molecule = convert_to_structure(MoleculeStructure(), smiles_string)
    current_fragment = convert_to_structure(MoleculeStructure(), fragment_string)

    def find_atom(fragment):
        for ele in ["Si", "P", "I", "Br", "Cl", "F", "B", "S", "O", "N", "C"]:
            for atom in fragment.atom_list:
                if atom.symbol == ele:
                    return atom
    base_atom = find_atom(current_fragment)

    def map_fragment(fragment, base_atom):
        fragment_branches = FragmentBranches()
        # keeps track of branch connectivity
        visited = {}
        for atom in fragment.atom_list:
            visited[atom] = False
        # keeps track of which atoms have been visited
        branch_num = 1  # TODO this might be unnecessary with daughter_branches

        def dfs(current_atom, previous_atom, branch):
            print(current_atom.symbol)
            visited[current_atom] = True

            current_atom_data = FragmentAtomData(current_atom, current_atom.symbol)
            # data object for current atom

            previous_bond = None
            for bond in current_atom.bonded_to:
                if bond.atom == previous_atom:
                    previous_bond = bond

            if len(current_atom.bonded_to) > 2:
                branch.sequence.append(current_atom_data)
                # data objects for atoms added sequentially to branch.sequence
                for bond in current_atom.bonded_to:
                    if bond != previous_bond and not visited[bond.atom]:
                        print("new branch")
                        new_branch(bond.atom, current_atom)
                        # current_atom will be the previous atom for new branch
            # if atom is a branch point, create new branch for each path

            elif len(current_atom.bonded_to) == 2:
                for bond in current_atom.bonded_to:
                    if bond != previous_bond:
                        current_atom_data.bonds.append(abbr_bond(bond))
                branch.sequence.append(current_atom_data)

                for bond in current_atom.bonded_to:
                    if not visited[bond.atom]:
                        dfs(bond.atom, current_atom, branch)
            # elif atom is contiguous point, continue to next atom

            else:
                branch.sequence.append(current_atom_data)
            # else dead end

        def new_branch(current_atom, previous_atom):
            nonlocal branch_num
            current_branch = Branch(branch_num)
            fragment_branches.map[branch_num] = current_branch
            branch_num += 1
            dfs(current_atom, previous_atom, current_branch)
            # start dfs on new branch and add info to fragment_branches


        # TODO maybe try starting with dfs on first atom instead??
        # TODO add ring detection eventually
        new_branch(base_atom, None)
        # start a branch at base atom

        for key in fragment_branches.map:
            print(key)
            branch = fragment_branches.map[key]
            for atom_info in branch.sequence:
                print(atom_info.symbol)
                for bond in atom_info.bonds:
                    print(bond)

    map_fragment(current_fragment, base_atom)

    base_atom_bonds = [abbr_bond(bond) for bond in base_atom.bonded_to]
    # find base atom (i.e. highest priority atom from find_atom) from fragment
    count_base_atom_bonds = Counter(base_atom_bonds)
    # counts number of unique bonds types and how many of each base atom has

    def find_starter_atoms(atom, count_atom_bonds, atom_list):
        # searches through all atoms in molecules in total_molecule to see if they match the fragment base atom
        if atom.symbol == base_atom.symbol:
            # check to see if atom is the same element
            bond_list = [abbr_bond(bond) for bond in atom.bonded_to]
            count_bond_list = Counter(bond_list)
            for key in count_atom_bonds:
                if key not in count_bond_list or count_atom_bonds[key] > count_bond_list[key]:
                    # check 1: are there bonds types in fragment base atom that current atom doesn't have
                    # check 2: does current atom have >= the amount of each bond type compared to fragment base atom
                    # i.e. are the bonds in fragment base atom a subset of the bonds of current atom
                    return
            atom_list.append(atom)
            # if all checks passed, atom is a potential base atom and is  stored in a list

    starting_atoms = []
    # keeping track of atoms that match fragment base atom
    for atom in total_molecule.atom_list:
        find_starter_atoms(atom, count_base_atom_bonds, starting_atoms)


fragmentize("O=C1OC(C2=C(C(N3C)=CC(N4C)=NC3=O)C(C5=NC6=C4C=C(C7=CC=C8C(OC=C8)=C7)C=C6N=C5)=C9N=COC9=C2)CC1C%10=CC(C%11=NC(C%12=C%13C(C(C%14=C%15C=CN=NC%15=CC=C%14)=CC=N%13)=CC(C%16=NC=CC=N%16)=C%12)=NN%11C)=NO%10", "ClCSi1(CCCBr)CC(CCF)CC1")
