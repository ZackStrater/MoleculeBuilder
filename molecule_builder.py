

class Atom:

    def __init__(self, symbol):
        self.symbol = symbol
        self.bonded_to = []
        self.can_bond = False


class Bond:

    def __init__(self, atom, bond_code):
        self.atom = atom
        self.bond_code = bond_code


bond_decoder = {   # TODO add chiral stereochemistry
    1: "",  # single bond
    2: "=",  # double bond
    3: "#",  # triple bond
    4: "",  # aromatic bond
    5: "$",  # ionic bond
    6: "/",  # vinyl up
    7: "\\",  # vinyl down
    8: ".",  # non-bond
}


def decode_bond(code):
    return bond_decoder[code]


class AtomSmiles:
    # data object for atoms.  keeps track of bonding information, ring closure bond info, and closure partners

    def __init__(self):
        self.bond_code = None
        self.ring_closures = []
        self.closure_partners = []


class Molecule:

    def __init__(self, seed_frag):
        self.seed_frag = seed_frag
        self.atom_list = []
        self.atom_list.extend(seed_frag.atom_list)

    def manual_add_frag(self, fragment, frag_atom, base_atom):
        self.atom_list.extend(fragment.atom_list)
        base_atom.bonded_to.append(Bond(frag_atom, 1))
        frag_atom.bonded_to.append(Bond(base_atom, 1))

    def manual_2point_add_frag(self, fragment, frag_atom1, frag_atom2, base_atom1, base_atom2):
        self.atom_list.extend(fragment.atom_list)
        base_atom1.bonded_to.append(Bond(frag_atom1, 1))
        base_atom2.bonded_to.append(Bond(frag_atom2, 1))
        frag_atom1.bonded_to.append(Bond(base_atom1, 1))
        frag_atom2.bonded_to.append(Bond(base_atom2, 1))


def convert_to_smiles(molecule, start):
    smiles_string = ""
    smiles_construction_list = []
    # ordered list of atoms and parentheses
    atom_info = {}
    # dictionary linking atoms to their AtomSmiles info
    visited = {}
    # keeping track of which atoms have already been visited
    ring_counter = 1
    # keeps track and numbers rings
    for a in molecule.atom_list:
        visited[a] = False
    # all atoms start unvisited

    def dfs(current_atom, previous_atom):
        nonlocal ring_counter
        visited[current_atom] = True
        smiles_construction_list.append(current_atom)

        previous_bond = None
        for bond in current_atom.bonded_to:
            if bond.atom == previous_atom:
                previous_bond = bond
        # this algorithm is backward looking, i.e. bond info

        atom_info[current_atom] = AtomSmiles()
        if previous_bond:
            atom_info[current_atom].bond_code = previous_bond.bond_code

        # ring detection
        for bond in current_atom.bonded_to:
            if bond != previous_bond and visited[bond.atom] and bond.atom not in atom_info[current_atom].closure_partners:
                # ring detected if the bond was not the bond we just came from
                # atom has been visited already
                # and atom has not previously marked as a closure partner with current atom (prevents double counting rings)
                atom_info[current_atom].ring_closures.append(decode_bond(bond.bond_code))
                if ring_counter > 9:
                    atom_info[current_atom].ring_closures.append("%")
                    atom_info[bond.atom].ring_closures.append("%")
                atom_info[current_atom].ring_closures.append(ring_counter)
                atom_info[bond.atom].ring_closures.append(ring_counter)
                # formatting (bond info)(% if ring # > 9)(ring #)
                # bonding info only applied to current atom, not both closure partners (accepted smiles notation)
                atom_info[current_atom].closure_partners.append(bond.atom)
                atom_info[bond.atom].closure_partners.append(current_atom)

                ring_counter += 1

            untraveled_paths = []
            for b in current_atom.bonded_to:
                if not visited[b.atom]:
                    untraveled_paths.append(bond.atom)
            if not visited[bond.atom]:
                if len(untraveled_paths) > 1:
                    smiles_construction_list.append("(")
                dfs(bond.atom, current_atom)
                if len(untraveled_paths) > 1:
                    smiles_construction_list.append(")")
                # parentheses added to denote fragments attached to current atom
                # parentheses not needed if only one path is left

    dfs(start, None)
    # begin dfs process

    # does this even do anything? maybe trying to reorder ring numbers so they go in order?
    atom_ring_list = []
    for ele in smiles_construction_list:
        if isinstance(ele, Atom):
            if len(atom_info[ele].ring_closures) > 0:
                atom_ring_list.append(smiles_construction_list.index(ele))

    # building the smiles string
    for ele in smiles_construction_list:
        if isinstance(ele, str):
            smiles_string += ele
        else:
            if atom_info[ele].bond_code:
                smiles_string += decode_bond(atom_info[ele].bond_code)
            smiles_string += ele.symbol
            for c in atom_info[ele].ring_closures:
                smiles_string += str(c)

    print(smiles_string)


class Benzene:

    def __init__(self):
        self.atom_list = []
        self.c1 = Atom("C")
        self.c2 = Atom("C")
        self.c3 = Atom("C")
        self.c4 = Atom("C")
        self.c5 = Atom("C")
        self.c6 = Atom("C")
        self.atom_list.extend([self.c1, self.c2, self.c3, self.c4, self.c5, self.c6])
        self.c1.bonded_to.append(Bond(self.c2, 2))
        self.c2.bonded_to.append(Bond(self.c3, 1))
        self.c3.bonded_to.append(Bond(self.c4, 2))
        self.c4.bonded_to.append(Bond(self.c5, 1))
        self.c5.bonded_to.append(Bond(self.c6, 2))
        self.c6.bonded_to.append(Bond(self.c1, 1))
        self.c2.bonded_to.append(Bond(self.c1, 2))
        self.c3.bonded_to.append(Bond(self.c2, 1))
        self.c4.bonded_to.append(Bond(self.c3, 2))
        self.c5.bonded_to.append(Bond(self.c4, 1))
        self.c6.bonded_to.append(Bond(self.c5, 2))
        self.c1.bonded_to.append(Bond(self.c6, 1))


class Naphthalene:

    def __init__(self):
        self.atom_list = []
        self.c1 = Atom("C")
        self.c2 = Atom("C")
        self.c3 = Atom("C")
        self.c4 = Atom("C")
        self.c5 = Atom("C")
        self.c6 = Atom("C")
        self.c7 = Atom("C")
        self.c8 = Atom("C")
        self.c9 = Atom("C")
        self.c10 = Atom("C")
        self.atom_list.extend([self.c1, self.c2, self.c3, self.c4, self.c5, self.c6, self.c7, self.c8, self.c9, self.c10])
        self.c1.bonded_to.append(Bond(self.c2, 1))
        self.c2.bonded_to.append(Bond(self.c3, 2))
        self.c3.bonded_to.append(Bond(self.c4, 1))
        self.c4.bonded_to.append(Bond(self.c5, 2))
        self.c5.bonded_to.append(Bond(self.c6, 1))
        self.c6.bonded_to.append(Bond(self.c7, 2))
        self.c7.bonded_to.append(Bond(self.c8, 1))
        self.c8.bonded_to.append(Bond(self.c9, 2))
        self.c9.bonded_to.append(Bond(self.c10, 1))
        self.c10.bonded_to.append(Bond(self.c1, 2))
        self.c2.bonded_to.append(Bond(self.c1, 1))
        self.c3.bonded_to.append(Bond(self.c2, 2))
        self.c4.bonded_to.append(Bond(self.c3, 1))
        self.c5.bonded_to.append(Bond(self.c4, 2))
        self.c6.bonded_to.append(Bond(self.c5, 1))
        self.c7.bonded_to.append(Bond(self.c6, 2))
        self.c8.bonded_to.append(Bond(self.c7, 1))
        self.c9.bonded_to.append(Bond(self.c8, 2))
        self.c10.bonded_to.append(Bond(self.c9, 1))
        self.c1.bonded_to.append(Bond(self.c10, 2))
        self.c10.bonded_to.append(Bond(self.c5, 1))
        self.c5.bonded_to.append(Bond(self.c10, 1))


class Propyl:

    def __init__(self):
        self.atom_list = []
        self.c1 = Atom("C")
        self.c2 = Atom("C")
        self.c3 = Atom("C")
        self.atom_list.extend([self.c1, self.c2, self.c3])
        self.c1.bonded_to.append(Bond(self.c2, 1))
        self.c2.bonded_to.append(Bond(self.c3, 1))
        self.c2.bonded_to.append(Bond(self.c1, 1))
        self.c3.bonded_to.append(Bond(self.c2, 1))


class Aminyl:

    def __init__(self):
        self.atom_list = []
        self.n1 = Atom("N")
        self.atom_list.extend([self.n1])


class Hydroxyl:

    def __init__(self):
        self.atom_list = []
        self.o1 = Atom("O")
        self.atom_list.extend([self.o1])


class Formyl:

    def __init__(self):
        self.atom_list = []
        self.c1 = Atom("C")
        self.o1 = Atom("O")
        self.atom_list.extend([self.c1, self.o1])
        self.c1.bonded_to.append(Bond(self.o1, 2))
        self.o1.bonded_to.append(Bond(self.c1, 2))


class Cyclopropyl:

    def __init__(self):
        self.atom_list = []
        self.c1 = Atom("C")
        self.c2 = Atom("C")
        self.c3 = Atom("C")
        self.atom_list.extend([self.c1, self.c2, self.c3])
        self.c1.bonded_to.append(Bond(self.c2, 1))
        self.c2.bonded_to.append(Bond(self.c3, 1))
        self.c3.bonded_to.append(Bond(self.c1, 1))
        self.c2.bonded_to.append(Bond(self.c1, 1))
        self.c3.bonded_to.append(Bond(self.c2, 1))
        self.c1.bonded_to.append(Bond(self.c3, 1))


benzene1 = Benzene()
benzene2 = Benzene()
benzene3 = Benzene()
propyl1 = Propyl()
propyl2 = Propyl()
propyl3 = Propyl()
propyl4 = Propyl()
propyl5 = Propyl()
formyl1 = Formyl()
aminyl1 = Aminyl()
hydroxyl1 = Hydroxyl()
naphthalene1 = Naphthalene()
cyclopropyl1 = Cyclopropyl()
cyclopropyl2 = Cyclopropyl()
aminyl3 = Aminyl()
molecule1 = Molecule(benzene1)
molecule1.manual_add_frag(benzene2, benzene2.c1, benzene1.c5)
molecule1.manual_add_frag(benzene3, benzene3.c1, benzene1.c3)
molecule1.manual_add_frag(propyl1, propyl1.c2, benzene1.c4)
molecule1.manual_add_frag(propyl2, propyl2.c2, propyl1.c1)
molecule1.manual_add_frag(propyl3, propyl3.c3, benzene2.c4)
molecule1.manual_add_frag(formyl1, formyl1.c1, benzene3.c5)
molecule1.manual_add_frag(aminyl1, aminyl1.n1, benzene3.c2)
molecule1.manual_add_frag(hydroxyl1, hydroxyl1.o1, propyl3.c1)
molecule1.manual_add_frag(propyl4, propyl4.c1, hydroxyl1.o1)
molecule1.manual_2point_add_frag(propyl5, propyl5.c1, propyl5.c3, aminyl1.n1, benzene3.c3)
molecule1.manual_2point_add_frag(naphthalene1, naphthalene1.c1, naphthalene1.c2, propyl4.c1, propyl4.c3)
molecule1.manual_2point_add_frag(cyclopropyl1, cyclopropyl1.c1, cyclopropyl1.c2, benzene1.c1, benzene1.c2)
molecule1.manual_2point_add_frag(cyclopropyl2, cyclopropyl2.c1, cyclopropyl2.c2, naphthalene1.c7, naphthalene1.c8)
molecule1.manual_add_frag(aminyl3, aminyl3.n1, propyl4.c1)


convert_to_smiles(molecule1, benzene1.c1)


benzene4 = Benzene()
propyl6 = Propyl()
propyl7 = Propyl()
aminyl3 = Aminyl()
hydroxyl2 = Hydroxyl()
hydroxyl3 = Hydroxyl()
hydroxyl4 = Hydroxyl()
molecule2 = Molecule(benzene4)
molecule2.manual_add_frag(aminyl3, aminyl3.n1, benzene4.c2)
molecule2.manual_add_frag(hydroxyl2, hydroxyl2.o1, benzene4.c5)
molecule2.manual_add_frag(propyl6, propyl6.c1, hydroxyl2.o1)
molecule2.manual_2point_add_frag(propyl7, propyl7.c1, propyl7.c3, aminyl3.n1, propyl6.c3)
molecule2.manual_add_frag(hydroxyl3, hydroxyl3.o1, propyl6.c3)
molecule2.manual_add_frag(hydroxyl4, hydroxyl4.o1, propyl6.c3)

convert_to_smiles(molecule2, benzene4.c4)





# things to program


# parens get added when you enter a ring, but this  is unecessary
# could make smiles string simpler - > when adding fragment, put fragment bond at beginning of atom.bonded to
# make a global map for bonds i.e 2 -> = 3 -> # etc...

# bicyclic fragments
# two point fragment addition
# triple bonds
# add attribute can_bond to Atom class (maybe tracks valency)
# create auto fragment addition
# E vs Z stereochem with \/
# add chirality with @ and @@
    # cares about fragment order

# change ring numbering system in convert_to_smiles
    # maybe have dictionary for visited ring fragments
        # when it's the first time you visit an atom in a ring fragment, mark in dictionary as visited and mark atom
        # then either mark adjacent atom, or mark atom once number of times fragment has been visited == num of atoms


