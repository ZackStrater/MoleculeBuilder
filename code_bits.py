
# fragments library smiles alterations functions
def remove_r_groups():
    from fragments_building_library import heterocycles, functionalized_arenes, functional_groups, hydrocarbons, amines, linkers, \
        amino_acids
    import re

    t = open("test_file.txt", "w")
    libs = [heterocycles, functionalized_arenes, functional_groups, hydrocarbons, amines, linkers, amino_acids]
    for i in libs:
        for key in i:
            t.write("\"" + key + "\":")
            t.write(" \"" + re.sub(r"\(\[R\]\)|\[R\]", "", i[key]) + "\",")
            t.write("\n")
        t.write("\n")

    t.close()


def aryl_smiles():
    from fragments_building_library import heterocycles, functionalized_arenes, functional_groups, hydrocarbons, amines, linkers, \
        amino_acids
    import re

    t = open("test_file.txt", "w")
    libs = [heterocycles, functionalized_arenes, functional_groups, hydrocarbons, amines, linkers, amino_acids]
    for i in libs:
        for key in i:
            if "=" in i[key] and "1" in i[key]:
                t.write("\"" + key + "\":")
                a = re.sub(r"\(\[R\]\)|\[R\]", "", i[key])
                b = re.sub(r"=", "", a)
                t.write(" \"" + b.lower() + "\",")
                t.write("\n")
        t.write("\n")

    t.close()



#find fragments hierarchy check
def hierarchy_check(*libraries):  # TODO need to move this somewhere
    import itertools
    from termcolor import cprint
    hierarchy_list = []

    def check_func(index, library, chem_name, structures):
        truncated_lib = itertools.islice(library.items(), index)
        for k, v in truncated_lib:
            for c in v:
                for z in structures:
                    if find_fragment(c, z) > 0:
                        hierarchy_list.append((chem_name, k))

    total_library = {}
    for lib in libraries:
        total_library.update(lib)
    for i, (key, value) in enumerate(total_library.items()):
        cprint("new entry:", "blue")
        print(i, key, value)
        cprint("checked entries:", "blue")
        check_func(i, total_library, key, value)
        print("\n")

    print(hierarchy_list)