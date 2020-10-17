
from fragments import heterocycles, functionalized_arenes, functional_groups, hydrocarbons, amines, linkers, amino_acids
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