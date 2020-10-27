
# TODO need to introduce a system to check anonymous heterocycles and polycyclic systems
# TODO have an element that can represent any element, have bond that can represent any bond

# TODO maybe add system where atom has a fragemnt attribute (so you can tell if a nitrogen is conjugated for example)

# TODO maybe add a funciton where it only checks lower case structure if it detects multiple lower case atom symbols?

# TODO add special instructions for certain structures
# TODO could replace current system of phantom bonds
# TODO could tell it to make sure X atom has however many bonds connected to whatever
# could also check for CIS/Trans, and check for conjugation

# TODO may need to distinguish between susbistuted amin heterocycles (i.e. methyl piperidine and piperidine)

# TODO add hierarchy function to library (searches all fragments above current fragment to see if any of them a a subgraph of current fragemnt)

amino_acids = {
    "Methionine": ["NC(C=O)CCSC"],
    "Tyrosine": ["NC(C=O)CC1=CC=C(O)C=C1"],
    "Phenylalanine": ["NC(C=O)CC1=CC=CC=C1", "NC(C=O)Cc1ccccc1"],
    "Tryptophan": ["NC(CC1=CNC2=CC=CC=C12)C=O", "NC(C=O)CC1=CNC2=C1C=CC=C2", "NC(Cc1cnc2ccccc12)C=O"],
    "Asparagine": ["NC(C=O)CC(N)=O"],
    "Glutamine": ["NC(C=O)CCC(N)=O"],
    "Cysteine": ["NC(C=O)CS"],
    "Lysine": ["NC(C=O)CCCCN"],
    "Aspartic Acid": ["NC(C=O)CC(O)=O"],
    "Glutamic Acid": ["NC(C=O)CCC(O)=O"],
    "Histidine": ["NC(C=O)CC1=CNC=N1", "NC(CC1=CN=CN1)C=O", "NC(C=O)Cc1cncc1"],
    "Arginine": ["NC(C=O)CCCNC(N)=N"],
    "Tert-Leucine": ["NC(C=O)C(C)(C)C"],
    "Cyclopropyl-Glycine": ["NC(C=O)C1CC1"],
    "Pyroglutamic acid": ["C(C(CC1)NC1=O)=O"],
    "Proline": ["C(C1NCCC1)=O"],
    "Isoleucine": ["NC(C(CC)C)C=O"],
    "Valine": ["NC(C=O)C(C)C"],
    "Leucine": ["NC(C=O)CC(C)C"],
    "Threonine": ["NC(C=O)[C@@H](C)O"],
    "Serine": ["NC(C=O)CO"],
    "Alanine": ["NC(C=O)CW"],
    "Glycine": ["NCC=O"],
    }

heterocycles = {
    "quaternary ammonium": ["NZ"],
    "hantzsch ester": ["CC(NC(C)=C1C(OX)=O)=C(C1)C(OX)=O"],
    "phthalimide": ["O=C1NC(C2=CC=CC=C21)=O", "O=C1NC(C2=C1C=CC=C2)=O"],
    "1,4-dioxane": ["C1COCCO1"],
    "morpholine": ["C1COCCN1"],
    "cytosine": ["O=C1NC=CC(N)=N1"],
    "Tetrahydrothiophene": ["C1CCCS1"],
    "thiophene": ["C1=CC=CS1", "c1cccs1"],
    "uracil": ["O=C(NC=C1)NC1=O"],
    "oxazole": ["C1=NC=CO1", "c1ncco1"],
    "thiazole": ["C1=CN=CS1", "c1cncs1"],
    "indoline": ["N1C2=CC=CC=C2CC1", "C1(CCN2)=C2C=CC=C1", "N1c2ccccc2CC1"],
    "pyrrolizidine": ["C1CCN2C1CCC2"],
    "indole": ["N1C2=CC=CC=C2C=C1", "C1(C=CN2)=C2C=CC=C1", "n1c2ccccc2cc1"],
    "adenine": ["NC1=C2C(N=CN2)=NC=N1", "NC1=NC=NC2=C1NC=N2", "NC1=C2C(NC=N2)=NC=N1", "NC1=NC=NC2=C1N=CN2", "Nc1ncnc2ncnc12"],
    "purine": ["C1=NC=C2C(N=CN2)=N1", "C1(N=CN2)=C2C=NC=N1", "C12=CN=CN=C1NC=N2", "C12=C(C=NC=N2)N=CN1", "c1ncc2c(ncn2)n1"],
    "benzofuran": ["C1=C2C(C=CO2)=CC=C1", "C1(C=CO2)=C2C=CC=C1", "c1c2c(cco2)ccc1"],
    "benzothiophene": ["C1=CC2=CC=CC=C2S1", "C1(C=CS2)=C2C=CC=C1", "c1cc2ccccc2s1"],
    "benzoisoxazole": ["C1=C(ON=C2)C2=CC=C1", "C12=C(C=NO2)C=CC=C1", "c1c(onc2)c2ccc1"],
    "benzoisothiazole": ["C1=NSC2=CC=CC=C21", "C1(C=NS2)=C2C=CC=C1", "c1nsc2ccccc21"],
    "benzoxazole": ["C1=NC2=CC=CC=C2O1", "C1(N=CO2)=C2C=CC=C1", "c1nc2ccccc2o1"],
    "benzothiazole": ["C1=NC2=CC=CC=C2S1", "C1(N=CS2)=C2C=CC=C1", "c1nc2ccccc2s1"],
    "isoxazole": ["C1=NOC=C1", "c1nocc1"],
    "quinoline": ["C1=C2C(C=CC=N2)=CC=C1", "C1(C=CC=N2)=C2C=CC=C1", "c1c2c(cccn2)ccc1"],
    "isoquinoline": ["C1=CC=C(C=CN=C2)C2=C1", "C1(C=CN=C2)=C2C=CC=C1", "c1ccc(ccnc2)c2c1"],
    "quinoxaline": ["C1=NC2=CC=CC=C2N=C1", "C1(N=CC=N2)=C2C=CC=C1", "c1nc2ccccc2nc1"],
    "quinazoline": ["C1=NC2=CC=CC=C2C=N1", "C1(C=NC=N2)=C2C=CC=C1", "c1nc2ccccc2cn1"],
    "cinnoline": ["C1=CC2=CC=CC=C2N=N1", "C1(C=CN=N2)=C2C=CC=C1", "c1cc2ccccc2nn1"],
    "pteridine": ["C1=NC2=NC=CN=C2C=N1", "C1(C=NC=N2)=C2N=CC=N1", "c1nc2nccnc2cn1"],
    "chromenone": ["O=C1C=COC2=CC=CC=C21", "O=C(C=CO1)C2=C1C=CC=C2"],
    "chromene": ["C1=COC2=CC=CC=C2C1", "C1(CC=CO2)=C2C=CC=C1", "c1(CC=CO2)c2cccc1"],
    "coumarin": ["O=C1C=CC2=CC=CC=C2O1", "O=C(O1)C=CC2=C1C=CC=C2", "O=C(C=C1)Oc2c1cccc2"],
    "quinolinone": ["O=C1C=CC2=CC=CC=C2N1", "O=C(N1)C=CC2=C1C=CC=C2", "O=C(N1)C=Cc2c1cccc2"],
    "isoquinolinone": ["O=C1NC=CC2=CC=CC=C21", "O=C1NC=CC2=C1C=CC=C2", "O=C1NC=Cc2c1cccc2"],
    "caprolactam": ["O=C1NCCCCC1"],
    "piperazine": ["N1CCNCC1"],
    "hydantoin": ["O=C1NCC(N1)=O"],
    "tetrazole": ["C1=NN=NN1", "C1=NNN=N1", "c1nnnn1"],
    "2-oxazolidone": ["O=C1NCCO1"],
    "guanine": ["O=C1C2C(NC=N2)N=C(N)N1", "O=C1C2C(N=C(N1)N)N=CN2", "O=C1C2C(NC(N)=N1)NC=N2", "O=C1C2C(NC(N)=N1)N=CN2"],
    "pyrolidinone": ["O=C1NCCC1"],
    "gamme-butyrrolactone": ["O=C1OCCC1"],
    "beta-lactone": ["O=C1CCO1"],
    "beta-lactam": ["O=C1CCN1"],
    "4-piperidinone": ["O=C1CCNCC1"],
    "2-piperidinone": ["O=C1NCCCC1"],
    "pyrrolidine": ["N1CCCC1"],
    "pyrrole": ["N1C=CC=C1", "n1cccc1"],
    "pyrazolidine": ["N1NCCC1"],
    "piperidine": ["C1CCNCC1"],
    "pyridine": ["C1=CC=CC=N1", "c1ccccn1"],
    "pyrazole": ["N1C=CC=N1", "n1cccn1"],
    "imidazole": ["C1=NC=CN1", "c1nccn1"],
    "pyridazine": ["C1=NN=CC=C1", "C1=CC=CN=N1", "c1nnccc1"],
    "pyrimidine": ["C1=NC=CC=N1", "c1ncccn1"],
    "pyrazine": ["C1=CN=CC=N1", "c1cnccn1"],
    "1,2,4-traizole": ["C1=NNC=N1", "C1=NN=CN1", "c1nncn1"],
    "1,2,3-triazole": ["N1N=NC=C1", "n1nncc1"],
    "1,2,3-triazine": ["C1=NN=CC=N1", "N1=NC=NC=C1", "c1nnccn1"],
    "1,3,5-triazine": ["C1=NC=NC=N1", "c1ncncn1"],
    "furan": ["C1=CC=CO1", "c1ccco1"],
    "chromane": ["C12=CC=CC=C1OCCC2", "C1(CCCO2)=C2C=CC=C1", "c12ccccc1OCCC2"],
    }

arenes = {
    "anthracene": ["C12=CC=CC=C1C=C3C(C=CC=C3)=C2", "C1(C=C(C=CC=C2)C2=C3)=C3C=CC=C1", "c12ccccc1cc3c(cccc3)c2"],
    "phenanthrene": ["C12=CC=CC=C1C3=CC=CC=C3C=C2", "C12=CC=CC=C1C(C=CC=C3)=C3C=C2", "c12ccccc1c3ccccc3cc2"],
    "fluorene": ["C12=CC=CC=C1CC3=CC=CC=C23", "C1(C(C=CC=C2)=C2C3)=C3C=CC=C1", "C12=CC=CC=C1CC3=C2C=CC=C3", "c12c(c3c(C2)cccc3)cccc1"],
    "naphthalene": ["C12=CC=CC=C1C=CC=C2", "C1(C=CC=C2)=C2C=CC=C1", "c1ccc2c(cccc2)c1"],
    "benzene": ["C1=CC=CC=C1", "c1ccccc1"],
    }

functional_groups = {
    "tetrahydropyranone": ["O=C1CCOCC1"],
    "tetrahydropyran": ["C1CCCCO1"],
    "tetrahydrofuran": ["C1CCCO1"],
    "cyclohexanone": ["O=C1CCCCC1"],
    "cyclohexenone": ["O=C1C=CCCC1"],
    "enone": ["C=CC=O"],
    "enol": ["C=COW"],
    "nitro": ["ON(=O)"],
    "1,3-dioxolane": ["C1OCCO1"],
    "dimethylamine": ["WCNCW"],
    "diethylamine": ["WC(CY)N(CY)CW"],
    "Urea": ["O=C(N)N"],
    "Thiourea": ["S=C(N)N"],
    "carbamate": ["O=C(N)O"],
    "Cyano": ["C#N"],
    "trifluoromethyl": ["C(F)(F)F"],
    "difluoromethyl": ["C(F)F"],
    "sulfonamide": ["O=S(N)=O"],
    "carboxylic acid": ["O=COW"],
    "acetoxy": ["WCC(OX)=O"],
    "methyl ester": ["O=COCW"],
    "ethyl ester": ["O=COC(X)CW"],
    "ester": ["O=COX"],
    "amide": ["NC(=O)Y"],
    "sulfoxide": ["O=SY"],
    "thiol": ["SW"],
    "thioether": ["SX"],
    "vinyl": ["C=C"],
    "alkyne": ["C#C"],
    "sulfonic acid": ["O=S(O)=O"],
    "Sulfone": ["O=S=O"],
    "guanidine": ["N=C(N)N"],
    "amidine": ["N=CN"],
    "imine": ["N=C"],
    "epoxide": ["C1OC1"],
    "cyclopropane": ["C1CC1"],
    "isopropyl": ["YC(CW)(CW)"],
    "tert-butyl": ["YC(CW)(CW)(CW)"],
    "aldehyde": ["O=CX"],
    "ketone": ["O=CY"],
    "methoxy": ["WCOX"],
    "ether": ["OX"],
    "hydroxyl": ["OW"],
    "oxo": ["O"],
    "primary amine": ["NW"],
    "secondary amine": ["NX"],
    "tertiary amine": ["NY"],
    "fluoro": ["F"],
    "chloro": ["Cl"],
    "bromo": ["Br"],
    "iodo": ["I"],
    }

hydrocarbons = {
    "1,3-cyclohexadiene": ["C1C=CC=CC1"],
    "1,4-cyclohexadiene": ["C1=CCC=CC1"],
    "cyclohexene": ["C1CCC=CC1"],
    "1,3-cyclopentadiene": ["C1C=CC=C1"],
    "cyclopentene": ["C1CCC=C1"],
    "cyclooctane": ["C1CCCCCCC1"],
    "cyclohpetane": ["C1CCCCCC1"],
    "cyclohexane": ["C1CCCCC1"],
    "cyclopentane": ["C1CCCC1"],
    "cyclobutane": ["C1CCC1"],
    "addamantyl": ["C1(CC2C3)CC(C2)CC3C1"],
    "octyl": ["CCCCCCCC"],
    "septyl": ["CCCCCCC"],
    "hexyl": ["CCCCCC"],
    "pentyl": ["CCCCC"],
    "butyl": ["CCCC"],
    "propyl": ["CCC"],
    "ethyl": ["CC"],
    "methyl": ["CW"],
    "methylene": ["CX"],
    "tertiary carbon": ["CY"],
    "quaternary carbon": ["CZ"],
    }





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

#aryl_smiles()

