fragments = {}

fragments["indole"] = "C12=CC=CC=C1NC=C2"
fragments["benzimidazole"] = "C12=CC=CC=C1NC=N2"
fragments["imidazole"] = "C1=CN=CN1"
fragments["benzene"] = "C1=CC=CC=C1"
fragments["naphthalene"] = "C12=CC=CC=C1C=CC=C2"
fragments["cyclohexyl"] = "C1CCCCC1"
fragments["methyl"] = "C"
fragments["ethyl"] = "CC"
fragments["butyl"] = "CCCC"
fragments["pyrrolidine"] = "C1CCNC1"
fragments["cyclopropyl"] = "C1CC1"

# program these with R groups
# write smiles_to_structure addendum that interprets R group as turn can bond to True
# under build_fragments, when bond get's added to, turn can bond to false if atom.bonded_to >= 3
# (that way primary atoms will be able to bond twice)

