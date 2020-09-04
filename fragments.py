fragments = {}

#fragments["indole"] = "[R]N1C2=C([R])C([R])=C([R])C([R])=C2C([R])=C1[R]"
#fragments["benzimidazole"] = "[R]N1C2=C(C([R])=C(C([R])=C2N=C1[R])[R])[R]"
#fragments["imidazole"] = "[R]C1=NC([R])=C([R])N1[R]"
#fragments["benzene"] = "[R]C1=C([R])C([R])=C([R])C([R])=C1[R]"
#fragments["naphthalene"] = "[R]C1=C2C(C([R])=C([R])C([R])=C2[R])=C([R])C([R])=C1[R]"
#fragments["cyclohexyl"] = "[R]C1C([R])C([R])C([R])C([R])C1[R]"
#fragments["methyl"] = "[R]C"
#fragments["ethyl"] = "[R]CC([R])"
fragments["butyl"] = "[R]CCCC([R])"
#fragments["pyrrolidine"] = "[R]C1C([R])C([R])N([R])C1([R])"
#fragments["cyclopropyl"] = "[R]C1C([R])C1([R])"
fragments["methylamine"] = "[R]NC[R]"

# program these with R groups
# write smiles_to_structure addendum that interprets R group as turn can bond to True
# under build_fragments, when bond get's added to, turn can bond to false if atom.bonded_to >= 3
# (that way primary atoms will be able to bond twice)

