from ete3 import Tree
import sys

filepath="/data/joscha/Data/TANGO1vsMIAstrict_dedup_genus_names_thresh100aa_150alt.nex.con_annotatedOC.tree"#sys.argv[1]/bioinf/data/user/joscha/Data/TANGO1vsMIAstrict_dedup_genus_names_thresh100aa_150alt.nex.con_annotatedOC.tree
tree=Tree(filepath,format=1)
#print(tree)
parent=tree.search_nodes(name="'48 'Metazoa |Ecdysozoa 79.3103448275862% Spiralia 18.96551724137931% Chordata 1.7241379310344827% |''")
print(parent)

tree.set_outgroup(parent[0])

output= filepath.replace(".tree", "_root.tree")
print(output)
tree.write(outfile=output, format=1)