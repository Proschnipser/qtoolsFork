from ete3 import Tree
import sys

filepath="/data/joscha/Data/TANGO1vsMIAvsOTORstrict_dedup_genus_names_thresh100aa_400.nex.con_annotatedname.tree"#sys.argv[1]
tree=Tree(filepath,format=1)

parent=tree.search_nodes(name="'TANGO1 60.64814814814815% MIA 39.351851851851855% '")[0]
print(parent.children)

tree.set_outgroup(parent.children[0].children[0])

output= filepath.replace(".tree", "_root.tree")
print(output)
tree.write(outfile=output, format=1)