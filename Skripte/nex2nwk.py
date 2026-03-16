from Bio import Phylo
import dendropy
import io
import os
import sys

filepath=sys.argv[1]
# Read the MrBayes NEXUS consensus tree
output= os.path.splitext(filepath)[0]+".nwk"
print(output)
tree = dendropy.Tree.get(
    path=filepath,
    schema="nexus"
)

tree.write(path=output, schema="newick")