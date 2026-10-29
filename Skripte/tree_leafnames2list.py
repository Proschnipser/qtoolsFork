from Bio import Phylo
import re
import pandas as pd
import sys

regex = re.compile("(?<=_)\d+")
df = pd.read_csv("/data/joscha/Data/uniprot_2025_01_08_MIAs_02_MOTHs.csv", index_col='no.')
filepath=sys.argv[1]
tree = Phylo.read(filepath, "newick")

taxa = [leaf.name for leaf in tree.get_terminals()]




output=filepath.replace(".txt", "_taxanames.txt")

with open(output,"a") as f:
    for name in taxa:
        m=re.search(regex, name)
        f.write(re.split(r"[.(]",df["OS"][int(m.group())])[0]+ "\n")
