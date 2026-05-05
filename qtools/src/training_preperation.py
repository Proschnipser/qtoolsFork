import pandas as pd
import  qtools.data_prepper as dp
import itertools
from pathlib import Path
from random import sample

treefile= "/data/joscha/Downloads/SRw3UZCwUl830gEIOhHRkw_newick.tree"
ref_tree = dp.read_tree(treefile)
outfile_prefix = '/data/joscha/output/qtools/'+ str(Path(treefile).stem).replace(".tree","")+"/"

Path(outfile_prefix).mkdir(parents=True, exist_ok=True)

# prepare names of species used for testing
train_names = sample([leaf.name for leaf in ref_tree.get_terminals()],10)
#_(\d+)_

print(train_names)
# prune the reference tree to names from testing and write to file 
tree = dp.prune_tree(ref_tree, train_names)
dp.write_tree(tree, outfile_prefix+'tree.ph')

minibatches= dp.get_quartets_from_tree(tree)
print(minibatches)
minibatches.to_csv(outfile_prefix  + 'minibatches.csv')

siamesebatches = itertools.combinations(train_names, 2)
pd.DataFrame(siamesebatches).to_csv(outfile_prefix  + 'siamesebatches.csv')

distances = dp.get_edgelength(tree, train_names)
distances.to_csv(outfile_prefix  + 'patristic.csv')