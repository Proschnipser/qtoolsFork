#!/usr/bin/env python3
from Bio import Phylo
import qtools as qt
import qtools.data_prepper as dp
import itertools
import pandas as pd
from random import sample
import numpy as np
from scipy.spatial import distance_matrix
from qtools.quartettroutines import siamesemodel, quartetmodel
from qtools.lossfunctions import  siamloss_siamnet, Xsq_SiamReg, siamloss, eloss

from tensorflow.keras.optimizers import Nadam

ref_tree = Phylo.read("/data/joscha/Downloads/SRw3UZCwUl830gEIOhHRkw_newick.tree","newick")

# prepare names of species used for testing
train_names = sample([leaf.name for leaf in ref_tree.get_terminals()],10)


print(train_names)
# prune the reference tree to names from testing and write to file 
tree = dp.prune_tree(ref_tree, train_names)

minibatches= dp.get_quartets_from_tree(tree)
print(minibatches)
siamesebatches = itertools.combinations(train_names, 2)
siamesebatches = pd.DataFrame(siamesebatches)

distances = dp.get_edgelength(tree, train_names)

batch_size=1
epochs = 100
mode='quartet'

# # # if running as siamese model 
# # if mode == 'siamese':
# #     Multimodel = siamesemodel
# #     minibatch_file = minibatches_file_siamese
# #     sigma = 'nan'
# #     loss_function = siamloss_siamnet
# #     metrics = None

# if running as quartet model, with or without siamese regulation
if mode == 'quartet':
    Multimodel = quartetmodel
    minibatch_file = minibatches_file_quartet 
    sigma = 0.1 
    loss_function = Xsq_SiamReg(sigma)
    metrics = [eloss, siamloss] 