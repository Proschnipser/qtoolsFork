#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
from tensorflow.keras.optimizers import Nadam
import qtools as qt
from qtools.lossfunctions import  siamloss_siamnet, Xsq_SiamReg, siamloss, eloss
from qtools.quartettroutines import siamesemodel, quartetmodel
from qtools.data_tracking import metadata, mutationscheme, Tracking

treefile= "/data/joscha/Downloads/SRw3UZCwUl830gEIOhHRkw_newick.tree"
outfile_prefix = '/data/joscha/output/qtools/'+ str(Path(treefile).stem).replace(".tree","")+"/"

batch_size=1
epochs = 100
mode='quartet'

# file with training sequences 
vector_file = ""


# file with edge length distance matrix 
edge_distance = outfile_prefix + 'patristic.csv'

# file with siamese minibatches 
minibatches_siamese = outfile_prefix + 'siamesebatches.csv' 


# file with quartet minibatches 
minibatches_quartet = outfile_prefix + 'minibatches.csv'


# if running as siamese model 
if mode == 'siamese':
    Multimodel = siamesemodel
    minibatch_file = minibatches_file_siamese
    sigma = 'nan'
    loss_function = siamloss_siamnet
    metrics = None

# if running as quartet model, with or without siamese regulation
if mode == 'quartet':
    Multimodel = quartetmodel
    minibatch_file = minibatches_file_quartet 
    sigma = 0.1 
    loss_function = Xsq_SiamReg(sigma)
    metrics = [eloss, siamloss] 

data = pd.read_csv(vector_file)
data = qt.qdata(data)






for e in range(epochs):
	# prepare training batches
    batches = data.batchmaker(minibatches, batch_size, edge_distance) 
           
    # train quartetnet (or siamesenet)
    multimodel.fit(batches, batch_size = batch_size)         

    # evaluate epoche
    losses = multimodel.history.history
    
    # calculate matrix with euclidean distances
    prediction = multimodel.predict(x_encoded)
    matrix_i = multimodel.get_distance_matrix(prediction)
    
    # calculate quartet scores (check how many quartets are in right split)
    scores = qt.get_qscores(matrix_i, x_species, scoring_batches)
    
    # make splitstree diagram from distance matrix
    qt.matrix2nexus(matrix=matrix_i, taxa=x_species, nexusfile='filename.nex', plot_now=True)   