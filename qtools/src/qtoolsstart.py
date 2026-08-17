#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib
import ast
import pandas as pd
import numpy as np
from pathlib import Path
from tensorflow.keras.optimizers import Nadam
import qtools as qt
from qtools.lossfunctions import  siamloss_siamnet, Xsq_SiamReg, siamloss, eloss
from qtools.quartettroutines import siamesemodel, quartetmodel
from qtools.data_tracking import metadata, mutationscheme, Tracking
import sys, os, traceback


treefile= "/data/joscha/Downloads/SRw3UZCwUl830gEIOhHRkw_newick.tree"
outfile_prefix = '/data/joscha/output/qtools/'+ str(Path(treefile).stem).replace(".tree","")+"/"
# path for writing results 
out_path = outfile_prefix+ '/trained_models_test/'

# create a new subdir in your out_dir from local time
out_dir = qt.update_dir(out_path)
qt.create_dir(out_dir)
qt.create_dir(out_dir+"/weights/")
features_dir = out_dir+"/features/"
qt.create_dir(features_dir)
# =============================================================================
# setting up training variables 
# =============================================================================


# variables for training 
learning_rate = 0.0001
batch_size = 1
epochs = 5
mode='quartet'

# file with training sequences 
vector_file = outfile_prefix+"vectors.csv"


# file with edge length distance matrix 
edge_distance = outfile_prefix + 'patristic.csv'

# file with siamese minibatches 
minibatches_siamese = outfile_prefix + 'siamesebatches.csv' 


# file with quartet minibatches 
minibatches_quartet = outfile_prefix + 'minibatches.csv'




# if running as siamese model 
if mode == 'siamese':
    Multimodel = siamesemodel
    minibatch_file = minibatches_siamese
    sigma = 'nan'
    loss_function = siamloss_siamnet
    metrics = None

# if running as quartet model, with or without siamese regulation
if mode == 'quartet':
    Multimodel = quartetmodel
    minibatch_file = minibatches_quartet
    sigma = 0.1 
    loss_function = Xsq_SiamReg(sigma)
    metrics = [eloss, siamloss] 

# =============================================================================
# prepare data 
# =============================================================================

# set up training data
data = pd.read_csv(vector_file)
data['mtxvector'] = data['mtxvector'].apply(ast.literal_eval)
print("Data mtxvector:", type(data['mtxvector'][0]), data['mtxvector'][0])
data = qt.qdata(data)

#  read minibatches 
minibatches = pd.read_csv(minibatch_file, index_col=0)
edge_distance = pd.read_csv(edge_distance, index_col=0)
scoring_batches = pd.read_csv(minibatches_quartet, index_col=0)


# =========================================================================
# set up model 
# =========================================================================
x_vector, x_species = data.get_data()
mean_vector = np.mean(x_vector, axis=0) 

print(len(mean_vector))
print(mean_vector)
# set up model 
vec_len = data.get_veclen()
output_dims=5
singlemodel = qt.fully_connected(vec_len,mean_vector,output_dims=output_dims)

# build quartetnet or siamesenet
multimodel = Multimodel.from_basemodel(singlemodel.model)
multimodel.compile(optimizer=Nadam(learning_rate=0.001), loss=loss_function, metrics=metrics)




# track the metadata
# you can inspect which variables are written to metadata with metadata.collected_keys
# you can modify which data are tracked in qtools.data_tracking.metadata
print(locals()['out_dir'])
metadata.record(locals()).write()



# =========================================================================
# evaluate model before first run
# =========================================================================

# get feature vectors
prediction = multimodel.predict(np.array(x_vector))

# calculate the distance matrix (euclidean distances) from the feature vectores
matrix_i = multimodel.get_distance_matrix(prediction)  

# create a splits diagram from the distance matrix with splitstree.
# the default directory for splitstree is '~/splitstree4/SplitsTree' 
# but you can change it with option 'splitstree_location'
qt.matrix2nexus(matrix_i, x_species,  out_dir + 'nexus/0.nex',
                plot_now=True)

# evaluate quartet scores
scores = qt.get_qscores(matrix_i, x_species, scoring_batches)


# evaluate loss
batches = data.batchmaker(minibatches, edge_distance)
losses = multimodel.evaluate(batches, return_dict=True)

# initialize the tracking before the first epoch 
t = Tracking(out_dir, y_names=x_species)
# track evaluation before first run 
#t.trackall(epoch=0, feature_vectors=prediction, score=scores, loss=losses)
#t.writeall()
t.write_species_names()
print("outdir",out_dir)
(Path(out_dir) / "nexus").mkdir(parents=True, exist_ok=True)

for e in range(epochs):
	# prepare training batches
    batches = data.batchmaker(minibatches, edge_distance, batch_size=batch_size) 
           
    # train quartetnet (or siamesenet)
    multimodel.fit(batches, batch_size = batch_size)
    #os.makedirs(os.path.dirname(weights_path), exist_ok=True)
    # save model weights 
    multimodel.basemodel.save_weights(f'{out_dir}/weights/m{e+1}_weights.h5')
    print("Saving to:", Path(f"{out_dir}/weights/m{e+1}_weights.h5").resolve())

    # evaluate epoche
    losses = multimodel.history.history
    
    # calculate matrix with euclidean distances
    prediction = multimodel.predict(x_vector)
    matrix_i = multimodel.get_distance_matrix(prediction)

    # save feature vectors and names
    np.savez(f"{features_dir}/epoch_{e+1}.npz", names=np.array(x_species), features=prediction)
    print("Saving to:", Path(f"{features_dir}/epoch_{e+1}.npz").resolve())
    
    # calculate quartet scores (check how many quartets are in right split)
    scores = qt.get_qscores(matrix_i, x_species, scoring_batches)
    
    # make splitstree diagram from distance matrix
    #qt.matrix2nexus(matrix=matrix_i, taxa=x_species, nexusfile='filename.nex', plot_now=True)   
    qt.matrix2nexus(matrix=matrix_i, taxa=x_species, nexusfile=f'{out_dir}nexus/{e+1}.nex', plot_now=True)

        
    
    
    
    # track the evaluated measures 
    t.trackall(e+1, feature_vectors=prediction, loss=losses, score=scores)
    t.writeall()

    # you can immediately plot the results if you want 
    #t.plotall(e+1, 'test', plot_live=True)
    
    
    
  