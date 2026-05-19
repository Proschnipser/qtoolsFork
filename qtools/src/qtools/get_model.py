# -*- coding: utf-8 -*-
from qtools.encoding import onehot_encoding
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Conv1D, AveragePooling1D, Flatten, Dense, Dropout

   
def CNN_onehot_model(seq_len, output_dims=30):
    model = Sequential(
        [Conv1D(filters=30, kernel_size=10, strides=1,
                activation='relu', input_shape=(seq_len,4)),
           AveragePooling1D(pool_size=5),
            Dropout(0.2),
         Flatten(),
         Dense(360, activation='relu'),
         Dense(output_dims, activation='relu')
         ])
    return model

def fully_connected_model(vector_length,output_dims=5):
    model = Sequential([
        Dense(vector_length, activation='relu', input_shape=(vector_length,)),
        Dense(output_dims)#, activation='linear')
        ])
    return model


class CNN_ONEHOT:
    def __init__(self, seq_length, output_dims=30):
        self.model = CNN_onehot_model(seq_length, output_dims=output_dims)
        self.encoding_function = onehot_encoding        
        self.seq_len = seq_length
        self.input_dims = (1, seq_length, 4)
        self.output_dims = output_dims

class fully_connected:
    def __init__(self, vector_length, output_dims=5):
        self.model = fully_connected_model(vector_length, output_dims=output_dims)
        self.vector_length = vector_length
        self.input_dims =(vector_length,)
        self.output_dims = output_dims




