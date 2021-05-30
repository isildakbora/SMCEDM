import math

import numpy as np
import pandas as pd

import ROOT, uproot

from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn import feature_selection
from sklearn.model_selection import StratifiedKFold
from sklearn.utils import class_weight
from sklearn import decomposition

# import necessary keras modules/functions
import keras
from keras.models import Sequential
from keras.models import Model
from keras.layers import Input
from keras.layers import Dense
from keras.layers import Dropout
from keras.layers import BatchNormalization
from keras.callbacks import ReduceLROnPlateau

from keras import backend as K

import tensorflow as tf

import itertools as it

import time

import pickle

def upsample(sample, scale):
    '''
    This function takes the sample and concatantes the same sample n times.
    '''
    if scale > 0:
        sample = np.vstack([sample]*scale)
    else:
        sample = sample
    return sample

# set the configuration
config = tf.ConfigProto( device_count = {'GPU': 1 , 'CPU': 36} ) 
sess = tf.Session(config=config) 
keras.backend.set_session(sess)

# Convert root files to pandas dataframes
input_branches = ['br_njets', 'br_nbjets', 'br_scalar_ht',
       'br_jet_pt_1', 'br_jet_pt_2', 'br_jet_pt_3', 'br_jet_pt_4', 'br_met',
       'br_met_phi', 'br_sphericity', 'br_aplanarity', 'br_fox_wolfram_1',
       'br_fox_wolfram_2', 'br_fox_wolfram_3', 'br_fox_wolfram_4', 'br_w_pt', 'br_w_eta', 'br_w_phi']

df_signal = uproot.open("/mnt/harddisk4/scratch/ttbar_no_dtG_flat.root")["outtree"].pandas.df(input_branches)
df_dyjets = uproot.open("/mnt/harddisk4/scratch/dyjets_flat_30M.root")["outtree"].pandas.df(input_branches)
df_wjets  = uproot.open("/mnt/harddisk4/scratch/wjets_flat_60M.root")["outtree"].pandas.df(input_branches)
df_single_top  = uproot.open("/mnt/harddisk4/scratch/single_top_merged_flat.root")["outtree"].pandas.df(input_branches)


x_signal     = df_signal.values
y_signal     = np.full(len(x_signal) , 1)

x_dyjets     = df_dyjets.values
x_wjets      = df_wjets.values
x_single_top = df_single_top.values

x_dyjets     = upsample(x_dyjets, int(1/3*len(x_signal)/len(x_dyjets)))
x_wjets      = upsample(x_wjets, int(1/3*len(x_signal)/len(x_wjets)))
x_single_top = x_single_top[np.random.randint(len(x_single_top), size=len(x_wjets))]

print(len(x_signal))
x_bkg    = np.concatenate((x_dyjets, x_wjets, x_single_top), axis=0)
y_bkg    = np.full(len(x_bkg) , 0)

f1 = 0.5 # train_fraction
f2 = f1+0.5*(1-f1) # train+val_fraction
idx_sig = np.arange(len(x_signal))
np.random.shuffle(idx_sig)
idx_sig = np.split(idx_sig, [int(f1*len(idx_sig)), int(f2*len(idx_sig))])
idx_sig_train, idx_sig_val, idx_sig_test = idx_sig[0], idx_sig[1], idx_sig[2]
x_sig_train = x_signal[idx_sig_train]
x_sig_val   = x_signal[idx_sig_val]
x_sig_test  = x_signal[idx_sig_test]
y_sig_train = np.full(len(x_sig_train), 1)
y_sig_test  = np.full(len(x_sig_test), 1)
y_sig_val   = np.full(len(x_sig_val), 1)

idx_dyjets = np.arange(len(x_dyjets))
np.random.shuffle(idx_dyjets)
idx_dyjets = np.split(idx_dyjets, [int(f1*len(idx_dyjets)), int(f2*len(idx_dyjets))])
idx_dyjets_train, idx_dyjets_val, idx_dyjets_test = idx_dyjets[0], idx_dyjets[1], idx_dyjets[2]
x_dyjets_train = x_dyjets[idx_dyjets_train]
x_dyjets_val   = x_dyjets[idx_dyjets_val]
x_dyjets_test  = x_dyjets[idx_dyjets_test]
y_dyjets_train = np.full(len(x_dyjets_train), 0)
y_dyjets_val   = np.full(len(x_dyjets_val), 0)
y_dyjets_test  = np.full(len(x_dyjets_test), 0)

idx_wjets = np.arange(len(x_wjets))
np.random.shuffle(idx_wjets)
idx_wjets = np.split(idx_wjets, [int(f1*len(idx_wjets)), int(f2*len(idx_wjets))])
idx_wjets_train, idx_wjets_val, idx_wjets_test = idx_wjets[0], idx_wjets[1], idx_wjets[2]
x_wjets_train = x_wjets[idx_wjets_train]
x_wjets_val   = x_wjets[idx_wjets_val]
x_wjets_test  = x_wjets[idx_wjets_test]
y_wjets_train = np.full(len(x_wjets_train), 0)
y_wjets_val   = np.full(len(x_wjets_val), 0)
y_wjets_test  = np.full(len(x_wjets_test), 0)

idx_single_top = np.arange(len(x_single_top))
np.random.shuffle(idx_single_top)
idx_single_top = np.split(idx_single_top, [int(f1*len(idx_single_top)), int(f2*len(idx_single_top))])
idx_single_top_train, idx_single_top_val, idx_single_top_test = idx_single_top[0], idx_single_top[1], idx_single_top[2]
x_single_top_train = x_single_top[idx_single_top_train]
x_single_top_val   = x_single_top[idx_single_top_val]
x_single_top_test  = x_single_top[idx_single_top_test]
y_single_top_train = np.full(len(x_single_top_train), 0)
y_single_top_val   = np.full(len(x_single_top_val), 0)
y_single_top_test  = np.full(len(x_single_top_test), 0)

print("number of signal events for the training :%d" % (len(x_sig_train)))
print("number of dyjets events for the training:%d" % (len(x_dyjets_train)))
print("number of wjets events for the training:%d" % (len(x_wjets_train)))
print("number of single_top events for the training:%d" % (len(x_single_top_train)))

from sklearn.utils import shuffle

x_train= np.concatenate((x_sig_train, x_dyjets_train, x_wjets_train, x_single_top_train), axis=0)
y_train = np.concatenate((y_sig_train, y_dyjets_train, y_wjets_train, y_single_top_train), axis=0)
x_train, y_train = shuffle(x_train, y_train, random_state=0)

x_val = np.concatenate((x_sig_val, x_dyjets_val, x_wjets_val, x_single_top_val), axis=0)
y_val = np.concatenate((y_sig_val, y_dyjets_val, y_wjets_val, y_single_top_val), axis=0)
x_val, y_val = shuffle(x_val, y_val, random_state=0)

x_test = np.concatenate((x_sig_test, x_dyjets_test, x_wjets_test, x_single_top_test), axis=0)
y_test = np.concatenate((y_sig_test, y_dyjets_test, y_wjets_test, y_single_top_test), axis=0)

scaler  = preprocessing.StandardScaler().fit(np.concatenate((x_train, x_val, x_test), axis=0))
pickle.dump(scaler, open('scaler.pkl', 'wb'))
x_train = scaler.transform(x_train)
x_val   = scaler.transform(x_val)

np.savez('test_data.npz', x_sig_test = x_sig_test, x_dyjets_test = x_dyjets_test, x_wjets_test = x_wjets_test,
         x_single_top_test = x_single_top_test)

print("training for", len(x_train), "events;", "testing for",  len(x_val), "events")

# define the model with variable hyperparameters
def dnn_model(x_tr, y_tr, x_ts, y_ts, callbacks, params):
    
    model = Sequential()
    
    model.add(Dense(params['first_hidden_layer'], input_dim=x_tr.shape[1],
                            activation=params['activation1'],
                            kernel_initializer=params['kernel_initializer']))
    
    #model.add(Dropout(params['dropout']))
    
    model.add(Dense(params['second_hidden_layer'], 
                            activation=params['activation2'],
                            use_bias=True,
                            kernel_regularizer=params['weight_regulizer']))

    model.add(BatchNormalization())
    
    model.add(Dense(params['third_hidden_layer'], 
                            activation=params['activation3'],
                            use_bias=True,
                            kernel_regularizer=params['weight_regulizer']))
    
    model.add(Dense(params['fourth_hidden_layer'], 
                            activation=params['activation4'],
                            use_bias=True,
                            kernel_regularizer=params['weight_regulizer']))
    
    model.add(Dense(1, activation=params['last_activation']))
    
    opt  = keras.optimizers.Adam(lr=params['learn_rate'])
    
    model.compile(loss=params['losses'],
                            optimizer=opt,
                            metrics=['binary_accuracy'])
    
    x_tr = np.asarray(x_tr).astype('float32')
    y_tr = np.asarray(y_tr).astype('float32').reshape((-1,1))
    x_ts = np.asarray(x_ts).astype('float32')
    y_ts = np.asarray(y_ts).astype('float32').reshape((-1,1))
    
    # Fit the Keras model on the dataset
    steps_per_epoch = int(np.ceil(x_tr.shape[0] / params['batch_size'])) - 1
    
    t_start = time.time()
    history = model.fit(x_tr, y_tr, 
                        validation_data=[x_ts, y_ts],
                        batch_size=params['batch_size'],
                        epochs=params['epochs'],
                        callbacks=[*callbacks],
                        verbose=1)
    t_stop = time.time()
    
    delta_t = t_stop - t_start
        
    return history, delta_t 

# construct the hyperparameter grid to be scanned

param_grid = {
     'activation1':["relu", "tanh"],
     'activation2':["relu", "tanh"],
     'activation3':["relu", "tanh"],
     'activation4':["relu", "tanh"],
     'first_hidden_layer': [64, 128],
     'second_hidden_layer': [64, 128],
     'third_hidden_layer': [64, 64],
     'fourth_hidden_layer': [32, 128],
     'batch_size': [2**9],
     'epochs': [50],
     'learn_rate': [0.00146],
     #'dropout': [0, 0.1],
     'weight_regulizer':[None],
     'emb_output_dims': [None],
     'optimizer': [ 'Adam'],   
     'losses': ['binary_crossentropy','logcosh'],
     'last_activation': ["sigmoid", "tanh"],
     'kernel_initializer':["glorot_normal"]}

# construct dicts for all possible hyperparameter permutations
keys, values = zip(*param_grid.items())
permutations_dicts = [dict(zip(keys, v)) for v in it.product(*values)]     

# train models with all possible hyperparameter sets and save models and accuracies

print("training for all possible " + str(len(permutations_dicts)) + " hyperparameter sets.")
history_arr = []

#callbacks = [ReduceLROnPlateau(monitor='val_binary_accuracy', factor=0.95, 
#                                    patience=5, verbose=1, min_delta=1e-4, mode='auto', min_lr=1.0e-5)]

callbacks = []

for i, param_set in enumerate(permutations_dicts):
    print('training for hyperparameter set ' + str(i) + '.')
    history, delta_t = dnn_model(x_train, y_train, x_val, y_val, callbacks, param_set) 
    history_arr.append([param_set, history.history['val_loss'], history.history['val_binary_accuracy'], history.history['loss'],
                        history.history['binary_accuracy'], delta_t/param_set['epochs']])
    
    tf.keras.backend.clear_session()

with open('history.log', 'wb') as fp:
        pickle.dump(history_arr, fp)