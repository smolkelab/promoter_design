from keras.models import Sequential, Model
from keras.callbacks import Callback, EarlyStopping, TensorBoard
from keras.layers.core import (
    Activation, Dense, Dropout, Flatten,
    Permute, Reshape)
from keras.layers.convolutional import Convolution1D
from keras.layers.pooling import MaxPooling1D, AveragePooling1D
from keras.layers.recurrent import GRU, SimpleRNN, LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.regularizers import l1, l2
from keras.layers.advanced_activations import LeakyReLU, PReLU
from keras.optimizers import Adam
from keras.layers.normalization import BatchNormalization
from keras.layers import Input, concatenate, Lambda, LocallyConnected1D
from keras.layers import GlobalAveragePooling1D
from keras import backend as K
from keras.initializers import RandomNormal
#import imp
#imp.load_source('Parametric_Blur','Parametric_Blur.py')
#from Parametric_Blur import Parametric_Blur

import numpy as np

SHIFT = 8

def primer_model():
  def add_conv(model):
    model.add(Convolution1D(filters=128, kernel_size=8, activation=None,strides=1,padding='same',
               kernel_regularizer = l2(5e-5), bias_regularizer = l2(5e-5)))
    model.add(BatchNormalization())
    model.add(Activation('relu'))
    model.add(MaxPooling1D(pool_size=2))
    return(model)

  model = Sequential()
  input_shape = (4, 362 - 8 + 1)
  model.add(Permute((2,1), input_shape = input_shape))
  model = add_conv(model)
  model = add_conv(model)
  model = add_conv(model)
  model.add(Flatten())
  model.add(Dense(units=64, kernel_regularizer = l2(5e-5), bias_regularizer = l2(5e-5)))
  model.add(BatchNormalization())
  model.add(Activation('relu'))
  model.add(Dense(units=128, kernel_regularizer = l2(5e-5), bias_regularizer = l2(5e-5)))
  model.add(BatchNormalization())
  model.add(Activation('relu'))
  model.add(Dense(units=1))
  return(model)

def huber_loss(y_true, y_pred):
  d = 0.15
  x = K.abs(y_true - y_pred)
  d_t = d*K.ones_like(x)
  quad = K.min(K.stack([x, d_t], axis = -1), axis = -1)
  return( 0.5*K.square(quad) + d*(x - quad) )

def shifting_batch_generator(dataset_X, dataset_y, batch_size, shift):
  window_size = dataset_X.shape[2] - shift + 1
  while True:
    X_out = []
    y_out = []
    for i in range(batch_size):
      sample_id = np.random.random_integers(0,dataset_X.shape[0]-1)
      offset = np.random.random_integers(0, shift-1)
      X_out.append(dataset_X[sample_id,:,offset:offset+window_size])
      y_out.append(dataset_y[sample_id])

    yield((np.stack(X_out,0), np.stack(y_out,0)))

def explicit_augmentor(X, y, shift):
  window_size = X.shape[2] - shift + 1
  ans = []
  for i in range(shift):
    tmp = X[...,i:(i+window_size)]
    ans.append((np.copy(tmp), np.copy(y)))
  X_out = np.concatenate([q[0] for q in ans], axis = 0)
  y_out = np.concatenate([q[1] for q in ans], axis = 0)
  return(X_out, y_out)

def do_model(dat_to_use, train = True, start_weights = None, **kwargs):

  shift = SHIFT
  # raises MemoryError on current machine :(
  #tr_X, tr_y = explicit_augmentor(dat_to_use[0][0], dat_to_use[0][1], shift)
  #va_X, va_y = explicit_augmentor(dat_to_use[1][0], dat_to_use[1][1], shift)

  batch_size = 128

  def add_conv(model, trainable = True):
    model.add(Convolution1D(filters=128, kernel_size=8, activation=None,strides=1,padding='same', kernel_regularizer = l2(1e-4), bias_regularizer = l2(1e-4), trainable = trainable)) # prev. had l1 of 1e-5
    model.add(BatchNormalization())
    model.add(Activation('relu'))
    #model.add(Parametric_Blur(width=7, sd_init=1.))
    model.add(MaxPooling1D(pool_size=2))
    return(model)

  model = Sequential()
  input_shape = (4, dat_to_use[3] - shift + 1)
  model.add(Permute((2,1), input_shape = input_shape))
  model = add_conv(model, trainable = True)
  model = add_conv(model, trainable = True)
  model = add_conv(model, trainable = True)
  model = add_conv(model)
  model = add_conv(model)
  model = add_conv(model)
  model.add(Flatten())
  model.add(Dense(units=128, kernel_regularizer = l2(1e-4), bias_regularizer = l2(1e-4)))
  model.add(BatchNormalization())
  model.add(Activation('relu'))
  model.add(Dense(units=128, kernel_regularizer = l2(1e-4), bias_regularizer = l2(1e-4)))
  model.add(BatchNormalization())
  model.add(Activation('relu'))
  model.add(Dense(units=2))

  #model.compile(optimizer=Adam(lr = 5e-5, decay = 0.), loss='logcosh')
  #model.compile(optimizer=Adam(lr = 5e-6, decay = 0.), loss=huber_loss)

  if start_weights != None:
    model.load_weights(start_weights) # start_weights is a filename

  if 'pretrained_weights' in kwargs.keys():
    pw, num_weights = kwargs.pop('pretrained_weights')
    if pw is not None and num_weights is not None:
      p_mod = primer_model() # hacky!
      p_mod.load_weights(pw + '.h5')
      print('Primer model loaded')
      print(p_mod.summary())

      print('p_mod 0 sum : ' + str(np.sum(p_mod.get_weights()[0])))
      print('model 0 pre-set sum: ' + str(np.sum(model.get_weights()[0])))

      for i in range(num_weights):
        model.layers[i].set_weights(p_mod.layers[i].get_weights())

      print('model 0 post-set sum: ' + str(np.sum(model.get_weights()[0])))

  model.compile(optimizer = Adam(lr=1e-5,decay=0.), loss=huber_loss)
  print('model 0 compiled sum: ' + str(np.sum(model.get_weights()[0])))
  print(model.summary())


  if train:
    earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=1)
    #model.fit(x = tr_X, y = tr_y, batch_size = batch_size, epochs = 100, verbose = 2,
    #          callbacks = [earlystopper], validation_data = (va_X, va_y))
    #model.fit(x = dat_to_use[0][0], y = dat_to_use[0][1], batch_size = batch_size, epochs = 100, verbose = 2,
    #          callbacks = [earlystopper], validation_data = (dat_to_use[1][0], dat_to_use[1][1]))
    model.fit_generator(generator =  shifting_batch_generator(dat_to_use[0][0], dat_to_use[0][1], batch_size, shift), steps_per_epoch =  8*dat_to_use[0][0].shape[0]/batch_size, # was 4!
      epochs = 100, callbacks = [earlystopper], validation_data = shifting_batch_generator(dat_to_use[1][0], dat_to_use[1][1], batch_size, shift), 
      validation_steps = 8*dat_to_use[1][0].shape[0]/batch_size, verbose = 2) # was 4!

  return(model)
