from keras.models import Sequential
from keras.callbacks import EarlyStopping
from keras.layers.core import Activation, Dense, Dropout, Flatten, Reshape
from keras.layers.convolutional import Convolution1D
from keras.layers.pooling import MaxPooling1D
from keras.regularizers import l2
from keras.optimizers import Adam
from keras.layers.normalization import BatchNormalization
from keras.layers import InputLayer
from keras import backend as K
from keras.initializers import RandomNormal
import numpy as np

import copy
from random import shuffle

SHIFT = 8
UNITS = 128
REG = 1e-4

def huber_loss(y_true, y_pred):
  d = 0.1
  x = K.abs(y_true - y_pred)
  d_t = d*K.ones_like(x)
  quad = K.min(K.stack([x, d_t], axis = -1), axis = -1)
  return( 0.5*K.square(quad) + d*(x - quad) )

class shifting_batch_generator(object):
  def __init__(self, dataset_X, dataset_y, batch_size, shift):
    self.dataset_X = dataset_X
    self.dataset_y = dataset_y
    self.batch_size = batch_size
    self.shift = shift
    self.window_size = dataset_X.shape[1] - shift + 1
    tuple_gen = ((p,q) for p in range(self.dataset_X.shape[0]) for q in range(self.shift))
    self.tuples = [i for i in tuple_gen]
    self._regen_tuples()

  def _regen_tuples(self):
    self.curr_tuples = copy.copy(self.tuples)
    shuffle(self.curr_tuples)

  def _get_batch(self):
    X_out = []
    y_out = []
    for i in range(self.batch_size):
      if len(self.curr_tuples) == 0:
        self._regen_tuples()
      sample_id, offset = self.curr_tuples.pop()
      X_out.append(self.dataset_X[sample_id,offset:offset+self.window_size,:])
      y_out.append(self.dataset_y[sample_id])
    return((np.stack(X_out,0), np.stack(y_out,0)))

  def iter(self):
    while True:
      yield(self._get_batch())


def do_model(dat_to_use, num_outputs, train = True):

  print(dat_to_use[0][0].shape)
  print(dat_to_use[0][1].shape)
  shift = SHIFT
  batch_size = 128

  def add_conv(model):
    model.add(Convolution1D(filters=UNITS, kernel_size=8, activation=None,strides=1,padding='same', kernel_regularizer = l2(REG), bias_regularizer = l2(REG)))
    model.add(BatchNormalization())
    model.add(Activation('relu'))
    model.add(MaxPooling1D(pool_size=2))
    return(model)

  model = Sequential()
  input_shape = (dat_to_use[3] - shift + 1, 4)
  print(input_shape)
  model.add(InputLayer(batch_input_shape=(None,) + input_shape))
  for i in range(6):
    model = add_conv(model)
  model.add(Flatten())
  model.add(Dense(units=UNITS, kernel_regularizer = l2(REG), bias_regularizer = l2(REG)))
  model.add(BatchNormalization())
  model.add(Activation('relu'))
  model.add(Dense(units=UNITS, kernel_regularizer = l2(REG), bias_regularizer = l2(REG)))
  model.add(BatchNormalization())
  model.add(Activation('relu'))
  model.add(Dense(units=num_outputs))

  model.compile(optimizer = Adam(lr=1e-5,decay=0.), loss=huber_loss)
  print(model.summary())

  if train:
    earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=1)
    gen_train = shifting_batch_generator(dat_to_use[0][0], dat_to_use[0][1], batch_size, shift)
    gen_valid = shifting_batch_generator(dat_to_use[1][0], dat_to_use[1][1], batch_size, shift)
    model.fit_generator(generator =  gen_train.iter(), steps_per_epoch =  shift*dat_to_use[0][0].shape[0]/batch_size, # was 4!
      epochs = 100, callbacks = [earlystopper], validation_data = gen_valid.iter(),
      validation_steps = shift*dat_to_use[1][0].shape[0]/batch_size, verbose = 2)
  return(model)
