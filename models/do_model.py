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

SHIFT = 8

def huber_loss(y_true, y_pred):
  d = 0.15
  x = K.abs(y_true - y_pred)
  d_t = d*K.ones_like(x)
  quad = K.min(K.stack([x, d_t], axis = -1), axis = -1)
  return( 0.5*K.square(quad) + d*(x - quad) )

def shifting_batch_generator(dataset_X, dataset_y, batch_size, shift):
  window_size = dataset_X.shape[1] - shift + 1
  while True:
    X_out = []
    y_out = []
    for i in range(batch_size):
      sample_id = np.random.random_integers(0,dataset_X.shape[0]-1)
      offset = np.random.random_integers(0, shift-1)
      X_out.append(dataset_X[sample_id,offset:offset+window_size,:])
      y_out.append(dataset_y[sample_id])
    yield((np.stack(X_out,0), np.stack(y_out,0)))

def do_model(dat_to_use, num_outputs, train = True):
  print(dat_to_use[0][0].shape)
  print(dat_to_use[0][1].shape)
  shift = SHIFT
  batch_size = 128

  def add_conv(model):
    model.add(Convolution1D(filters=128, kernel_size=8, activation=None,strides=1,padding='same', kernel_regularizer = l2(1e-4), bias_regularizer = l2(1e-4)))
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
  model.add(Dense(units=128, kernel_regularizer = l2(1e-4), bias_regularizer = l2(1e-4)))
  model.add(BatchNormalization())
  model.add(Activation('relu'))
  model.add(Dense(units=128, kernel_regularizer = l2(1e-4), bias_regularizer = l2(1e-4)))
  model.add(BatchNormalization())
  model.add(Activation('relu'))
  model.add(Dense(units=num_outputs))

  model.compile(optimizer = Adam(lr=1e-5,decay=0.), loss=huber_loss)
  print(model.summary())

  if train:
    earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=1)
    model.fit_generator(generator =  shifting_batch_generator(dat_to_use[0][0], dat_to_use[0][1], batch_size, shift), steps_per_epoch =  shift*dat_to_use[0][0].shape[0]/batch_size, # was 4!
      epochs = 100, callbacks = [earlystopper], validation_data = shifting_batch_generator(dat_to_use[1][0], dat_to_use[1][1], batch_size, shift), 
      validation_steps = shift*dat_to_use[1][0].shape[0]/batch_size, verbose = 2)
  return(model)
