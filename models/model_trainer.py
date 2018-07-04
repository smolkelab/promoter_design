import sys
import os
import numpy as np
import copy
import pandas
import random
from sklearn.cross_validation import train_test_split
import ConfigParser
import do_model
import scipy.stats as ss

from sklearn.preprocessing import LabelEncoder, OneHotEncoder

#def one_hot_encode(sequences):
#    sequence_length = len(sequences[0])
#    integer_type = np.int8 if sys.version_info[
#        0] == 2 else np.int32  # depends on Python version
#    integer_array = LabelEncoder().fit(np.array(('ACGTN',)).view(integer_type)).transform(
#        sequences.view(integer_type)).reshape(len(sequences), sequence_length)
#    one_hot_encoding = OneHotEncoder(
#        sparse=False, n_values=5).fit_transform(integer_array)
#    return one_hot_encoding.reshape(len(sequences), 1, sequence_length, 5).swapaxes(2, 3)[:, :, [0, 1, 2, 4], :]
def main_method(config):

  #######################################################################################
  #Preprocessing
  #######################################################################################
  filename = config.get('Files', 'means')
  val_frac = float(config.get('Params', 'val_frac'))
  test_frac = float(config.get('Params', 'test_frac'))
  random_seed = int(config.get('Params', 'random_seed'))
  train_file = config.get('Files', 'train')
  valid_file = config.get('Files', 'valid')
  test_file = config.get('Files', 'test')

  train_file, valid_file, test_file = [q + '_' + str(random_seed) for q in [train_file, valid_file, test_file]]

  if not all([os.path.isfile(q + '.npy') for q in [train_file, valid_file, test_file]]):

    pre_text, post_text = (config.get('Params', 'pad_left'), config.get('Params', 'pad_right'))

    alldat = pandas.read_csv(filename, sep = ',', header=None)
    y_a = np.array(alldat[2])
    y_b = np.array(alldat[3])
    y = np.stack([y_a, y_b], axis = 1)
	num_outputs = int(config.get('Params','num_outputs'))
	assert(num_outputs in [1,2])
	if num_outputs == 1:
	  y = np.apply_along_axis(np.mean, 1, y)
	  y_min = float(config.get('Params','y_min'))
	  y_max = float(config.get('Params','y_max'))
	  useable = np.logical_and(y > y_min, y < y_max)
	else:
	  y_a_min = float(config.get('Params','y_a_min'))
	  y_a_max = float(config.get('Params','y_a_max'))
	  y_b_min = float(config.get('Params','y_b_min'))
	  y_b_max = float(config.get('Params','y_b_max'))
	  useable = np.logical_and(np.logical_and(y_a > y_a_min, y_a < y_a_max), np.logical_and(y_b > y_b_min, y_b < y_b_max))
    X = np.array([ pre_text + q + post_text for q in np.array(alldat['Seq']) ])
    X = X[useable]
    y = y[useable,...]

    test_set_size = np.floor(len(X)*test_frac).astype(int)
    validation_set_size = np.floor(len(X)*val_frac).astype(int)

    if random_seed > 0:
      random.seed(random_seed); np.random.seed(random_seed)
    train_sequences, test_sequences, y_train_val, y_test_val = train_test_split(X, y, test_size=test_set_size)
    if random_seed > 0:
      random.seed(random_seed); np.random.seed(random_seed)

    train_sequences, valid_sequences, y_train_val, y_valid_val = train_test_split(train_sequences, y_train_val, test_size=validation_set_size)

    X_train = one_hot_encode(train_sequences).squeeze()  # get rid of the unwanted extra dimension
    X_valid = one_hot_encode(valid_sequences).squeeze()
    X_test  = one_hot_encode(test_sequences).squeeze()
    np.save(train_file, X_train)
    np.save(valid_file, X_valid)
    np.save(test_file, X_test)
    np.save(train_file + '_y', y_train_val)
    np.save(valid_file + '_y', y_valid_val)
    np.save(test_file + '_y', y_test_val)

  else:
    print('Loading saved datasets...')
    X_train = np.load(train_file + '.npy')
    X_valid = np.load(valid_file + '.npy')
    X_test = np.load(test_file + '.npy')
    y_train_val = np.load(train_file + '_y.npy')
    y_valid_val = np.load(valid_file + '_y.npy')
    y_test_val = np.load(test_file + '_y.npy')
    print('Datasets loaded.')

  dat_to_use_all = [[X_train, y_train_val],[X_valid, y_valid_val], [X_test, y_test_val], X_train.shape[2]]

  #######################################################################################
  #Train the model
  #######################################################################################
  do_train = config.get('Mode', 'do_train').strip() == 'True'
  filename_sim = config.get('Files', 'model_output')

  sim = do_model.do_model(dat_to_use_all, num_outputs, train = do_train)

  if do_train:
    model_json = sim.to_json()
    with open(filename_sim+'.json', "w") as json_file:
      json_file.write(model_json)
    sim.save_weights(filename_sim+'.h5')

  else:
    sim.load_weights(filename_sim+'.h5')
    print('Weights loaded.')
  #######################################################################################
  #Test the model
  #######################################################################################

  preds = sim.predict(X_test[...,0:X_test.shape[2]-do_model.SHIFT+1]).squeeze()
  if num_outputs == 1:
    output = {'Means': y_test_val.squeeze(), 'Preds':preds}
  else:
    output = {'Means_A': y_test_val[:,0].squeeze(),
              'Means_B': y_test_val[:,1].squeeze(),
              'Preds_A': preds[:,0].squeeze(),
              'Preds_B': preds[:,1].squeeze()}    

  pandas.DataFrame(output).to_csv(config.get('Files','preds'))

if __name__ == '__main__':
  config = ConfigParser.RawConfigParser()
  config.read(sys.argv[1])
  sim, preds = main_method(config, True)