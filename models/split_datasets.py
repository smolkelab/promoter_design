# Divide raw data into partitions; maintain backwards compatibility with prior training
# test data will be saved separately, validation partition 0 will be original validation data

import sys
import os
import pandas
import numpy as np
import ConfigParser
import random
from sklearn.cross_validation import train_test_split
from model_trainer import one_hot_encode

def main_method(config):

  #######################################################################################
  #Preprocessing
  #######################################################################################
  filename = os.path.expanduser(config.get('Files', 'means'))
  val_frac = float(config.get('Params', 'val_frac'))
  test_frac = float(config.get('Params', 'test_frac'))
  random_seed = int(config.get('Params', 'random_seed'))
  train_file = os.path.expanduser(config.get('Files', 'train'))
  valid_file = os.path.expanduser(config.get('Files', 'valid'))
  test_file = os.path.expanduser(config.get('Files', 'test'))

  train_file, valid_file, test_file = [q + '_' + str(random_seed) for q in [train_file, valid_file, test_file]]

  pre_text, post_text = (config.get('Params', 'pad_left'), config.get('Params', 'pad_right'))
  alldat = pandas.read_csv(filename, sep = ',')
  X = np.array([ pre_text + q + post_text for q in np.array(alldat['Seq']) ])
  output_names = [q for q in list(alldat) if q != 'Seq']
  output_names.sort()
  y = np.array([alldat[q] for q in output_names])
  y = np.swapaxes(y, 0, 1)

  test_set_size = np.floor(len(X)*test_frac).astype(int)
  validation_set_size = np.floor(len(X)*val_frac).astype(int)

  # Select the test set
  if random_seed > 0:
    random.seed(random_seed); np.random.seed(random_seed)
  train_sequences, test_sequences, y_train_val, y_test_val = train_test_split(X, y, test_size=test_set_size)
  X_test = one_hot_encode(test_sequences)
  np.save(test_file, X_test)
  np.save(test_file + '_y', y_test_val)

  # Select validation sets
  i = 0; going = True
  while going: # will break when we run out of sequences
    if random_seed > 0:
      random.seed(random_seed); np.random.seed(random_seed)
    try:
      train_sequences, valid_sequences, y_train_val, y_valid_val = train_test_split(train_sequences, y_train_val, test_size=validation_set_size)
      X_valid = one_hot_encode(valid_sequences)
      np.save(valid_file + '_' + str(i), X_valid)
      np.save(valid_file + '_' + str(i) + '_y', y_valid_val)
    except ValueError: # save whatever's left in 'train_sequences' and break
      X_valid = one_hot_encode(train_sequences)
      np.save(valid_file + '_' + str(i), X_valid)
      np.save(valid_file + '_' + str(i) + '_y', y_train_val)
      going = False
    i += 1

  # save output_names if reloading
  with open(test_file + '_names.txt', 'w') as fn:
    for n in output_names:
      fn.write(n + '\n')

if __name__ == '__main__':
  config = ConfigParser.RawConfigParser()
  config.read(sys.argv[1])
  main_method(config)
