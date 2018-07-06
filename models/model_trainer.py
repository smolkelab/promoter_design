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

DNA = ['A','C','G','T']

def one_hot_encode(sequences):
  dna_dict = {}
  for i,q in enumerate(DNA):
    dna_dict[q] = i
  def _encode_one(seq, dna_dict):
    ans = np.zeros(shape=(len(DNA), len(seq)), dtype = 'int')
    for i,q in enumerate(seq):
      ans[i, dna_dict[q]] = 1
  ans = [_encode_one(q, dna_dict) for q in sequences]
  return(np.stack(ans, axis = 0))

def de_onehot(seq_arr):
  x = np.argmax(seq_arr, axis = 1)
  seqs = []
  for i in x:
    seq = ''.join([DNA[q] for q in i])
    seqs.append(seq)
  return(seqs)

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

  if not all([os.path.isfile(q + '.npy') for q in [train_file, valid_file, test_file]]):

    pre_text, post_text = (config.get('Params', 'pad_left'), config.get('Params', 'pad_right'))

    alldat = pandas.read_csv(filename, sep = ',', header=None)
    X = np.array([ pre_text + q + post_text for q in np.array(alldat['Seq']) ])
	output_names = [q for q in list(alldat) if q != 'Seq']
	output_names.sort()
    y = np.array([alldat[q] for q in output_names])
	num_outputs = y.shape[1]
	
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
	test_sequences = de_onehot(X_test)
    y_train_val = np.load(train_file + '_y.npy')
    y_valid_val = np.load(valid_file + '_y.npy')
    y_test_val = np.load(test_file + '_y.npy')
    print('Datasets loaded.')

  dat_to_use_all = [[X_train, y_train_val],[X_valid, y_valid_val], [X_test, y_test_val], X_train.shape[2]]

  #######################################################################################
  #Train the model
  #######################################################################################
  do_train = config.get('Mode', 'do_train').strip() == 'True'
  filename_sim = os.path.expanduser(config.get('Files', 'model_output'))

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
  output = {'Seqs': test_sequences}
  for i,q in enumerate(output_names):
    output[i] = y_test_val[i,:].squeeze()
	output['Pred_' + i] = preds[i,:].squeeze()

  pandas.DataFrame(output).to_csv(os.path.expanduser(config.get('Files','preds')))

if __name__ == '__main__':
  config = ConfigParser.RawConfigParser()
  config.read(sys.argv[1])
  sim, preds = main_method(config, True)