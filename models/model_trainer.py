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
    ans = np.zeros(shape=(len(seq), len(DNA)), dtype = 'int')
    for i,q in enumerate(seq):
      ans[i, dna_dict[q]] = 1
    return(ans)
  ans = [_encode_one(q, dna_dict) for q in sequences]
  return(np.stack(ans, axis = 0))

def de_onehot(seq_arr):
  x = np.argmax(seq_arr, axis = 2)
  seqs = []
  for i in x:
    seq = ''.join([DNA[q] for q in i])
    seqs.append(seq)
  return(seqs)

def main_method(config, valid_split = None):

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

  if not all([os.path.isfile(q + '.npy') for q in [train_file, valid_file, test_file]]) and valid_split == None:

    pre_text, post_text = (config.get('Params', 'pad_left'), config.get('Params', 'pad_right'))
    alldat = pandas.read_csv(filename, sep = ',')
    X = np.array([ pre_text + q + post_text for q in np.array(alldat['Seq']) ])
    output_names = [q for q in list(alldat) if q != 'Seq']
    output_names.sort()
    y = np.array([alldat[q] for q in output_names])
    y = np.swapaxes(y, 0, 1)

    test_set_size = np.floor(len(X)*test_frac).astype(int)
    validation_set_size = np.floor(len(X)*val_frac).astype(int)

    if random_seed > 0:
      random.seed(random_seed); np.random.seed(random_seed)
    train_sequences, test_sequences, y_train_val, y_test_val = train_test_split(X, y, test_size=test_set_size)
    if random_seed > 0:
      random.seed(random_seed); np.random.seed(random_seed)
    train_sequences, valid_sequences, y_train_val, y_valid_val = train_test_split(train_sequences, y_train_val, test_size=validation_set_size)

    X_train = one_hot_encode(train_sequences)
    X_valid = one_hot_encode(valid_sequences)
    X_test  = one_hot_encode(test_sequences)
    np.save(train_file, X_train)
    np.save(valid_file, X_valid)
    np.save(test_file, X_test)
    np.save(train_file + '_y', y_train_val)
    np.save(valid_file + '_y', y_valid_val)
    np.save(test_file + '_y', y_test_val)
    # save output_names if reloading
    with open(test_file + '_names.txt', 'w') as fn:
      for n in output_names:
        fn.write(n + '\n')

  else:
    print('Loading saved datasets...')
    X_test = np.load(test_file + '.npy')
    y_test_val = np.load(test_file + '_y.npy')
    test_sequences = de_onehot(X_test)
    output_names = []
    with open(test_file + '_names.txt', 'r') as fn:
      for l in fn:
        output_names.append(l.strip())
    print('Output names: ' + str(output_names))

    if valid_split == None:
      X_train = np.load(train_file + '.npy')
      X_valid = np.load(valid_file + '.npy')
      y_train_val = np.load(train_file + '_y.npy')
      y_valid_val = np.load(valid_file + '_y.npy')

    # Assemble the training dataset from cross-validation splits
    else:
      num_valid = int(config.get('Params','num_valid'))
      X_all_valid = [np.load(valid_file + '_' + str(q) + '.npy') for q in range(num_valid)]
      y_all_valid = [np.load(valid_file + '_' + str(q) + '_y.npy') for q in range(num_valid)]
      X_valid = X_all_valid.pop(valid_split)
      y_valid_val = y_all_valid.pop(valid_split)
      X_train = np.concatenate(X_all_valid, axis = 0)
      y_train_val = np.concatenate(y_all_valid, axis = 0)

    print('Datasets loaded.')

  dat_to_use_all = [[X_train, y_train_val],[X_valid, y_valid_val], [X_test, y_test_val], X_train.shape[1]]

  #######################################################################################
  #Train the model
  #######################################################################################
  do_train = config.get('Mode', 'do_train').strip() == 'True'
  filename_sim = os.path.expanduser(config.get('Files', 'model_output'))
  num_outputs = y_test_val.shape[1]
  sim = do_model.do_model(dat_to_use_all, num_outputs, train = do_train)

  if do_train:
    model_json = sim.to_json()
    with open(filename_sim+'.json', "w") as json_file:
      json_file.write(model_json)
    if valid_split == None:
      sim.save_weights(filename_sim+'.h5')
    else:
      sim.save_weights(filename_sim + '_' + str(valid_split) + '.h5')

  else:
    if valid_split == None:
      sim.load_weights(filename_sim+'.h5')
    else:
      sim.load_weights(filename_sim + '_' + str(valid_split) + '.h5')

    print('Weights loaded.')
  #######################################################################################
  #Test the model
  #######################################################################################

  preds = sim.predict(X_test[:,0:X_test.shape[1]-do_model.SHIFT+1,:])
  output = {'Seqs': test_sequences}
  for i,q in enumerate(output_names):
    output[q] = y_test_val[:,i].squeeze()
    output['Pred_' + q] = preds[:,i].squeeze()

  output = pandas.DataFrame(output)
  if valid_split == None:
    output.to_csv(os.path.expanduser(config.get('Files','preds')), index = False)
  else:
    fn_output = os.path.expanduser(config.get('Files','preds')).split('.')
    fn_output = fn_output[0] + '_' + str(valid_split) + '.' + fn_output[1]
    output.to_csv(fn_output, index = False)

  print('Predictions written')

if __name__ == '__main__':
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.read(sys.argv[1])
  if len(sys.argv) > 2:
    valid_split = int(sys.argv[2])
  else:
    valid_split = None
  main_method(config, valid_split)
