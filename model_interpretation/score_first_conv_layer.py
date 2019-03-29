# given a file of weights, find weights for the first convolutional layer and corresponding BatchNorm.
# For each possible 8-mer, get the output from each filter after convolution and BatchNorm,
# but before ReLU.
# Structure: get a function that takes in the sequence - represented as a list of the indices corresponding
# to each base in the sequence - and outputs the score.
# Execute this function for each filter/8-mer combination. Save the result as a table.

import sys
import os
import numpy as np
import h5py
from functools import partial

DNA = ['A','C','G','T']
BASEFOUR_DICT = {'00':0, '01':1, '10':2, '11':3}
NUM_FILTERS = 128

CONV_KERNEL_KEYS = ['conv1d_1', 'conv1d_1', 'kernel:0']
CONV_BIAS_KEYS = ['conv1d_1', 'conv1d_1', 'bias:0']
BATCHNORM_BETA_KEYS = ['batch_normalization_1','batch_normalization_1','beta:0']
BATCHNORM_GAMMA_KEYS = ['batch_normalization_1','batch_normalization_1','gamma:0']
BATCHNORM_MEAN_KEYS = ['batch_normalization_1','batch_normalization_1','moving_mean:0']
BATCHNORM_VAR_KEYS = ['batch_normalization_1','batch_normalization_1','moving_variance:0']

def get_weights_helper(f, key_list):
  for q in key_list:
    f = f[q]
  return f.value

def score_function_skeleton(f, seq_idxes):
  wts_conv_kernel = get_weights_helper(f, CONV_KERNEL_KEYS)
  assert len(seq_idxes) == wts_conv_kernel.shape[0]
  wts_conv_bias = get_weights_helper(f, CONV_BIAS_KEYS)
  b = get_weights_helper(f, BATCHNORM_BETA_KEYS)
  g = get_weights_helper(f, BATCHNORM_GAMMA_KEYS)
  m = get_weights_helper(f, BATCHNORM_MEAN_KEYS)
  v = get_weights_helper(f, BATCHNORM_VAR_KEYS)

  # "fancy indexing" to get a 1D array of each (row_id, column_id) desired - cf. https://jakevdp.github.io/PythonDataScienceHandbook/02.07-fancy-indexing.html
  wts_use = wts_conv_kernel[np.arange(wts_conv_kernel.shape[0]), np.array(seq_idxes),:]
  # apply convolution
  score = np.sum(wts_use, axis = 0) + wts_conv_bias
  # apply batchnorm and return - cf. https://github.com/aleju/papers/blob/master/neural-nets/Batch_Normalization.md

  return g*( (score - m)/np.sqrt(v) ) + b

# given an integer,
# get the corresponding list of indices 'score_fx' will want
# e.g. int_to_indices_seq(1000) -> [0, 0, 0, 3, 3, 2, 2, 0]
def int_to_indices_seq(idx):
  assert idx >= 0 and idx < 4**8
  # convert to list of binary digits, zero-padding as needed
  idx = [q for q in bin(idx)[2:] ]
  idx = ['0']*(16 - len(idx)) + idx
  idxes_out = []
  while len(idx) > 0:
    idxes_out.append(BASEFOUR_DICT[''.join([idx.pop(0), idx.pop(0)])])
  return idxes_out

def main(fn_in, fn_out):
  f = h5py.File(fn_in, 'r')
  score_fx = partial(score_function_skeleton, f)
  idxes_list = []
  for i in range(4**8):
    idxes_list.append(int_to_indices_seq(i))
  seqs = [ ''.join([DNA[p] for p in q]) for q in idxes_list ]
  with open(fn_out, 'w') as fo:
    for (s,l) in zip(seqs, idxes_list):
      fo.write( s + ',' + ','.join([str(q) for q in score_fx(l)]) + '\n' )

if __name__ == '__main__':
  [fn_in, fn_out] = [os.path.expanduser(q) for q in sys.argv[1:]]
  main(fn_in, fn_out)
