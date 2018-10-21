# Python-based FACS-Seq means estimator
# Added 2018-10: if the output file already exists, pre-load the means;
# create a column 'RMSD' in the score file, for the RMSD of each mean estimate
# w.r.t. to the final means.

import sys
import os
from numba import cuda, float32
import numpy as np
import ConfigParser
import pandas
from math import erf, sqrt, log

REPLICATE_NAMES = ['A','B']
A_SMALL_FLOAT = -1000000.
ALLOWED_MODES = ['NUMPY','CUDA','DEBUG']

def phi(x,mu,sigma):
  # standardarize (BK)
  assert(sigma >= 0.)
  x = (x-mu)/float(sigma)
  #'Cumulative distribution function for the standard normal distribution'
  # https://stackoverflow.com/questions/809362/how-to-calculate-cumulative-normal-distribution-in-python/12955360
  return (1.0 + erf(x / sqrt(2.0))) / 2.0

def log_zeroable(x):
  if x == 0.:
    return(float('-inf'))
  return(log(x))

def isint(x):
  try:
    int(x)
  except ValueError:
    return(False)
  return(True)

# https://numba.pydata.org/numba-doc/dev/cuda/reduction.html
@cuda.reduce
def sum_reduce(a, b):
  return a + b

def progress(count, total, status=''):
  bar_len = 60
  filled_len = int(round(bar_len * count / float(total)))
  percents = round(100.0 * count / float(total), 1)
  bar = '=' * filled_len + '-' * (bar_len - filled_len)
  sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
  sys.stdout.flush() # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)

@cuda.jit
def mle_score_kernel_simple(A, B, C):
  i, j = cuda.grid(2)
  if i < C.shape[0] and j < C.shape[1]:
    tmp = 0.
    for k in range(A.shape[1]):
      tmp += A[i, k] * B[k, j]
    C[i, j] = tmp

def mle_score_wrapper(d_cell_arr, d_probs_arr, threads_per_side):
  assert(d_cell_arr.dtype == d_probs_arr.dtype)
  d_score_arr = cuda.cudadrv.devicearray.DeviceNDArray(shape = (d_cell_arr.shape[0], d_probs_arr.shape[1]), strides = d_cell_arr.strides, dtype = d_cell_arr.dtype)
  blockspergrid_score = ( ((d_score_arr.shape[0] + (threads_per_side - 1)) // threads_per_side), ((d_score_arr.shape[1] + (threads_per_side - 1)) // threads_per_side) )
  mle_score_kernel_simple[blockspergrid_score, (threads_per_side, threads_per_side)](d_cell_arr, d_probs_arr, d_score_arr)
  return(d_score_arr)

@cuda.jit
def mle_argmax_kernel(score_arr, argmax_arr, argscore_arr):
  x = cuda.grid(1)
  if x < argmax_arr.shape[0]:
    mval = A_SMALL_FLOAT
    midx = -1
    for i in range(score_arr.shape[1]):
      if score_arr[x,i] > mval:
        midx = i
        mval = score_arr[x,i]
    argmax_arr[x] = midx
    argscore_arr[x] = mval

def mle_argmax_wrapper(d_score_arr, threads_per_block):
  argmax_arr = np.ascontiguousarray(np.zeros(shape = (d_score_arr.shape[0],), dtype = 'int'))
  argscore_arr = np.ascontiguousarray(np.zeros(shape = (d_score_arr.shape[0],), dtype = 'float'))
  d_argmax_arr = cuda.to_device(argmax_arr)
  d_argscore_arr = cuda.to_device(argscore_arr)
  blockspergrid_argmax = (argmax_arr.shape[0] + (threads_per_block - 1)) // threads_per_block
  mle_argmax_kernel[blockspergrid_argmax, threads_per_block](d_score_arr, d_argmax_arr, d_argscore_arr)
  d_argmax_arr.copy_to_host(argmax_arr)
  d_argscore_arr.copy_to_host(argscore_arr)
  return(argmax_arr, argscore_arr)

def mle_score_numpy(cell_arr, probs_arr):
  score_arr = np.matmul(cell_arr, probs_arr)
  return(score_arr)

def mle_argmax_numpy(score_arr):
  argmax_arr = np.argmax(score_arr, axis = 1)
  argscore_arr = np.ascontiguousarray(np.zeros(shape = (score_arr.shape[0],), dtype = 'float'))
  for i in range(argmax_arr.shape[0]):
    argscore_arr[i] = score_arr[i, argmax_arr[i]]
  return(argmax_arr, argscore_arr)

def mle_wrapper(probs_arr, threads_per_side, threads_per_block, mode, cell_arr_np = None, d_cell_arr = None):
  if mode == 'CUDA' or mode == 'DEBUG':
    assert(d_cell_arr is not None)
    d_probs_arr = cuda.to_device(probs_arr)

  assert(np.all(probs_arr > A_SMALL_FLOAT))
  if mode == 'CUDA' or mode == 'DEBUG':
    d_score_arr = mle_score_wrapper(d_cell_arr, d_probs_arr, threads_per_side)
  if mode == 'NUMPY' or mode == 'DEBUG':
    assert(cell_arr_np is not None)
    score_arr_np = mle_score_numpy(cell_arr_np, probs_arr)

  if mode == 'CUDA' or mode == 'DEBUG':
    argmax_arr, argscore_arr = mle_argmax_wrapper(d_score_arr, threads_per_block)
  if mode == 'NUMPY' or mode == 'DEBUG':
    argmax_arr_np, argscore_arr_np = mle_argmax_numpy(score_arr_np)
  if mode == 'DEBUG':
    print('argmax_diff: ' + str(np.sum(argmax_arr - argmax_arr_np)))
    print('argscore_diff: ' + str(np.sum(argscore_arr - argscore_arr_np)))

  if mode == 'CUDA' or mode == 'DEBUG':
    score = sum_reduce(argscore_arr)
  if mode == 'NUMPY' or mode == 'DEBUG':
    score_np = np.sum(argscore_arr_np)
  if mode == 'DEBUG':
    print('final_score_diff: ' + str(score - score_np))

  if mode == 'NUMPY' or mode == 'DEBUG':
    argmax_arr = argmax_arr_np
    score = score_np

  return(argmax_arr, score)

def get_prob_arr(mean_range, bin_edges, **kwargs):
  sigma = kwargs['sigma'] # maybe eventually include other parameters to be globally optimized over, e.g. barcode mis-assignment probabilities
  fuzz = kwargs['fuzz']
  bin_edges = np.array([float('-inf')] + list(bin_edges) + [float('inf')])
  prob_arr = np.zeros(shape = (bin_edges.shape[0]-1, mean_range.shape[0]), dtype = 'float')
  for i in range(prob_arr.shape[1]): # for each choice of mean
    m = mean_range[i]
    for j in range(prob_arr.shape[0]): # for each bin
      try:
        prob = phi(bin_edges[j+1], m, sigma) - phi(bin_edges[j], m, sigma)
        prob = prob*(1-fuzz) + fuzz/prob_arr.shape[0]
        prob_arr[j,i] = log_zeroable(prob) # put the 'log' in log-likelihood
      except ValueError:
        print(phi(bin_edges[j+1], m, sigma))
        print(phi(bin_edges[j], m, sigma))
        print(phi(bin_edges[j+1], m, sigma) - phi(bin_edges[j], m, sigma))
        raise ValueError
  return(prob_arr)

def process_one_param_choice(cell_arr, mean_range, bin_edges, threads_per_side, threads_per_block, mode, verbose, **kwargs):
  d_cell_arr = None
  cell_arr_np = None
  if mode == 'CUDA' or mode == 'DEBUG':
    d_cell_arr = cuda.to_device(cell_arr)
  if mode == 'NUMPY' or mode == 'DEBUG':
    cell_arr_np = cell_arr

  prob_arr = get_prob_arr(mean_range, bin_edges, **kwargs)
  argmaxes, score = mle_wrapper(prob_arr, threads_per_side, threads_per_block, mode, d_cell_arr = d_cell_arr, cell_arr_np = cell_arr_np)
  mean_ests = mean_range[argmaxes.astype('int')]
  if verbose == 2:
    print('kwargs: ' + str(kwargs))
    print('score: ' + str(score))
  return(mean_ests, score)

def main_estimator(cell_arr, mean_range, sigma_range, fuzz_range, bin_edges, table_fn, threads_per_side, threads_per_block, mode, verbose, means_final):
  best_score = float('-inf')
  mean_ests = None
  global_ests = None
  num_choices = len(sigma_range)*len(fuzz_range)
  curr = 0
  if table_fn != None:
    tf = open(table_fn, 'w')
    if not means_final is None:
      tf.write('sigma,fuzz,score,rmsd\n')
      print('tracking RMSD')
    else:
      tf.write('sigma,fuzz,score\n')
      print('not tracking RMSD')
  try:
    for i in sigma_range:
      for j in fuzz_range:
        new_mean_ests, score = process_one_param_choice(cell_arr, mean_range, bin_edges, threads_per_side, threads_per_block, mode, verbose, sigma = i, fuzz = j)
        if table_fn != None:
          if not means_final is None:
            rmsd = np.sqrt(np.mean(np.square(new_mean_ests - means_final)))
            tf.write(str(i) + ',' + str(j) + ',' + str(score) + ',' + str(rmsd) + '\n')
          else:
            tf.write(str(i) + ',' + str(j) + ',' + str(score) + '\n')
        if score > best_score:
          mean_ests = new_mean_ests
          global_ests = {'sigma':i, 'fuzz':j}
          best_score = score
        curr += 1
        if verbose == 1:
          progress(curr, num_choices)
  finally:
    if table_fn != None:
      tf.close()
  return(mean_ests, global_ests)

def main_input_handler(read_tuple, config):
  read_arr, rep_name, means_final = read_tuple

  cell_counts = np.array([int(q) for q in config.get('sort_' + rep_name,'cell_counts').strip().split(',')])
  bin_edges = np.array([float(q) for q in config.get('sort_' + rep_name,'bin_edges').strip().split(',')])
  assert(bin_edges.shape[0] == cell_counts.shape[0] - 1)
  mean_range = np.array([float(q) for q in config.get('fitting_' + rep_name,'mean_range').strip().split(',')])
  sigma_range = np.array([float(q) for q in config.get('fitting_' + rep_name,'sigma_range').strip().split(',')])
  fuzz_range = np.array([float(q) for q in config.get('fitting_' + rep_name, 'fuzz_range').strip().split(',')])
  mean_range = np.arange(mean_range[0], mean_range[1], mean_range[2])
  sigma_range = np.arange(sigma_range[0], sigma_range[1], sigma_range[2])
  fuzz_range = np.arange(fuzz_range[0], fuzz_range[1], fuzz_range[2])
  table_fn = os.path.expanduser(config.get('fitting_' + rep_name, 'score_file'))

  threads_per_side = int(config.get('Params', 'THREADS_PER_SIDE'))
  threads_per_block = int(config.get('Params', 'THREADS_PER_BLOCK')) # mle_argmax_wrapper -> mle_wrapper
  mode = config.get('Params','MODE')
  assert mode in ALLOWED_MODES
  verbose = int(config.get('Params','VERBOSE'))

  for i in range(read_arr.shape[1]):
    read_arr[:,i] = read_arr[:,i]*cell_counts[i]/np.sum(read_arr[:,i]) # convert read counts to cell counts
  (mean_ests, global_ests) = main_estimator(read_arr, mean_range, sigma_range, fuzz_range, bin_edges, table_fn, threads_per_side, threads_per_block, mode, verbose, means_final)
  return(mean_ests, global_ests)

if __name__ == '__main__':
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.read(sys.argv[1])

  reads_fn = os.path.expanduser(config.get('Input','in_file'))
  read_df = pandas.read_csv(reads_fn)
  seqs = [q for q in read_df['Seq']]
  read_df = read_df.drop('Seq',axis=1) # drop label 0 on axis 1 (columns)

  # gotta do some data munging
  for q in read_df.columns.values:
    if not isint(q):
      read_df = read_df.drop(q, axis=1)

  colnames = np.array([int(q) for q in read_df.columns.values])
  read_arr = np.array(read_df).astype('float')
  read_arr = read_arr[:,np.argsort(colnames)]

  read_arrs = [read_arr[:,:12], read_arr[:,12:]]
  
  # added 2018-10: check whether the means already exist. If yes, each time a hyperparameter choice is tested,
  # get the RMSD of the predicted means w.r.t. the final means.

  if os.path.isfile(os.path.expanduser(config.get('Output','out_file'))):
    means_final = pandas.read_csv(os.path.expanduser(config.get('Output','out_file')))
    means_final = [means_final['Means_A'], means_final['Means_B']]
  else:
    means_final = [None, None]

  read_tuples = zip(read_arrs, REPLICATE_NAMES, means_final)
  
  mg_list = [main_input_handler(q, config) for q in read_tuples]
  means = [q[0] for q in mg_list]
  global_ests = [q[1] for q in mg_list]
  print('global_ests: ' + str(global_ests))
  out = pandas.DataFrame({'Seqs': seqs, 'Means_A': means[0], 'Means_B':means[1]})
  out.to_csv(os.path.expanduser(config.get('Output','out_file')), index = False)
