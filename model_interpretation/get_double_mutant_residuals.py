# Given a matrix of 'raw' double mutant score diffs as from 'characterize_by_mutagenesis.py',
# get a matrix reducing the 16 values from each pair of bases tested to 1.
# Reduction approach: assume that for independent bases, score of double mutant = sum of single mutants.
# approach 0:
# get (double mut - sum(single muts)) for the 9 double mutants in the grid; return the one with the greatest absolute value.
# approach 1: get cov(true, predicted) for the 9 double mutants
# approach 2: get correlation(true, predicted) for the 9 double mutants
# approach 3: get L2 norm of true - predicted for the 9 double mutants

import sys
import os
sys.path.append('../seq_design')
import seq_evolution
sys.path.append('../models')
from model_trainer import one_hot_encode
import build_design_cfgs
import ConfigParser
import numpy as np
import pandas as pd
from get_seq_gradients import parse_fn_in

def merge_greatest(grid_in, grid_ref, x_true, y_true):
  return np.max(grid_in - grid_ref)

def merge_cov(grid_in, grid_ref, x_true, y_true):
  grid_in = np.delete(grid_in, x_true, axis = 0); grid_in = np.delete(grid_in, y_true, axis = 1); grid_in = grid_in.flatten()
  grid_ref = np.delete(grid_ref, x_true, axis = 0); grid_ref = np.delete(grid_ref, y_true, axis = 1); grid_ref = grid_ref.flatten()
  return np.cov(grid_in, grid_ref)[0,1]

def merge_corr(grid_in, grid_ref, x_true, y_true):
  grid_in = np.delete(grid_in, x_true, axis = 0); grid_in = np.delete(grid_in, y_true, axis = 1); grid_in = grid_in.flatten()
  grid_ref = np.delete(grid_ref, x_true, axis = 0); grid_ref = np.delete(grid_ref, y_true, axis = 1); grid_ref = grid_ref.flatten()
  return np.corrcoef(grid_in, grid_ref)[0,1]

def merge_l2(grid_in, grid_ref, x_true, y_true):
  grid_in = np.delete(grid_in, x_true, axis = 0); grid_in = np.delete(grid_in, y_true, axis = 1); grid_in = grid_in.flatten()
  grid_ref = np.delete(grid_ref, x_true, axis = 0); grid_ref = np.delete(grid_ref, y_true, axis = 1); grid_ref = grid_ref.flatten()
  return np.sqrt(np.sum((grid_in - grid_ref)**2))

MERGE_DICT = {0: merge_greatest, 1: merge_cov, 2: merge_corr, 3: merge_l2}

def reduce_grid(grid_in, x_true, y_true, merge_fx):
  grid_ref = np.zeros(grid_in.shape)
  for i in range(grid_ref.shape[0]):
    for j in range(grid_ref.shape[1]):
      grid_ref[i,j] = grid_in[x_true,j] + grid_in[i,y_true]
  return merge_fx(grid_in, grid_ref, x_true, y_true)

def main(seq_in, fn_in, fn_out, merge_fx):
  seq_in_oh = one_hot_encode(seq_in).squeeze()
  dat_in = np.loadtxt(fn_in, delimiter = ',')
  ans = np.zeros(shape = (seq_in_oh.shape[0], seq_in_oh.shape[0]))
  for i in range(ans.shape[0]):
    for j in range(ans.shape[1]):
      grid_in = dat_in[(4*i):(4*(i+1)), (4*j):(4*(j+1))]
      i_coor = np.where(seq_in_oh[i,:])[0][0]
      j_coor = np.where(seq_in_oh[j,:])[0][0]
      ans[i,j] = reduce_grid(grid_in, i_coor, j_coor, merge_fx)

  np.savetxt(fn_out, ans, delimiter = ',')

if __name__ == '__main__':
  seq_in = 'TACGTAAATAATTAATAGTAGTGACCGGGCCGATAGATGAGTCATTAGGGATTCCGCTCGCCCTGTGTCTGGGTGTTGCGGCATCCGGCATCCAGTAGTGGGTGTAGAATTGTGTGATAGGCATCCAGTGTCTGCATCCCAGCCACACCCCACTTTAGCACTATTTTCACCAGTGCGCCGCTCCCGTTGTCAATGGGTCTACCCCCTGTTTTCCAGGAGGTATATAAAGGAATGGTTTTTCGCGTTATCGATTTATATTATATGTTAATAAAAAATGGTATTTAATTTTTATTTCACCAAGTCCAATTCTCAATTCTCTCATAACTACATTTACTCAATGTCTAAAGGTGAAGAATTATTCAC'
  fn_in = 'char_mut_0_0_double.csv'
  fns_out = ('char_mut_0_0_reduced_greatest.csv', 'char_mut_0_0_reduced_cov.csv', 'char_mut_0_0_reduced_corr.csv', 'char_mut_0_0_reduced_l2.csv')
  merge_id = (0,1,2,3)
  for (p,q) in zip(fns_out, merge_id):
    main(seq_in, fn_in, p, MERGE_DICT[q])

