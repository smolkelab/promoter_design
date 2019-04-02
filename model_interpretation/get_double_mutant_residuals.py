# Given a matrix of 'raw' double mutant score diffs as from 'characterize_by_mutagenesis.py',
# get a matrix reducing the 16 values from each pair of bases tested to 1.
# Reduction approach: assume that for independent bases, score of double mutant = sum of single mutants.
# get (double mut - sum(single muts)) for the 9 double mutants in the grid; return the one with the greatest absolute value.

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

def reduce_grid(grid_in, x_true, y_true):
  grid_ref = np.zeros(grid_in.shape)
  for i in range(grid_ref.shape[0]):
    for j in range(grid_ref.shape[1]):
      grid_ref[i,j] = grid_in[x_true,j] + grid_in[i,y_true]
  return np.max(grid_in - grid_ref)

def main(seq_in, fn_in, fn_out):
  seq_in_oh = one_hot_encode(seq_in).squeeze()
  dat_in = np.loadtxt(fn_in, delimiter = ',')
  ans = np.zeros(shape = (seq_in_oh.shape[0], seq_in_oh.shape[0]))
  for i in range(ans.shape[0]):
    for j in range(ans.shape[1]):
      grid_in = dat_in[(4*i):(4*(i+1)), (4*j):(4*(j+1))]
      i_coor = np.where(seq_in_oh[i,:])[0][0]
      j_coor = np.where(seq_in_oh[j,:])[0][0]
      ans[i,j] = reduce_grid(grid_in, i_coor, j_coor)

  np.savetxt(fn_out, ans, delimiter = ',')

if __name__ == '__main__':
  seq_in = 'TACGTAAATAATTAATAGTAGTGACCGGGCCGATAGATGAGTCATTAGGGATTCCGCTCGCCCTGTGTCTGGGTGTTGCGGCATCCGGCATCCAGTAGTGGGTGTAGAATTGTGTGATAGGCATCCAGTGTCTGCATCCCAGCCACACCCCACTTTAGCACTATTTTCACCAGTGCGCCGCTCCCGTTGTCAATGGGTCTACCCCCTGTTTTCCAGGAGGTATATAAAGGAATGGTTTTTCGCGTTATCGATTTATATTATATGTTAATAAAAAATGGTATTTAATTTTTATTTCACCAAGTCCAATTCTCAATTCTCTCATAACTACATTTACTCAATGTCTAAAGGTGAAGAATTATTCAC'
  fn_in = 'char_mut_0_0_double.csv'
  fn_out = 'char_mut_0_0_reduced.csv'
  main(seq_in, fn_in, fn_out)

