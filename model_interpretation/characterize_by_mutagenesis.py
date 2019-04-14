# Given a sequence, models, and mutable bases,
# generate all possible single mutants of the sequence.
# Test the strength of these mutants;
# output is the original sequence and strength, then a 4xN matrix with 0 in the original positions and in masks,
# and the difference (new - old) in strengths elsewhere. Rownames are ('A','C','G','T').

# models, mutable bases, and other parameters are given by a config file extracted from a filename,
# as in 'get_seq_gradients.py'.

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

def get_single_mutants(seq_in_oh, mutable):
  mutants = []
  mutants_key = []
  for m in mutable: # numerical positions of bases that can be changed
    for i in range(seq_in_oh.shape[-1]):
      if seq_in_oh[m,i] == 0.:
        mutant = np.copy(seq_in_oh)
        mutant[m,:] = 0.
        mutant[m,i] = 1.
        mutants_key.append((m,i))
        mutants.append(mutant)

  mutants_key.append((None, None))
  mutants.append(seq_in_oh)
  mutants = np.stack(mutants, axis = 0)
  return mutants, mutants_key

# Output: array 'mutants' as in get_single_mutants();
# mutants_key will now be list of (m1, i1, m2, i2) tuples.
def get_double_mutants(seq_in_oh, mutable):
  single_mutants, single_mutants_key = get_single_mutants(seq_in_oh, mutable)
  double_mutants_out = []
  double_mutants_key_out = []
  for (single, (m,i)) in zip(single_mutants, single_mutants_key):
    if m is not None and i is not None:
      mutable_inner = [q for q in mutable if q > m]
      double_mutants, double_mutants_key = get_single_mutants(single, mutable_inner)
      double_mutants_key = [(m,i,p,q) for (p,q) in double_mutants_key]
      double_mutants_out.append(double_mutants)
      double_mutants_key_out.extend(double_mutants_key)

  double_mutants_out = np.concatenate(double_mutants_out, axis = 0)
  return single_mutants, single_mutants_key, double_mutants_out, double_mutants_key_out

# Get the results of prediction (all as differences from the original 'true strength');
# get them as a large square grid with sides 0_A, 0_C, 0_G, 0_T, 1_A, etc (base values interleaved within position id's)
def get_double_mutant_pred_array(seq_orig, mutable, single_mutant_preds, single_mutants_key, double_mutant_preds, double_mutants_key):
  str_orig = single_mutant_preds[-1]
  single_dict = {p:q for (p,q) in zip(single_mutants_key, single_mutant_preds)}
  double_dict = {p:q for (p,q) in zip(double_mutants_key, double_mutant_preds)}
  # order doesn't matter for double mutants
  double_dict.update({(m2,i2,m1,i1):q for ((m1,i1,m2,i2),q) in zip(double_mutants_key, double_mutant_preds)})

  (num_pos, num_bases) = seq_orig.shape; ans_size = num_pos*num_bases
  ans = np.zeros(shape = (ans_size, ans_size))
  for i in range(num_pos):
    for j in range(num_pos):
      if i in mutable and j in mutable and i != j: # changing the same base twice is meaningless
        base_i = np.where(seq_orig[i,:])[0][0]; base_j = np.where(seq_orig[j,:])[0][0] # what where the original bases at these locations?
        for a in range(num_bases):
          for b in range(num_bases):
            pos_x = i*num_bases + a; pos_y = j*num_bases + b
            if a == base_i:
              if b == base_j:
                ans[pos_x, pos_y] = 0. # this is the WT sequence - just for clarity
              else: # b is the (single) mutant
                ans[pos_x, pos_y] = single_dict[(j,b)] - str_orig
            else:
              if b == base_j: # a is the (single) mutant
                ans[pos_x, pos_y] = single_dict[(i,a)] - str_orig
              else: # this is a true double mutant
                ans[pos_x, pos_y] = double_dict[(i,a,j,b)] - str_orig
  return ans

def get_preds(mutants, evolver):
  preds = evolver._test_sequences(np.swapaxes(mutants,1,2)) # actually reversing another 'swapaxes' in this function
  #print preds.shape # (778, 2, 9)
  preds = np.apply_along_axis(evolver.params['merge_outputs'], 1, preds)
  preds = np.apply_along_axis(evolver.params['merge_models'], 1, preds)
  return preds


def main_single(seq_in, fn_template, fn_out, loaded_models = None):
  
  seq_in_oh = one_hot_encode(seq_in).squeeze()
  #print seq_in_oh.shape # (363,4)
  (cfg, _) = build_design_cfgs.build_one_cfg(*parse_fn_in(fn_template))
  evolver = seq_evolution.seq_evolution_class(cfg, loaded_models)
  loaded_models = evolver.models

  mutants, mutants_key = get_single_mutants(seq_in_oh, evolver.mutable)
  #print mutants.shape # (778, 363, 4)
  preds = get_preds(mutants, evolver)
  
  ans = np.zeros(shape = seq_in_oh.shape)
  for ((m,i), pred) in zip(mutants_key, preds):
    if m is None or i is None:
      pass #str_orig = pred
    else:
      ans[m,i] = pred - preds[-1]
  np.savetxt(fn_out, ans, delimiter = ',')
  return loaded_models

def main_double(seq_in, fn_template, fn_out, loaded_models = None):
  seq_in_oh = one_hot_encode(seq_in).squeeze()
  (cfg, _) = build_design_cfgs.build_one_cfg(*parse_fn_in(fn_template))
  evolver = seq_evolution.seq_evolution_class(cfg, loaded_models)
  loaded_models = evolver.models
  #print seq_in_oh.shape # (363,4)
  single_mutants, single_mutants_key, double_mutants, double_mutants_key = get_double_mutants(seq_in_oh, evolver.mutable)
  single_mutant_preds = get_preds(single_mutants, evolver)
  double_mutant_preds = get_preds(double_mutants, evolver)
  arr_out = get_double_mutant_pred_array(seq_in_oh, evolver.mutable, single_mutant_preds, single_mutants_key, double_mutant_preds, double_mutants_key)
  np.savetxt(fn_out, arr_out, delimiter = ',')
  return loaded_models


if __name__ == '__main__':
  fn_template = '~/facs-seq_test/seq_design/all_FS9_seqs/0_GPD_strong_nofilter_mean_screen_0.45_selected.txt'
  seq_in = 'TACGTAAATAATTAATAGTAGTGACCGGGCCGATAGATGAGTCATTAGGGATTCCGCTCGCCCTGTGTCTGGGTGTTGCGGCATCCGGCATCCAGTAGTGGGTGTAGAATTGTGTGATAGGCATCCAGTGTCTGCATCCCAGCCACACCCCACTTTAGCACTATTTTCACCAGTGCGCCGCTCCCGTTGTCAATGGGTCTACCCCCTGTTTTCCAGGAGGTATATAAAGGAATGGTTTTTCGCGTTATCGATTTATATTATATGTTAATAAAAAATGGTATTTAATTTTTATTTCACCAAGTCCAATTCTCAATTCTCTCATAACTACATTTACTCAATGTCTAAAGGTGAAGAATTATTCAC'
  #fn_out = 'char_mut_0_0.csv'
  #main_single(seq_in, fn_template, fn_out)
  fn_out = 'char_mut_0_0_double.csv'
  main_double(seq_in, fn_template, fn_out)
