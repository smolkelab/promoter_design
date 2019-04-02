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

def get_single_mutants(seq_in_oh, mutable, keep_orig = True):
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
  if keep_orig:
    mutants_key.append((None, None))
    mutants.append(seq_in_oh)
  mutants = np.stack(mutants, axis = 0)
  return mutants, mutants_key

def main(seq_in, fn_template, fn_out, as_diff = True):
  (cfg, _) = build_design_cfgs.build_one_cfg(*parse_fn_in(fn_template))
  evolver = seq_evolution.seq_evolution_class(cfg)
  seq_in_oh = one_hot_encode(seq_in).squeeze()
  #print seq_in_oh.shape # (363,4)
  mutants, mutants_key = get_single_mutants(seq_in_oh, evolver.mutable, keep_orig = as_diff)
  #print mutants.shape # (778, 363, 4)

  preds = evolver._test_sequences(np.swapaxes(mutants,1,2)) # actually reversing another 'swapaxes' in this function
  #print preds.shape # (778, 2, 9)
  preds = np.apply_along_axis(evolver.params['merge_outputs'], 1, preds)
  preds = np.apply_along_axis(evolver.params['merge_models'], 1, preds)

  if as_diff:
    str_orig = preds[-1]

  ans = np.zeros(shape = seq_in_oh.shape)
  for ((m,i), pred) in zip(mutants_key, preds):
    if m is None or i is None:
      pass #str_orig = pred
    else:
      if as_diff:
        ans[m,i] = pred - str_orig
      else:
        ans[m,i] = pred

  np.savetxt(fn_out, ans, delimiter = ',')

if __name__ == '__main__':
  fn_template = '~/facs-seq_test/seq_design/all_FS9_seqs/0_GPD_strong_nofilter_mean_screen_0.45_selected.txt'
  seq_in = 'TACGTAAATAATTAATAGTAGTGACCGGGCCGATAGATGAGTCATTAGGGATTCCGCTCGCCCTGTGTCTGGGTGTTGCGGCATCCGGCATCCAGTAGTGGGTGTAGAATTGTGTGATAGGCATCCAGTGTCTGCATCCCAGCCACACCCCACTTTAGCACTATTTTCACCAGTGCGCCGCTCCCGTTGTCAATGGGTCTACCCCCTGTTTTCCAGGAGGTATATAAAGGAATGGTTTTTCGCGTTATCGATTTATATTATATGTTAATAAAAAATGGTATTTAATTTTTATTTCACCAAGTCCAATTCTCAATTCTCTCATAACTACATTTACTCAATGTCTAAAGGTGAAGAATTATTCAC'
  fn_out = 'char_mut_0_0.csv'
  main(seq_in, fn_template, fn_out)
