# Given a CSV with header 'Seqs,Scores',
# and the parameters used to generate the config file used to evolve those sequences,
# output a csv with header 'Seqs,Scores,0...(len[seqs]-1)',
# where new fields are the L2 norm of the gradient for each sequence at that position,
# and values outside the 'mutable mask' are zeroed out.

# assume we're executing from 'facs-seq/model_interpretation'

import sys
sys.path.append('../seq_design')
import seq_gradient_evolution
sys.path.append('../models')
from model_trainer import one_hot_encode
import build_design_cfgs
import ConfigParser
import numpy as np
import pandas as pd

GET_OBJECTIVE = {'strong':'Strength', 'induced':'Induced strength', 'AR':'AR'}
GET_FILTER = {'gcfilter': True, 'nofilter': False}
GET_FUNCTION = {'mean': 'Mean', '1sd': 'Mean-sd'}
GET_STRATEGY = {'screen': 'Screening', 'evolve-thresh': 'Evolution to threshold', 'evolve-cycle': 'Evolution: cycle-limited', 
                'gradient-thresh': 'Gradient to threshold', 'gradient-cycle': 'Gradient: cycle-limited'}

# fn_in should encode design settings needed to generate the config: exploit this
def parse_fn_in(fn_in):
  fn_in = fn_in.strip().split('/')[-1]
  fn_in = '.'.join(fn_in.split('.')[:-1])
  [i, promoter, objective, filter, function, strategy, threshold] = fn_in.split('_')[:7]
  i = int(i)
  objective = GET_OBJECTIVE[objective]
  filter = GET_FILTER[filter]
  function = GET_FUNCTION[function]
  strategy = GET_STRATEGY[strategy]
  threshold = float(threshold)
  return i, promoter, objective, filter, function, strategy, threshold

def main(fn_in, fn_out):
  # seq_gradient_evolution will expect a ConfigParser - generate this here
  (cfg, _) = build_design_cfgs.build_one_cfg(*parse_fn_in(fn_in))
  df_in = pd.read_csv(fn_in)

  # import the sequences from fn_in
  seqs_oh = one_hot_encode(df_in['Seqs'])
  #print seqs_oh.shape # (112, 363, 4)
  evolver = seq_gradient_evolution.seq_evolution_class_gradient(cfg)
  seqs_oh = seqs_oh[:,:seqs_oh.shape[1] - evolver.shift + 1,:]
  _, grads = evolver.losses_and_grads([seqs_oh])

  #print grads.shape # (112, 363, 4)
  grads = np.apply_along_axis(np.linalg.norm, -1, grads)

  # mask out the constant bases
  grads = np.multiply(grads, evolver.mutable_mask[...,0])

  np.save(fn_out, grads)

if __name__ == '__main__':
  #fn_in = '/home/benkotopka/facs-seq_test/seq_designs/all_FS9_seqs/16_ZEV_induced_gcfilter_1sd_evolve-thresh_1.6_selected.txt'
  fn_in = '/home/benkotopka/facs-seq_test/seq_designs/all_FS9_seqs/0_GPD_strong_nofilter_mean_screen_0.45_selected.txt'
  fn_out = 'tmp'
  main(fn_in, fn_out)
