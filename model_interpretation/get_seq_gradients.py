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
import pandas as pd

GET_OBJECTIVE = {'strong':'Strength', 'induced':'Induced strength', 'AR':'AR'}
GET_FILTER = {'gcfilter': True, 'nofilter': False}
GET_FUNCTION = {'mean': 'Mean', '1sd': 'Mean-sd'}
GET_STRATEGY = {'screen': 'Screening', 'evolve-thresh': 'Evolution to threshold', 'evolve-cycle': 'Evolution: cycle-limited', 
                'gradient-thresh': 'Gradient to threshold', 'gradient-cycle': 'Gradient: cycle-limited'}

# fn_in should encode design settings needed to generate the config: exploit this
def parse_fn_in(fn_in):
  fn_in = fn_in.strip().split('/')[-1]
  fn_in = fn_in.split('.')[0]
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
  cfg = build_design_cfgs.build_one_cfg(parse_fn_in(fn_in))
  df_in = pd.read_csv(fn_in)

  # import the sequences from fn_in
  seqs_oh = one_hot_encode(df_in['Seqs'])
  evolver = seq_gradient_evolution.seq_evolution_class_gradient(cfg)
  evolver.seqs = seqs_oh
  _, grads = evolver.losses_and_grads()
  np.save(grads, fn_out)