# Workflow is Desktop/pred_test -> Data/2018-07/Select_control_seqs.R -> pred_test/merge_means_predictions.py (-> Data/2018-07/Visualize_control_seqs.R)
# Use the output 'means_controls_merged.csv', instead of the individual control CSVs - the merged output accounts for the same sequence appearing in multiple
# sets, eliminating redundancy. If sequences are used more than once, just synthesize it with the last pool it's represented in (the column 'Last.Control.Use').
# Finish control design by running 'oligo_design' on subsets:
# Break the controls into subsets: each control pool can get its own set of oligos.
# Do each pool in duplicate; try to have different GG sites in each as well.
# For each pool, generate a config as used by 'oligo_design.py', and run oligo_design.main()
# Merge the outputs.

import sys
import os
import numpy as np
import pandas as pd
import ConfigParser
import copy
import oligo_design
import primer_gen

# All sequences will be padded for modeling purposes. We need to remove this pads; different libraries have different lengths.
PAD_LENS = { 'GPD': '25,25', 'ZEV': '58,58' }
# How long should the pool-specific primer sites be? For GPD, constrained by the lengtb of the promoters vs. the oligos.
POOL_TOEHOLD_LENS = { 'GPD': '17,17', 'ZEV': '25,25' }

# Given a single "master" config, modify it for one specific oligo pool.
# Keys to reset:
# pad_lens = 25,25
# fwd_pools = CCGGCGTGCGCGATACC,CAAGGCGGGTGAAGCTC,TCAATTCTTGCTACGGT
# rev_pools = GCGAGCCTTCCCAACGT,TCGGCGAGCCCAGAAGT,ATATCATCTGTCCTGAC
# assembly_id = test

def customize_cfg(cfg, table, keyname, primer_gen):
  c = copy.copy(cfg)
  lib_id = table['Library'].unique(); assert(len(lib_id) == 1)
  c.set('Params','pad_lens',PAD_LENS[lib_id[0]])
  c.set('Params','assembly_id',c.get('Params','assembly_id') + '_' + str(keyname))
  fwd_pools = # TBD: use a 'primer_gen' object to generate these primers
  rev_pools = 
  c.set('Params','fwd_pools',fwd_pools)
  c.set('Params','rev_pools',rev_pools)
  return(c)

if __name__ == '__main__':
  # load data
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  fn_in = os.path.expanduser(cfg.get('Files', 'table_in'))
  dat = pd.read_csv(fn_in)
  dat = dat.loc[pd.notna(dat['Last.Control.Use']),:]
  pool_names = dat['Last.Control.Use'].unique()
  dfs = { q: dat.loc[dat['Last.Control.Use'] == q,:] for q in pool_names }
  cfgs = { q: customize_cfg(cfg, dfs[q], q) for q in pool_names }
  finals = []; rejects = []
  for q in pool_names:
    d = dfs[q]; c = cfgs[q]
    seqs = d['Seq']
	final_oligos, rejected = oligo_design.seqs_to_oligos(seqs, c)
    finals.extend(final_oligos)
	rejects.extend(rejected)

  with open(os.path.expanduser(cfg.get('Files','oligo_fn')),'w') as fo:
    for name, oligo in final_oligos:
      fo.write(name + '\n')
      fo.write(oligo + '\n')
  with open(os.path.expanduser(cfg.get('Files','rejected_fn')),'w') as fr:
    for seq in rejected:
      fr.write(seq + '\n')