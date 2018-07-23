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
import StringIO
import oligo_design
import primer_gen

# All sequences will be padded for modeling purposes. We need to remove this pads; different libraries have different lengths.
PAD_LENS = { 'GPD': '25,25', 'ZEV': '58,58' }
# How long should the pool-specific primer sites be? For GPD, constrained by the lengtb of the promoters vs. the oligos.
POOL_TOEHOLD_LENS = { 'GPD': '17,17', 'ZEV': '25,25' }

# Copying ConfigParsers properly is a little tricky.
# Solution cf. https://stackoverflow.com/questions/23416370/manually-building-a-deep-copy-of-a-configparser-in-python-2-7
def copy_cfg(cfg):
  cfg_str = StringIO.StringIO()
  cfg.write(cfg_str)
  cfg_str.seek(0)
  cfg_out = ConfigParser.RawConfigParser()
  cfg_out.readfp(cfg_str)
  return(cfg_out)

# Given a single "master" config, modify it for one specific oligo pool.
# Keys to reset:
# pad_lens = 25,25
# fwd_pools = CCGGCGTGCGCGATACC,CAAGGCGGGTGAAGCTC,TCAATTCTTGCTACGGT
# rev_pools = GCGAGCCTTCCCAACGT,TCGGCGAGCCCAGAAGT,ATATCATCTGTCCTGAC
# assembly_id = test

def customize_cfg(cfg, table, keyname, fwd_pool_gen, rev_pool_gen, num_fwd, num_rev):
  c = copy_cfg(cfg)
  lib_id = table['Library'].unique(); assert(len(lib_id) == 1)
  c.set('Params','pad_lens',PAD_LENS[lib_id[0]])
  c.set('Params','assembly_id',c.get('Params','assembly_id') + '_' + str(keyname))
  fwd_pools = []
  for i in range(num_fwd):
    fwd_pools.append(fwd_pool_gen.next())
  rev_pools = []
  for i in range(num_rev):
    rev_pools.append(rev_pool_gen.next())
  c.set('Params','fwd_pools',','.join(fwd_pools))
  c.set('Params','rev_pools',','.join(rev_pools))
  return(c)

# generator; yields lines starting with 'start_line' (0-indexed, of course)
def start_file_at_line(fn_in, start_line):
  with open(fn_in, 'r') as fi:
    curr_line = 0
    while curr_line < start_line:
      x = fi.next()
    while True:
      yield(fi.next().strip())

if __name__ == '__main__':
  # load data
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  fn_in = os.path.expanduser(cfg.get('Files', 'table_in'))
  dat = pd.read_csv(fn_in)
  dat = dat.loc[pd.notna(dat['Last.Control.Use']),:]
  pool_names = dat['Last.Control.Use'].unique()
  dfs = { q: dat.loc[dat['Last.Control.Use'] == q,:] for q in pool_names }
  cfgs = {}
  fn_fwd_pool = os.path.expanduser(cfg.get('Pools', 'fwd_pool_fn'))
  start_fwd_pool = int(cfg.get('Pools','fwd_pool_start'))
  fn_rev_pool = os.path.expanduser(cfg.get('Pools', 'rev_pool_fn'))
  start_rev_pool = int(cfg.get('Pools','rev_pool_start'))
  fwd_pool_gen = start_file_at_line(fn_fwd_pool, start_fwd_pool)  # open the file; start at a specified line
  rev_pool_gen = start_file_at_line(fn_rev_pool, start_rev_pool)
  for q in pool_names:
    cfgs[q] = customize_cfg(copy_cfg(cfg), dfs[q], q, fwd_pool_gen, rev_pool_gen, num_fwd = 2, num_rev = 2)
  finals = []; rejects = []
  for q in pool_names:
    d = dfs[q]; c = cfgs[q]
    seqs = d['Seq'].tolist()
    final_oligos, rejected = oligo_design.seqs_to_oligos(seqs, c)
    finals.extend(final_oligos)
    rejects.extend(rejected)

  with open(os.path.expanduser(cfg.get('Files','oligo_fn')),'w') as fo:
    for name, oligo in finals:
      fo.write(name + '\n')
      fo.write(oligo + '\n')
  with open(os.path.expanduser(cfg.get('Files','rejected_fn')),'w') as fr:
    for seq in rejected:
      fr.write(seq + '\n')
