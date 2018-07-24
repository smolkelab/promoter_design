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

# All sequences will be padded for modeling purposes. We need to remove this pads; different libraries have different lengths.
PAD_LENS = { 'GPD': '25,25', 'ZEV': '58,58' }
# How long should the pool-specific primer sites be? For GPD, constrained by the lengtb of the promoters vs. the oligos.
POOL_TOEHOLD_LENS = { 'GPD': '17,17', 'ZEV': '25,25' }

DNA_C = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def rc(seq):
  return(''.join([DNA_C[q] for q in seq])[::-1])

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
def start_toeholds_at_line(fn_in, start_toe):
  df = pd.read_csv(fn_in)
  toeholds = df['Toehold'].tolist()
  curr_toe = start_toe
  while True:
    yield(toeholds[curr_toe])
    curr_toe += 1
	
def toe_to_full_dict(fn_in):
  df = pd.read_csv(fn_in)
  toeholds = df['Toehold'].tolist()
  fulls = df['Full'].tolist()
  return(dict(zip(toeholds, fulls)))
	
# Given a Pandas DF, write the oligos and primers to be synthesized to files.
# Take reverse complements where needed.
def table_to_oligos(table, fn_oligo, fn_primer, fn_fwd_pool, fn_rev_pool):
  # setup to map oligo toeholds to full sequences
  fwd_dict = toe_to_full_dict(fn_fwd_pool)
  rev_dict = toe_to_full_dict(fn_rev_pool)
  table_zipped = zip(table['Experiment'],table['pool_id'],table['seq_id'],table['fwd_oligos'],table['rev_oligos'], table['fwd_pool'], table['rev_pool'] )
  prev_exp = None; prev_pool = None
  with open(fn_oligo, 'w') as fo, open(fn_primer, 'w') as fp:
    # handle the constant primers
    fp.write('>Const|F\n')
    fp.write(table['fwd_const'].tolist()[0] + '\n')
    fp.write('>Const|R\n')
    fp.write(rc(table['rev_const'].tolist()[0]) + '\n')
    for (exp, pool_id, seq_id, fwd_oligo, rev_oligo, fwd_pool, rev_pool) in table_zipped:
      name_stem = '>' + '|'.join([exp, str(pool_id), str(seq_id)])
      fo.write(name_stem + '|F\n')
      fo.write(fwd_oligo + '\n')
      fo.write(name_stem + '|R\n')
      fo.write(rev_oligo + '\n')
      if exp != prev_exp or pool_id != prev_pool:
        oligo_name_stem = '>' + '|'.join([exp, str(pool_id)])
        fp.write(oligo_name_stem + '|F\n')
        fp.write(rc(fwd_dict[fwd_pool]) + '\n')
        fp.write(oligo_name_stem + '|R\n')
        fp.write(rev_dict[rev_pool] + '\n')
      if exp != prev_exp:
        prev_exp = exp
      if prev_pool != pool_id:
        prev_pool = pool_id

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
  fwd_pool_gen = start_toeholds_at_line(fn_fwd_pool, start_fwd_pool)
  rev_pool_gen = start_toeholds_at_line(fn_rev_pool, start_rev_pool)
  for q in pool_names:
    cfgs[q] = customize_cfg(copy_cfg(cfg), dfs[q], q, fwd_pool_gen, rev_pool_gen, num_fwd = 2, num_rev = 2)
  tables = []
  rejects = []
  for q in pool_names:
    d = dfs[q]; c = cfgs[q]
    seqs = d['Seq'].tolist()
    table, rejected = oligo_design.seqs_to_df(seqs, c)
    tables.append(table)
    rejects.extend(rejected)

  table_out = pd.concat(tables)
  table_out.to_csv(os.path.expanduser(cfg.get('Files','table_out')))
  with open(os.path.expanduser(cfg.get('Files','rejected_fn')),'w') as fr:
    for seq in rejected:
      fr.write(seq + '\n')
  table_to_oligos(table_out, os.path.expanduser(cfg.get('Files','oligo_fn')), os.path.expanduser(cfg.get('Files','primer_fn')), fn_fwd_pool, fn_rev_pool)
