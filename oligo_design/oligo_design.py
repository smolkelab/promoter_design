# given list of target seqs (col. headed 'Seqs' in CSV), 
# config containing: constant pad positions to remove from both ends, required len. of oligo, number of pools, 
# get DF with columns:
# pool_id: name of pool
# seq_id: name of this sequence within the pool
# oligo_f: 'fwd' oligo to synthesize
# oligo_r: 'rev' oligo to synthesize
# pool_primer_f: 'fwd' primer for 'this' pool (binds oligo_r) (+)
# pool_primer_r: 'rev' primer for 'this' pool (binds oligo_f) (+) (!)
# gg_site: golden gate site for this construct
# gg_prod: final output of the construct
# const_primer_f: constant forward primer (binds oligo_f)
# const_primer_r: constant reverse primer (binds oligo_r)

# Process this DF with methods in 'validate_and_get_primers.py'

import sys
import os
import pandas as pd
import numpy as np
import random
from numpy.random import choice
import ConfigParser

DNA = ['A','C','G','T']
DNA_C = {'A':'T', 'C':'G', 'G':'C','T':'A'}

# Represent a set of sequences, divided into two oligos, to be joined by Golden Gate assembly.
# Track the use of Golden Gate sites so that no assembly shares a site.

class sequence_pool(object):
  def __init__(self, num_seqs, site_start, site_end, gg_len):
    # the extra column is so that one col. will always have 0 - covers a potential corner case
    self.gg_mat = np.zeros(shape = (num_seqs, len(DNA)**gg_len + 1), dtype = 'bool')
    self.ss = site_start
    self.se = site_end
    self.gg_len = gg_len
    self.gg_to_bin_id = {} # map 4-base GG sequences to columns in self.gg_mat, and vice versa
    self.curr_col = 0
    self.seq_to_useable_ggs = {} # map a sequence to a (gg_seq -> gg_start) dictionary
    self.seqs = []

  def add(self, seq):
    seq_id = len(self.seqs)
    self.seqs.append(seq)
    useable_ggs = []
    useable_starts = []
    gg_start = self.ss
    gg_end = gg_start + self.gg_len
    # what GG sites could work for this sequence?
    while gg_end <= self.se:
      useable_ggs.append(seq[gg_start:gg_end])
      useable_starts.append(gg_start)
      gg_start += 1; gg_end += 1
    self.seq_to_useable_ggs[seq] = dict(zip(useable_ggs, useable_starts))
    # record this sequence's useable GG sites
    for g in useable_ggs:
      if g not in self.gg_to_bin_id:
        self.gg_to_bin_id[g] = self.curr_col; self.gg_to_bin_id[self.curr_col] = g
        self.curr_col += 1
      self.gg_mat[seq_id, self.gg_to_bin_id[g]] = True

  # assign GGs to sequences; try to assign as many as possible.
  # Heuristic: find the GG used by the fewest sequences, and assign to that one.
  # Take that sequence out of consideration; continue until no remaining bins can be used.
  def assign_seqs(self):
    accepted = np.zeros(shape = (len(self.seqs),), dtype = 'bool')
    ans = []
    while True:
      bin_counts = np.apply_along_axis(np.sum, 0, self.gg_mat)
      u_cts, u_idxes = np.unique(bin_counts, return_index = True)
      if len(u_cts) == 1: # they're all 0
          break
      bin_id = u_idxes[1] # u_cts[0] will be 0
      seq_id = np.min(np.where(self.gg_mat[:,bin_id])[0])
      accepted[seq_id] = True
      self.gg_mat[seq_id,:] = False # this sequence can be removed
      self.gg_mat[:,bin_id] = False # this bin is used
      g = self.gg_to_bin_id[bin_id] # this is a bidirectional dictionary; get the GG sequence
      # get the actual seq and starting position
      seq = self.seqs[seq_id]
      g_s = self.seq_to_useable_ggs[seq][g]
      ans.append((seq, g_s))
    rejected = np.array(self.seqs)[np.logical_not(accepted)]
    return(ans, rejected.tolist())
    
# generate a random pad sequence that won't contain cut sites.
# Context is required to ensure we don't finish a partial site that started before the pad.
def safe_pad(seq, pad_len, forbidden_site_list, is_fwd):
  max_forbid_len = max([len(q) for q in forbidden_site_list])
  for (i,f) in enumerate(forbidden_site_list): # strip N's from binding sites
    f = f.strip('N'); assert('N' not in f)
    forbidden_site_list[i] = f
  if is_fwd:
    test_seq = seq + ''.join( np.random.choice(np.array(DNA), size = pad_len, replace = True).tolist() )
    test_seq_sitecheck = test_seq[-(pad_len + max_forbid_len):]
  else:
    test_seq = ''.join( np.random.choice(np.array(DNA), size = pad_len, replace = True).tolist() ) + seq
    test_seq_sitecheck = test_seq[:(pad_len + max_forbid_len)]
  if any([q in test_seq_sitecheck for q in forbidden_site_list]):
    # recursion approach. Keep trying until a forbidden-site-free sequence is generated.
    return(safe_pad(seq, pad_len, forbidden_site_list, is_fwd))
  return(test_seq)

# Given a list of sequences and a dictionary of params, get the output table as a Pandas DataFrame.
def build_pools_table(seqs, params):
  oligo_len = int(params['oligo_len'])
  gg_site_len = int(params['gg_site_len'])
  fwd_pools = params['fwd_pools']
  rev_pools = params['rev_pools']
  len_fwd_seq = oligo_len - len(params['fwd_toehold']) - gg_site_len - len(fwd_pools[0])
  len_rev_seq = oligo_len - len(params['rev_toehold']) - gg_site_len - len(rev_pools[0])
  # golden gate site can start at len_seqs - len_rev_seq - len_gg_site; must end at/before len_fwd_seq + len_gg_site
  site_start = len(seqs[0]) - len_rev_seq - gg_site_len
  site_end = len_fwd_seq + gg_site_len
  num_pools = len(fwd_pools)
  pools = []
  for i in range(num_pools):
    if len(seqs) == 0:
      break
    this_pool = sequence_pool(len(seqs), site_start, site_end, gg_site_len)
    for s in seqs:
      this_pool.add(s)
    ans, seqs = this_pool.assign_seqs()
    pools.append(ans)
  ### Build the output DataFrame from 'pools' ###
  # 'pools' is a list of lists of tuples; each inner ('pool') list corresponds to a pool.
  # Each pool list contains tuples (seq, g_s); each 'seq' is an original sequence
  # (the thing to be constructed, without any padding), each 'g_s' is the starting position for the GG site.
  # Create a table with columns 'Design', 'gg_start', 'fwd_pool', 'rev_pool'
  dfs = []
  for (i,q) in enumerate(pools):
    pool_dict = {'Design':[], 'gg_start':[], 'fwd_pool':[],'rev_pool':[], 'Experiment':[],'pool_id':[]}
    for (seq, g_s) in q:
      pool_dict['Design'].append(seq)
      pool_dict['gg_start'].append(g_s)
      pool_dict['fwd_pool'].append(fwd_pools[i])
      pool_dict['rev_pool'].append(rev_pools[i])
      pool_dict['Experiment'].append(params['assembly_id'])
      pool_dict['pool_id'].append(i)
    dfs.append(pd.DataFrame(pool_dict))
  table = pd.concat(dfs)
  
  return(table, seqs) # 'seqs' is now the rejected sequences

def fill_one_oligo(tuple_in, params):
  seq, gg_start, fwd_pool, rev_pool = tuple_in
  oligo_len = int(params['oligo_len'])
  gg_site_len = int(params['gg_site_len'])
  forbidden_site_list = [params['fwd_gg_cut'], params['rev_gg_cut']]
  fwd_prefix = params['fwd_toehold']
  fwd_postfix = fwd_pool # NB: GG site has to be included here!
  rev_prefix = rev_pool # see above
  rev_postfix = params['rev_toehold']
  fwd_seq = seq[:(gg_start+gg_site_len)]
  rev_seq = seq[gg_start:]
  assert(len(fwd_seq) <= oligo_len)
  assert(len(rev_seq) <= oligo_len)
  fwd_seq = safe_pad(fwd_seq, oligo_len - len(fwd_seq), forbidden_site_list, True)
  rev_seq = safe_pad(rev_seq, oligo_len - len(rev_seq), forbidden_site_list, False)
  return(fwd_seq, rev_seq)

# create final oligos; add PCR regions and padding if needed
def fill_pools_table(table, params):
  forbidden_site_list = [params['fwd_gg_cut'], params['rev_gg_cut']]
  oligo_len = int(params['oligo_len'])
  gg_site_len = int(params['gg_site_len'])
  table_zipped = zip(table['Design'], table['gg_start'], table['fwd_pool'], table['rev_pool'])
  fwd_seqs = []; rev_seqs = []
  for q in table_zipped:
    f, r = fill_one_oligo(q, params)
    fwd_seqs.append(f); rev_seqs.append(r)
  table['fwd_oligos'] = fwd_seqs
  table['rev_oligos'] = rev_seqs
  return(table)
  
'''
      # naming convention for oligos: assembly id, pool id (number in this assembly), seq id (number in this pool), F/R 
      name_stem = '>' + '|'.join([params['assembly_id'], str(i), str(j)])
      final_oligos.append((name_stem + '|F', fwd_seq))
      final_oligos.append((name_stem + '|R', rev_seq))
  return(final_oligos)'''
    
# For modularity: given a list of sequences and a config, get the output table as a Pandas DataFrame.
def seqs_to_df(seqs, cfg):
  assert(all([len(q) == len(seqs[0]) for q in seqs]))
  params = dict(cfg.items('Params'))
  random_seed = int(params['random_seed'])
  random.seed(random_seed); np.random.seed(random_seed)
  # remove padding from sequences
  [pad_len_f, pad_len_r] = [int(q) for q in params['pad_lens'].strip().split(',')]
  seqs = [q[pad_len_f:-pad_len_r] for q in seqs]
  # get distinguishing sequences for pool PCRs
  fwd_pools = params['fwd_pools'].strip().split(',')
  rev_pools = params['rev_pools'].strip().split(',')
  assert(len(fwd_pools) == len(rev_pools))
  assert(all([len(q) == len(fwd_pools[0]) for q in fwd_pools]))
  assert(all([len(q) == len(rev_pools[0]) for q in rev_pools]))
  params['fwd_pools'] = fwd_pools; params['rev_pools'] = rev_pools
  #pools, rejected = build_pools(seqs, params)
  #final_oligos = fill_oligos(pools, params)
  oligo_table, rejected = build_pools_table(seqs, params)
  oligo_table = fill_pools_table(oligo_table, params)
  return(oligo_table, rejected)
  #return(final_oligos, rejected)

def main(cfg):
  seqs = list(pd.read_csv(os.path.expanduser(cfg.get('Files','selected_fn')))['Seqs'])
  oligo_table, rejected = seqs_to_oligos(seqs, cfg)
  fo = os.path.expanduser(cfg.get('Files','oligo_fn'))
  oligo_table.to_csv(fo)
  
  '''with open(os.path.expanduser(cfg.get('Files','oligo_fn')),'w') as fo:
    for name, oligo in final_oligos:
      fo.write(name + '\n')
      fo.write(oligo + '\n')'''
  with open(os.path.expanduser(cfg.get('Files','rejected_fn')),'w') as fr:
    for seq in rejected:
      fr.write(seq + '\n')
  
if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  main(cfg)
  
