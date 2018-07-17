# given list of target seqs (col. headed 'Seqs' in CSV), 
# config containing: constant pad positions to remove from both ends, required len. of oligo, number of pools, 
# get .fa file with oligos needed to build seq, .fa file with oligos needed to order (e.g. PCR), table of PCR primers for each pool
import sys
import os
import pandas
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

def build_pools(seqs, params):
  oligo_len = int(params['oligo_len'])
  gg_site_len = int(params['gg_site_len'])
  fwd_pools = params['fwd_pools']
  rev_pools = params['rev_pools']
  len_fwd_seq = oligo_len - len(params['fwd_toehold']) - gg_site_len - len(params['fwd_gg_cut']) - len(fwd_pools[0])
  len_rev_seq = oligo_len - len(params['rev_toehold']) - gg_site_len - len(params['rev_gg_cut']) - len(rev_pools[0])
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
  return(pools, seqs) # 'seqs' is now the rejected sequences

# create final oligos; add PCR regions and padding if needed
def fill_oligos(pools, params):
  fwd_pools = params['fwd_pools'][:len(pools)]
  rev_pools = params['rev_pools'][:len(pools)]
  oligo_len = int(params['oligo_len'])
  gg_site_len = int(params['gg_site_len'])
  final_oligos = [] # fill with (name, seq) tuples
  for i, (pool, fwd_unique, rev_unique) in enumerate(zip(pools, fwd_pools, rev_pools)):
    fwd_prefix = params['fwd_toehold']
    fwd_postfix = params['fwd_gg_cut'] + fwd_unique
    rev_prefix = rev_unique + params['rev_gg_cut']
    rev_postfix = params['rev_toehold']
    for j, (seq, g_s) in enumerate(pool):
      # split the sequence by its GG site
      fwd_seq = seq[:(g_s+gg_site_len)]
      rev_seq = seq[g_s:]
      fwd_seq = fwd_prefix + fwd_seq + fwd_postfix
      rev_seq = rev_prefix + rev_seq + rev_postfix
      assert(len(fwd_seq) <= oligo_len)
      assert(len(rev_seq) <= oligo_len)
      # pad to full length if needed
      fwd_seq = fwd_seq + ''.join( np.random.choice(np.array(DNA), size = oligo_len - len(fwd_seq), replace = True).tolist() )
      rev_seq = ''.join( np.random.choice(np.array(DNA), size = oligo_len - len(rev_seq), replace = True).tolist() ) + rev_seq
      # naming convention for oligos: assembly id, pool id (number in this assembly), seq id (number in this pool), F/R 
      name_stem = '>' + '|'.join([params['assembly_id'], str(i), str(j)])
      final_oligos.append((name_stem + '|F', fwd_seq))
      final_oligos.append((name_stem + '|R', rev_seq))
  return(final_oligos)

def main(cfg):
  params = dict(cfg.items('Params'))
  random_seed = int(params['random_seed'])
  random.seed(random_seed); np.random.seed(random_seed)
  seqs = list(pandas.read_csv(os.path.expanduser(cfg.get('Files','selected_fn')))['Seqs'])
  assert(all([len(q) == len(seqs[0]) for q in seqs]))
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
  pools, rejected = build_pools(seqs, params)
  final_oligos = fill_oligos(pools, params)
  with open(os.path.expanduser(cfg.get('Files','oligo_fn')),'w') as fo:
    for name, oligo in final_oligos:
      fo.write(name + '\n')
      fo.write(oligo + '\n')
  with open(os.path.expanduser(cfg.get('Files','rejected_fn')),'w') as fr:
    for seq in rejected:
      fr.write(seq + '\n')
  
if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  main(cfg)
  