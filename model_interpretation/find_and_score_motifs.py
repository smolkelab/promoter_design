
# Given a sequence, find the lowest-scoring motif of a given length, and find out what the motif is
# define position as the first base of the motif
# define strength as the mean loss in strength due to a mutation in the motif
# if there are multiple TATA boxes, return the strongest
import sys
import os
import numpy as np
import pandas as pd

DNA = ['A','C','G','T']
DNA_POS = {q:i for (i,q) in enumerate(DNA)}
BIG_POS = 100000.

# IUPAC standard base codes: sometimes multiple positions equal 0.
IUPAC = {i:q for (i,q) in enumerate(['-','T','G','K','C','Y','S','B','A','W','R','D','M','H','V','N'])}

# 'seq_mat' is a score matrix; row for each position, column for each base (A,C,G,T)
# for each possible window, find if the subsequence matches the motif; if yes, get the mean score diff of all mutations
# return the index and score of the biggest motif

def mat_to_score(mat):
  muts = np.sum(mat)
  num_muts = float(mat.shape[0]*(mat.shape[1] - 1))
  return muts/num_muts

'''def matches_seq(seq_mat, seq):
  if seq_mat.shape[0] != len(seq):
    return False
  for (row, base) in zip(seq_mat, seq):
    base_id = DNA_POS[base]
    if np.any(row[base_id] != 0.):
      return False
  print seq_mat
  print seq
  return True'''

def get_seq(mat):
  ans = []
  pow_arr = [2**q for q in range(mat.shape[1]-1,-1,-1)]
  for x in mat:
    x = x == 0.
    x = int(np.sum(x*pow_arr))
    x = IUPAC[x]
    ans.append(x)
  return(x)

def find_strength_pos_seq(seq_mat, motif_len):
  min_score = BIG_POS
  min_pos = -1
  min_seq = ''
  for i in range(seq_mat.shape[0] - motif_len + 1):
    window = seq_mat[i:i+motif_len,]
    score = mat_to_score(window)
    if score < min_score:
      score = min_score
      min_pos = i
      min_seq = get_seq(window)

  return(min_score, min_pos, min_seq)

def process_one_file(fn_in, fn_out, dir_single_muts, motif_len, min_pos, max_pos):
  dat_in = pd.read_csv(fn_in)
  fn_stem = fn_in.split('/')[-1]
  fn_stem = '.'.join(fn_stem.split('.')[:-1])
  scores = []
  poses = []
  motifs = []
  for i in range(dat_in.shape[0]):
    fn_muts = os.path.join(dir_single_muts, fn_stem + '_' + str(i) + '_single.csv')
    mat_muts = np.loadtxt(fn_muts, delimiter = ',')[min_pos:max_pos,]
    (score, pos, seq) = find_strength_pos_seq(mat_muts, motif_len)
    # express 'pos' in terms of the original sequence
    pos = pos + min_pos
    scores.append(score); poses.append(pos); motifs.append(seq)
  dat_in['Mut_score'] = scores
  dat_in['Pos'] = poses
  dat_in['Motif'] = motifs
  dat_in.to_csv(fn_out)

def main(fn_csv):
  with open(fn_csv, 'r') as fi:
    for l in fi:
      [fn_in, fn_out, dir_single_muts, motif_len, min_pos, max_pos] = l.strip().split(',')
      [fn_in, fn_out, dir_single_muts] = [os.path.expanduser(q) for q in [fn_in, fn_out, dir_single_muts]]
      min_pos = int(min_pos); max_pos = int(max_pos); motif_len = int(motif_len)
      process_one_file(fn_in, fn_out, dir_single_muts, motif_len, min_pos, max_pos)

if __name__ == '__main__':
  main(sys.argv[1])
