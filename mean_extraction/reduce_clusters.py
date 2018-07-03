# Do the final reduction of read groups.
# First, exclude singletons.
# For each remaining group, do a muliple sequence alignment; try to generate a consensus sequence.
# Acceptable consensus sequence has a majority-vote base call at each position; N's and gaps don't vote.
# Save reads with all N's at a position in a separate file; maybe some of these can be recovered by imputing constant sequences.

# To do: use the most common group as a test for autocorrelation in the data (there are 88 reads in that group - possible that some reads were miscalled to another group)

# Constant matrices for DNA-to-integer conversion
NTLIST = ["A","C","G","T","N","-"]
NTSTR = ''.join(NTLIST)
NTDICT = dict((j,i) for i,j in enumerate(NTLIST))
LEN_DNA = 4 # there are 4 bases

import sys
import numpy as np
from numba import jit

from time import time
def timer(text, mark):
  print(text + ': ' + "{:.3g}".format(time() - mark))

def count_lines(file_in):
  with open(file_in) as f:
    ans = sum(1 for _ in f)
  return(ans)

# h/t Lukasz Kidzinski
# Function draws a progress bar in consol. It's especially useful for CPU verison which takes hours
# and an empty screen is frightening to people
def progress(count, total, status=''):
  pass
#    bar_len = 60
#    filled_len = int(round(bar_len * count / float(total)))

#    percents = round(100.0 * count / float(total), 1)
#    bar = '=' * filled_len + '-' * (bar_len - filled_len)

#    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
#    sys.stdout.flush() # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)

# Eliminate singletons (and other unclassifiable sequences)
def clean_clusters(file_in, file_out):
  cluster_id = 0
  curr_clust = []
  with open(file_out, 'w') as fo:
    with open(file_in) as fi:
      for line in fi:
        seq, bin_id, gp = line.split(',')
        gp = int(gp)
        if gp != cluster_id: # start a new cluster
          if len(curr_clust) > 1: # eliminate singletons
            for cl in curr_clust:
              fo.write(cl) # cl should still have its trailing newline
          curr_clust = []
          cluster_id = gp
        curr_clust.append(line)
    # write the last cluster
    for cl in curr_clust:
      fo.write(cl)

@jit(nopython = True)
def onehot_one_seq_jit(base_ids, num_bases):
  ans = np.zeros(shape = (num_bases, LEN_DNA), dtype = np.int64)
  for i in range(len(base_ids)):
    basecall = base_ids[i]
    if basecall < LEN_DNA:
      ans[i, basecall] = 1
  return(ans)

def onehot_one_seq(seq, num_bases):
  #ans = np.zeros((num_bases, LEN_DNA), dtype = int)
  base_ids = [NTDICT[b] for b in seq]
  return(onehot_one_seq_jit(base_ids, num_bases))

@jit(nopython = True)
def call_base_from_cts_jit(cts):
  if np.sum(cts) == 0:
    return(4) # 'N' in NTLIST
  if np.sum(cts == np.max(cts)) > 1:
    return(None) # no unique max
  return(np.argmax(cts))

def call_base_from_cts(cts):
  c_idx = call_base_from_cts_jit(cts)
  if c_idx == None:
    return(None)
  return(NTSTR[c_idx]) # maybe *marginally* faster than the list lookup. Probably not.

@jit(nopython = True)
def array_sum(aln_arr):
  ans = np.zeros(shape = aln_arr.shape[1:], dtype = aln_arr.dtype)
  for i in range(ans.shape[0]):
    for j in range(ans.shape[1]):
      ans[i,j] = np.sum(aln_arr[:,i,j])
  return(ans)

# Given an aligned cluster, generate a consensus sequence, and indicate whether it's useable.
def get_consensus_from_alignment(aln):
  # sequences should all be the same length
  if not [len(aln[0])]*len(aln) == [len(q) for q in aln]:
    raise Exception('Sequence lengths in alignment differ: ' + str([len(q) for q in aln]))
  num_bases = len(aln[0])
  has_n = False

  aln_arr = np.zeros(shape = (len(aln), num_bases, LEN_DNA), dtype = int)
  for i in range((len(aln))):
    aln_arr[i,...] = onehot_one_seq(aln[i], num_bases)
  aln_arr = array_sum(aln_arr)

  seq_out = [0]*num_bases
  has_n = False
  for i in range(num_bases):
    call = call_base_from_cts(aln_arr[i,:])
    if call == 'N':
      has_n = True
    if call == None:
      return(None, None)
    seq_out[i] = call

  return(''.join(seq_out), has_n)

def collect_seqs_bins(reads, num_bins):
  seqs_in = []
  bins_collected = np.zeros(num_bins, dtype = int)
  for read in reads:
    s, b = read.split(',')
    seqs_in.append(s)
    bins_collected[int(b)] += 1
  return(seqs_in, bins_collected)

def lazy_alignment(reads, len_out, num_bins, debug = False):
  if debug:
    print(reads)
    print(len_out)
  seqs_in, bins_collected = collect_seqs_bins(reads, num_bins)
  if debug:
    print(seqs_in)
    print(bins_collected)

  seqs_in = [q[:len_out] for q in seqs_in]
  seqs_in = [q + '-'*(len_out - len(q)) for q in seqs_in]
  if debug:
    print(seqs_in)
  aln, has_n = get_consensus_from_alignment(seqs_in)
  if debug:
    print(aln)
  return(aln, bins_collected, has_n)

# given an open file object and the first line of a cluster, read in the next cluster. Return the cluster and the first line of the next cluster (or None if the file ended).
# NB: file is not numerically sorted by cluster_id, but clusters are grouped together.
def read_next_cluster(file_in, cluster_start):
  clust_id = int(cluster_start.split(',')[2])
  clust_out = [cluster_start]
  while(True):
    try:
      l = file_in.next()
      n_id = int(l.split(',')[2])
    except StopIteration:
      l = None
      n_id = None
    if n_id == clust_id:
      clust_out.append(l)
    else:
      start_out = l
      break
  return([clust_out, start_out])

def write_cluster(clust, fa, fN, fam, len_aln, num_bins):
  clust_id = clust[0].split(',')[2].strip()
  debug = False #clust_id == '1000013'
  aln, bins_collected, has_n = lazy_alignment([','.join(q.split(',')[:2]) for q in clust], len_aln, num_bins, debug) # aln is [None, bins_collected] or [seq_out, has_n, bins_collected]
  c_class = 'good'
  if aln == None: # no real alignment
    c_class = 'ambig'
  elif has_n: # alignment has N's
    c_class = 'N'
  if c_class == 'good':
    bins = [str(q) for q in bins_collected.tolist()]
    seq_out = [aln, str(clust_id)]
    seq_out.extend(bins)
    seq_out = ','.join(seq_out)
    fa.write("%s\n" % seq_out)
  elif c_class == 'ambig':
    for i in clust:
      fam.write(i)
  else:
    for i in clust:
      fN.write(i)

# make sure all the reads in a cluster are actually from the same group
def test_cluster(clust):
  clust_ids = [int(q.split(',')[2]) for q in clust]
  if not all([q == clust_ids[0] for q in clust_ids]):
    print(clust)
    raise Exception('Cluster extraction is bugged - cluster IDs inconsistent')


# Try to use 'lazy_alignment' to align clusters.
def try_lazy_alignment(file_in, file_aln, file_N, file_ambig, len_aln, num_bins, n_lines = None):
  with open(file_aln, 'w') as fa, open(file_N, 'w') as fN, open(file_ambig, 'w') as fam:
    with open(file_in) as fi:
      start_next = fi.next()
      lines_done = 0
      while(start_next != None):
        clust, start_next = read_next_cluster(fi, start_next)
        write_cluster(clust, fa, fN, fam, len_aln, num_bins)
        lines_done += len(clust)
        if n_lines != None:
          progress(lines_done, int(n_lines))

def main_method(file_in, file_cleaned, file_aln, file_N, file_ambig, len_aln, n_lines, num_bins):
  t = time()
  clean_clusters(file_in, file_cleaned)
  timer('Clusters cleaned', t)
  t = time()
  try_lazy_alignment(file_cleaned, file_aln, file_N, file_ambig, len_aln, num_bins, n_lines)
  timer('Lazy alignment done', t)

if __name__ == '__main__':
  file_in, file_cleaned, file_aln, file_N, file_ambig, len_aln, num_bins = sys.argv[1:]
  len_aln = int(len_aln)
  num_bins = int(num_bins)
  n_lines = count_lines(file_in)
  main_method(file_in, file_cleaned, file_aln, file_N, file_ambig, len_aln, n_lines, num_bins)
