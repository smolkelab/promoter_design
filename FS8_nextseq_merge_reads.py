
# NB: now this is 0-indexed, and reading from one giant file. Different from FS7.
import sys
import os
import numpy as np
import pandas

class read_tracker(object):

  def __init__(self, n_bins):
    self.keys = set()
    self.vals = {}
    self.n_bins = n_bins

  def add_line(self, line):
    seq, bin_id = line.strip().split(',')
    bin_id = int(bin_id)
    assert(bin_id >= 0)
    assert(bin_id < self.n_bins)
    # create a new seq entry, if needed
    if not seq in self.keys:
      x = np.zeros(shape = (self.n_bins,), dtype = 'int32')
      self.keys.add(seq)
      self.vals[seq] = x

    # increment the appropriate bin
    self.vals[seq][bin_id] += 1

  def summarize(self):
    ans_seqs = list(self.keys)
    ans_seqs.sort()
    ans_cts = np.zeros(shape = (len(ans_seqs), self.n_bins), dtype = 'int32')
    for i,k in enumerate(ans_seqs):
      ans_cts[i,:] = self.vals[k]
    return(ans_seqs, ans_cts)

def process_one_file(fn, n_bins):
  print('Reading file: ' + fn)
  with open(fn, 'r') as fi:
    lines = [q for q in fi] # any faster to read in a whole file at once?

  print('Building tracker: ' + fn)
  reads = read_tracker(n_bins)
  for l in lines:
    reads.add_line(l)
  return(reads)

if __name__ == '__main__':
  print('Scanning read file\n')
  in_file = sys.argv[1]
  n_bins = int(sys.argv[2])
  out_file = sys.argv[3]
  reads = process_one_file(in_file, n_bins)

  print('Extracting output\n')
  ans_seqs, ans_cts = reads.summarize()
  print('Seqs: ' + str(len(ans_seqs)))
  dict_out = {}
  dict_out['Seq'] = ans_seqs
  for i in range(ans_cts.shape[1]):
    dict_out[str(i)] = ans_cts[:,i]
  df_out = pandas.DataFrame(dict_out)
  df_out.to_csv(out_file)
