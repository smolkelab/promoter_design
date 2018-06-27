# Original brittle read merge tool for FS6_proc.py lost almost all of my reads - 
# from 7.7M down to ~60K!
# Visual inspection of files suggests more reads can be recovered.
# Approach: start with file of reads with 5' and 3' constants stripped and bin ID'd,
# sorted alphabetically by sequence.
# Filter by sequence length to sub-files (within a specified range);
# will not try to deal with indels.
# First, use sort to cluster related sequences together by allowing a specified # of mismatches,
# with 'N' treated as a wildcard.
# Then reverse sequences and capture any reads that were mistakenly lost due to mutations in the 5' region.
# Generate a consensus sequence for each read set, or reject them (with a record.)
# Do the final merge.

import sys
import os

from time import time
def timer(text, mark):
  print(text + ': ' + "{:.3g}".format(time() - mark))

import hgbrian_GPU_NW

MAX_GAPS = 5
CHUNK_SIZE = 50000
NEW_CLUSTER_THRESH = 197 # see ~/Dropbox/Lab/Data/2018-01 for a histogram of NW selection thresholds with settings as in hgbrian_GPU_NW_score_test.py

def affix_read_id(file_in, file_out):
  read_id = 0
  with open(file_in) as fi, open(file_out,'w') as fo:
    for line in fi:
      line_out = line.strip() + ',' + str(read_id)
      fo.write("%s\n" % line_out)
      read_id += 1

# Associate each read with an ID number and a cluster number. (NB: there may be other fields already - preserve these, append new fields to the end of the line.)
# But assume CSV lines.
# Compare each read to the immediately following read.
def group_adjacent_reads(file_in, file_out, file_tmp, has_seq_and_bin, score_file):
  hgbrian_GPU_NW.score_file(file_in, file_tmp, max_gaps = MAX_GAPS, chunk_size = CHUNK_SIZE)

  cluster_id = 0
  #read_id = 0
  last_line_found = False
  with open(file_in, 'r') as fi, open(file_tmp, 'r') as fs, open(file_out, 'w') as fo, open(score_file, 'w') as f_scores:
    for line in fi:
      # drop the sequence and bin, if needed; append the cluster ID
      line_split = line.strip().split(',')
      if has_seq_and_bin:
        line_split = line_split[2:]
      #line_split.extend([str(read_id), str(cluster_id)])
      line_split.extend([str(cluster_id)])
      line_out = ','.join(line_split)
      fo.write("%s\n" % line_out)
      # difference between "this" line and its *successor*
      # will not exist for the last line, which has no successor,
      # so allow a 'StopIteration' exactly once.
      try:
        score = float(fs.next())
        f_scores.write(str(score) + '\n')
        if score < NEW_CLUSTER_THRESH:
          cluster_id += 1 # <- update the cluster_id of the next read, if needed
      except StopIteration: # originally just had a 'pass', which allowed too-short distance files to be accepted.
        if last_line_found:
          raise StopIteration
        else:
          last_line_found = True
  #os.remove(file_tmp)

def reverse_and_resort(file_in, file_out):
  with open(file_in, 'r') as fi:
    with open(file_out, 'w') as fo:
      for line in fi:
        line_split = line.split(',')
        line_split[0] = line_split[0][::-1] # <-- Reversal: https://stackoverflow.com/questions/931092/reverse-a-string-in-python
        line_out = ','.join(line_split)
        fo.write(line_out)
  # Sort the new file:
  os.system('sort -s -t , -k1,1 ' + file_out + ' -o ' + file_out)

# These files should each have a read_id and cluster_id; output one file in read_id order, with both clusters for each read.
def merge_grouped_reads(file_1, file_2, file_out):
  # sort the files by the ID
  os.system('sort -s -t , -k1,1 ' + file_1 + ' -o ' + file_1)
  os.system('sort -s -t , -k1,1 ' + file_2 + ' -o ' + file_2)

  with open(file_1, 'r') as f1, open(file_2, 'r') as f2, open(file_out, 'w') as fo:
    for line_pair in zip(f1,f2):
      lines_split = [q.strip().split(',') for q in line_pair]
      reads = [q[0] for q in lines_split]
      if reads[0] != reads[1]:
        raise Exception('Mismatch while merging grouped reads' + str(reads[0]) + ' ' + str(reads[1]) )
      x = [reads[0]]
      x.extend([q[1] for q in lines_split])
      line_out = ','.join(x)
      fo.write("%s\n" % line_out)

def main_method(file_in, dir_intermed, file_out, file_idd, file_scores_fwd, file_scores_rev):
  # affix read IDs; file_in has fields (seq, bin)
  # file_idd is provided, as it'll be needed later
  t = time()
  affix_read_id(file_in, file_idd)
  timer('Affix read ID',t)

  # get the forward grouping
  t = time()
  file_tmp = os.path.join(dir_intermed, file_in + '_fwd_group_tmp')
  file_fwd = os.path.join(dir_intermed, file_in + '_fwd_group')
  group_adjacent_reads(file_idd, file_fwd, file_tmp, True, file_scores_fwd) # read_id, fwd_cluster_id
  timer('Forward grouping',t)

  # reverse and resort the original file
  t = time()
  file_reversed = os.path.join(dir_intermed, file_in + '_id_reversed')
  reverse_and_resort(file_idd, file_reversed) # reversed_seq, bin, read_id
  timer('Reverse and resort',t)

  # get the reverse grouping
  t = time()
  file_tmp = os.path.join(dir_intermed, file_in + '_rev_group_tmp')
  file_rev = os.path.join(dir_intermed, file_in + '_rev_group')
  group_adjacent_reads(file_reversed, file_rev, file_tmp, True, file_scores_rev) # read_id, rev_cluster_id
  timer('Reverse grouping',t)

  # merge the group files; file_out has fields (read_id, fwd_group, rev_group)
  t = time()
  merge_grouped_reads(file_fwd, file_rev, file_out)
  timer('Resort forward and reverse groups',t)

if __name__ == '__main__':
  [file_in, dir_intermed, file_out] = sys.argv[1:]
  main_method(file_in, dir_intermed, file_out)
