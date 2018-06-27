import sys
import os
import time
import read_merge_v2 as read_merge
import read_merge_connected_components

import reduce_clusters
import ConfigParser


def timer(text, mark):
  print(text + ': ' + "{:.3g}".format(time.time() - mark))

def main_method(file_in, dir_intermed, cluster_output, file_scores_fwd, file_scores_rev,
    max_gaps, chunk_size, new_cluster_thresh, threads_per_block,
    file_cleaned, file_aln, file_N, file_ambig, len_aln, n_lines, num_bins):
  # generate a group file with fields (read_id, fwd_group, rev_group)
  file_idd = os.path.join(dir_intermed, file_in + '_id') #seq, bin, read_id
  file_orig_groups = os.path.join(dir_intermed, file_in + '_orig_groups') #read_id, fwd_id, rev_id
  merge_start = time.time()
  read_merge.main_method(file_in_dir, file_in, dir_intermed, file_orig_groups, file_idd, file_scores_fwd, file_scores_rev, max_gaps, chunk_size, new_cluster_thresh, threads_per_block)
  timer('Merge runtime', merge_start)
  # use connected components to generate a file whose lines are comma-separated lists of grouped reads
  file_group_list = os.path.join(dir_intermed, file_in + '_group_list')
  # output a final file, with fields (seq, bin, final_group)
  group_start = time.time()
  read_merge_connected_components.main_method(file_idd, file_orig_groups,file_group_list, cluster_output)
  timer('Grouping runtime', group_start)

  reduce_clusters.main_method(cluster_output, file_cleaned, file_aln, file_N, file_ambig, len_aln, n_lines, num_bins)

if __name__=='__main__':
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.read(sys.argv[1])
  file_in = config.get('Input', 'file_in')
  dir_intermed = config.get('Dirs', 'dir_intermed')
  cluster_output = config.get('Files_Intermediate', 'cluster_output')
  file_scores_fwd = config.get('Output', 'file_scores_fwd')
  file_scores_rev = config.get('Output', 'file_scores_rev')
  max_gaps = int(config.get('Params','MAX_GAPS'))
  chunk_size = int(config.get('Params','CHUNK_SIZE'))
  new_cluster_thresh = int(config.get('Params','NEW_CLUSTER_THRESH'))
  threads_per_block = int(config.get('Params','THREADS_PER_BLOCK'))
  file_cleaned = config.get('Files_Intermediate', 'file_cleaned')
  file_aln = config.get('Output', 'file_aligned')
  file_N = config.get('Output', 'file_N')
  file_ambig = config.get('Output', 'file_ambig')
  len_aln = int(config.get('Params','LEN_ALN'))
  #n_lines = int(config.get('Params','N_LINES'))
  n_lines = None
  num_bins = int(config.get('Params','NUM_BINS'))

  main_method(file_in, dir_intermed, cluster_output, file_scores_fwd, file_scores_rev,
    max_gaps, chunk_size, new_cluster_thresh, threads_per_block,
    file_cleaned, file_aln, file_N, file_ambig, len_aln, n_lines, num_bins)
