import sys
import os
import time
import read_merge_v2 as read_merge
import read_merge_connected_components
import reduce_clusters
import filter_seqs
import ConfigParser

def timer(text, mark):
  print(text + ': ' + "{:.3g}".format(time.time() - mark))

def main_method(file_in, dir_intermed, cluster_output, file_scores_fwd, file_scores_rev,
    max_gaps, chunk_size, new_cluster_thresh, threads_per_block,
    file_cleaned, file_aln, file_N, file_ambig, len_aln, n_lines, num_bins, **kwargs):
  # generate a group file with fields (read_id, fwd_group, rev_group)
  file_in_stub = file_in.split('/')[-1]
  file_idd = os.path.join(dir_intermed, file_in_stub + '_id') #seq, bin, read_id
  file_orig_groups = os.path.join(dir_intermed, file_in_stub + '_orig_groups') #read_id, fwd_id, rev_id
  merge_start = time.time()
  read_merge.main_method(file_in, dir_intermed, file_orig_groups, file_idd, file_scores_fwd, file_scores_rev, max_gaps, chunk_size, new_cluster_thresh, threads_per_block)
  timer('Merge runtime', merge_start)
  # use connected components to generate a file whose lines are comma-separated lists of grouped reads
  file_group_list = os.path.join(dir_intermed, file_in_stub + '_group_list')
  # output a final file, with fields (seq, bin, final_group)
  group_start = time.time()
  read_merge_connected_components.main_method(file_idd, file_orig_groups,file_group_list, cluster_output)
  timer('Grouping runtime', group_start)

  reduce_clusters.main_method(cluster_output, file_cleaned, file_aln, file_N, file_ambig, len_aln, n_lines, num_bins)

  filter_seqs.main_method(kwargs['cfg'])

if __name__=='__main__':
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.read(sys.argv[1])
  file_in = os.path.expanduser(config.get('Input', 'file_in'))
  dir_intermed = os.path.expanduser(config.get('Dirs', 'dir_intermed'))
  if not os.path.exists(dir_intermed):
    os.makedirs(dir_intermed)

  cluster_output = os.path.expanduser(config.get('Files_Intermediate', 'cluster_output'))
  file_scores_fwd = os.path.expanduser(config.get('Output', 'file_scores_fwd'))
  file_scores_rev = os.path.expanduser(config.get('Output', 'file_scores_rev'))
  max_gaps = int(config.get('Params','MAX_GAPS'))
  chunk_size = int(config.get('Params','CHUNK_SIZE'))
  new_cluster_thresh = int(config.get('Params','NEW_CLUSTER_THRESH'))
  threads_per_block = int(config.get('Params','THREADS_PER_BLOCK'))
  file_cleaned = os.path.expanduser(config.get('Files_Intermediate', 'file_cleaned'))
  file_aln = os.path.expanduser(config.get('Output', 'file_aligned'))
  file_aln_filtered = os.path.expanduser(config.get('Output','file_aligned_filtered'))
  file_N = os.path.expanduser(config.get('Output', 'file_N'))
  file_ambig = os.path.expanduser(config.get('Output', 'file_ambig'))
  #len_aln = int(config.get('Params','LEN_ALN'))
  len_aln = len(config.get('Params','SEQ').strip())
  n_lines = None
  num_bins = int(config.get('Params','NUM_BINS'))

  main_method(file_in, dir_intermed, cluster_output, file_scores_fwd, file_scores_rev,
    max_gaps, chunk_size, new_cluster_thresh, threads_per_block,
    file_cleaned, file_aln, file_N, file_ambig, len_aln, n_lines, num_bins, cfg = config)
