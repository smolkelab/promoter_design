# Figure SH: need to show that read clustering works and eliminates related sequence problem

# focus on 'aligned' clusters (successfully built, etc.)
# specifically, on the larger ones (more likely to have daughters somewhere?)
# get representative intra-cluster distance
# get representative consensus-to-consensus distance: show that nothing is nearby!!!
# there's roughly 1M clusters
# compare just spacers? e.g. upstream to all upstreams, 5' UTR to all 5' UTRs?
# sequence-to-sequence Hamming distance
# spacer-to-spacer gapped alignment? would need to be an in- and a del- in the same spacer... should be super rare

import os
import ConfigParser
import hgbrian_GPU_NW

FN_TMP_SEQS = os.path.expanduser('~/facs-seq_test/seqs_tmp.txt')
FN_TMP_SCORES = os.path.expanduser('~/facs-seq_test/scores_tmp.txt')
MAX_GAPS = 5
CHUNK_SIZE = 50000
THREADS_PER_BLOCK = 64


# given a file of \n-delimited seqs and a reference sequence seq_in (a string),
# get a list where entry i is the NW distance between seq_in and line i in fn_in
def file_vs_ref_dist(fn_in, seq_in, fn_tmp_seqs = FN_TMP_SEQS, fn_tmp_scores = FN_TMP_SCORES, max_gaps = MAX_GAPS, chunk_size = CHUNK_SIZE, threads_per_block = THREADS_PER_BLOCK):
  with open(fn_in, 'r') as fi, open(fn_tmp_seqs, 'w') as ft:
    for l in fi:
      ft.write(l)
      ft.write(seq_in + '\n')

  hgbrian_GPU_NW.score_file(fn_tmp_seqs, fn_tmp_scores, max_gaps, chunk_size, threads_per_block)
  ans = []
  with open(fn_tmp_scores, 'r') as fi: # record every other line, starting with the first (line 0)
    for (i,l) in enumerate(fi):
      if i % 2 == 0:
        ans.append(l)
  os.remove(fn_tmp_seqs)
  os.remove(fn_tmp_scores)
  return(ans)

def list_vs_ref_dist(seqs_target_list, seq_in, fn_in_tmp):
  with open(fn_in_tmp, 'w') as fo:
    for s in seqs_target_list:
      fo.write(s + '\n')
  ans = file_vs_ref_dist(fn_in_tmp, seq_in)
  os.remove(fn_in_tmp)
  return(ans)

def list_vs_self_dist(seqs_target_list, fn_in_tmp):
  ans = {}
  for (i,s) in enumerate(seqs_target_list):
    dists = list_vs_ref_dist(seqs_target_list[i+1:], s, fn_in_tmp)
    ans[s] = dists
  return(ans)
  
# given a CSV with header Group,Pre,Ambig,N,Aln, get the group IDs with 'Aln' >= some threshold
def get_biggest_groups(fn_in, aln_thresh):
  ans = []
  with open(fn_in, 'r') as fi:
    is_header = True
    for l in fi:
      if is_header:
        is_header = False
      else:
        [gp, _, _, _, aln] = [int(q) for q in l.strip().split(',')]
        if aln >= aln_thresh:
          ans.append(gp)
  return(ans)

# given a list of reads, get the longest distance within the reads
def longest_dist_within_group(reads, fn_in_tmp):
  '''reads = []
  with open(clustered_read_fn, 'r') as fi:
    for l in fi:
      [seq, _, this_gp_id] = l.strip().split(',')
      if int(this_gp_id) == target_gp_id:
        reads.append(seq)'''
  dist_dict = list_vs_self_dist(reads, fn_in_tmp)
  max_dist = -1
  for k in dist_dict.keys():
    tmp = max(dist_dict[k])
    if tmp > max_dist:
      max_dist = tmp
  return(max_dist)

# given a seq and the filename of a read table, get the minimum distance to another group, and a list of all group IDs with that distance
# read table has no header; first CSV field is sequence, next is group ID
def min_dist_in_read_table_to_seq(seq, read_table_fn, fn_in_tmp):
  seqs = []
  gp_ids = []
  with open(read_table_fn, 'r') as fi:
    for l in fi:
      [s, g] = l.strip().split(',')[:2]
      g = int(g)
      if s != seq:
        seqs.append(s)
        gp_ids.append(g)
  
  # get the distance to each line in the read table
  dists = list_vs_ref_dist(seqs, seq, fn_in_tmp)
  target_dist = min(dists)
  
  # get all group IDs with this distance  
  gps_out = []
  for (g,d) in zip(gp_ids, dists):
    if d == target_dist:
      gps_out.append(g)

  return(target_dist, gps_out)

# given a list of group IDs and a read table, get a list of the sequence corresponding to each group ID
def get_seqs_from_group_ids(group_ids, read_table_fn):
  ans_dict = {}
  for g in group_ids:
    ans_dict[g] = None
  with open(read_table_fn, 'r') as fi:
    for l in fi:
      [s, g] = l.strip().split(',')[:2]
      if int(g) in ans.keys():
        ans_dict[g] = s
  ans = [ans_dict[q] for q in group_ids]
  return(ans)
  
# given a config object, get the following:
# a table gp_id, longest internal distance for all groups above a certain number of reads
# a 'table' gp_id, # reads, minimum distance to another group, comma-separated list of gp_ids with that distance, for all groups above a certain number of reads

def main_method(cfg, reads_thresh, internal_dist_fn_out, intergroup_dist_fn_out, fn_in_tmp = FN_TMP_SEQS):
  target_groups = get_biggest_groups(fn_in = os.path.expanduser(cfg.get('Output', 'file_fates')), aln_thresh = reads_thresh)
  target_seqs = get_seqs_from_group_ids(group_ids = target_groups, read_table_fn = os.path.expanduser(cfg.get('Output', 'file_aligned_filtered')))
  # drop groups without a seq found (meaning seq was dropped in filtering)
  final_groups = []
  final_seqs = []
  for (g,s) in zip(target_groups, target_seqs):
    if s is not None:
      final_groups.append(g)
      final_seqs.append(s)

  # get a dictionary group_id:reads
  reads_clustered_fn = os.path.expanduser(cfg.get('Files_Intermediate', 'cluster_output'))
  reads_dict = {}
  for g in final_groups:
    reads_dict[g] = []
  with open(reads_clustered_fn, 'r') as fi:
    for l in fi:
      [seq, _, g] = l.strip().split(',')
      g = int(g)
      if g in reads_dict.keys():
        reads_dict[g].append(seq)
  
  # write a table gp_id,longest internal distance to internal_dist_fn_out
  # modify longest_dist_within_group!
  with open(internal_dist_fn_out, 'w') as fo:
    for g in final_groups:
      d = longest_dist_within_group(reads_dict[g], fn_in_tmp)
      fo.write(str(g) + ',' + str(d) + '\n')
  
  # write a table gp_id, shortest distance to another seq, CSV list of group ids with that distance
  with open(intergroup_dist_fn_out, 'w') as fo:
    for (g,s) in zip(final_groups, final_seqs):
      target_dist, gps_out = min_dist_in_read_table_to_seq(s, os.path.expanduser(cfg.get('Output', 'file_aligned_filtered')), fn_in_tmp)
    ans = str(g) + ',' + str(target_dist) + ',' + ','.join(gps_out)
    fo.write(ans + '\n')

if __name__ == '__main__': # debug: change this!
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(os.path.expanduser('~/facs-seq/GPD/miseq/build_table_align_GPD.cfg'))
  main_method(cfg, reads_thresh = 50, 
    internal_dist_fn_out = os.path.expanduser('~/facs-seq_test/internal_test.csv'), 
    intergroup_dist_fn_out = os.path.expanduser('~/facs-seq_test/intergroup_test.csv'))