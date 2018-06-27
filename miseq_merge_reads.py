# NB for NextSeq data, can take a simpler approach and just count sequences without allowing for miscalled bases - 
# those reads are shorter and there are enough that the losses won't be missed.

# Run read_merge.py, read_merge_connected_components.py, and read_merge_final_grouping.py in sequence,
# on a file containing extracted reads from FACS-Seq (in FS6_proc.py, the 'reads_merged_filename')

import sys
import os
import time
import read_merge_v2 as read_merge
import read_merge_connected_components

def timer(text, mark):
  print(text + ': ' + "{:.3g}".format(time.time() - mark))

def main_method(file_in, dir_intermed, file_out, file_scores_fwd, file_scores_rev):
  # generate a group file with fields (read_id, fwd_group, rev_group)
  file_idd = os.path.join(dir_intermed, file_in + '_id') #seq, bin, read_id
  file_orig_groups = os.path.join(dir_intermed, file_in + '_orig_groups') #read_id, fwd_id, rev_id
  merge_start = time.time()
  read_merge.main_method(file_in_dir, file_in, dir_intermed, file_orig_groups, file_idd, file_scores_fwd, file_scores_rev)
  timer('Merge runtime', merge_start)
  # use connected components to generate a file whose lines are comma-separated lists of grouped reads
  file_group_list = os.path.join(dir_intermed, file_in + '_group_list')
  # output a final file, with fields (seq, bin, final_group)
  group_start = time.time()
  read_merge_connected_components.main_method(file_idd, file_orig_groups,file_group_list, file_out)
  timer('Grouping runtime', group_start)

if __name__=='__main__':
  [file_in, dir_intermed, file_out, file_scores_fwd, file_scores_rev] = sys.argv[1:]
  main_method(file_in, dir_intermed, file_out, file_scores_fwd, file_scores_rev)
