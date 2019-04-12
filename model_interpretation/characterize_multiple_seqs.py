# Given a header-less CSV with lines of format (sequence_number,source_file)
# and a target directory,
# for each sequence, create a directory containing the output from characterize_by_mutagenesis.main_single(),
# the output from get_double_mutant_residuals.main() with 'corr' merge,
# and the output from get_double_mutant_residuals.main() with 'l2' merge

import sys
import os
import characterize_by_mutagenesis
import get_double_mutant_residuals

def process_one_sequence(seq_id, source_file, dir_out):
  seq_dir_name = source_file.split('/')[-1]
  seq_dir_name = '.'.join(seq_dir_name.split('.')[:-1]) + '_' + str(seq_id)
  os.mkdir(os.path.join(dir_out, seq_dir_name))
  with open(source_file, 'r') as fi:
    lines = fi.readlines()[1:]
  seq = lines[seq_id].strip().split(',')[0]
  characterize_by_mutagenesis.main_single(seq, source_file, os.path.join(dir_out, seq_dir_name, seq_dir_name + '_single.csv'))
  fn_tmp = os.path.join(dir_out, seq_dir_name, seq_dir_name + '_tmp.csv')
  characterize_by_mutagenesis.main_double(seq, source_file, fn_tmp)
  get_double_mutant_residuals.main(seq, fn_tmp, os.path.join(dir_out, seq_dir_name, seq_dir_name + '_double_corr.csv'), get_double_mutant_residuals.MERGE_DICT[2])
  get_double_mutant_residuals.main(seq, fn_tmp, os.path.join(dir_out, seq_dir_name, seq_dir_name + '_double_l2.csv'), get_double_mutant_residuals.MERGE_DICT[3])
  os.remove(fn_tmp)

def main(table_fn, dir_out):
  with open(table_fn, 'r') as fi:
    for l in fi:
      [source_file, seq_id] = l.strip().split(',')
      source_file = os.path.expanduser(source_file)
      seq_id = int(seq_id)
      process_one_sequence(seq_id, source_file, dir_out)

if __name__ == '__main__':
  [table_fn, dir_out] = [os.path.expanduser(q) for q in sys.argv[1:]]
  main(table_fn, dir_out)
