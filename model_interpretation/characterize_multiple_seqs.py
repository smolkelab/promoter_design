# Given a header-less CSV with lines of format (sequence_number,source_file)
# and a target directory,
# for each sequence, create a directory containing the output from characterize_by_mutagenesis.main_single(),
# the output from get_double_mutant_residuals.main() with 'corr' merge,
# and the output from get_double_mutant_residuals.main() with 'l2' merge

import sys
import os
import characterize_by_mutagenesis
import get_double_mutant_residuals

def process_one_sequence(seq_id, source_file, dir_out, loaded_models = None, do_double = True):
  seq_dir_name = source_file.split('/')[-1]
  seq_dir_name = '.'.join(seq_dir_name.split('.')[:-1]) + '_' + str(seq_id)
  os.mkdir(os.path.join(dir_out, seq_dir_name))
  with open(source_file, 'r') as fi:
    lines = fi.readlines()[1:]
  seq = lines[seq_id].strip().split(',')[0]
  loaded_models = characterize_by_mutagenesis.main_single(seq, source_file, os.path.join(dir_out, seq_dir_name, seq_dir_name + '_single.csv'), loaded_models = loaded_models)
  if do_double:
    fn_tmp = os.path.join(dir_out, seq_dir_name, seq_dir_name + '_tmp.csv')
    loaded_models = characterize_by_mutagenesis.main_double(seq, source_file, fn_tmp, loaded_models = loaded_models)
    get_double_mutant_residuals.main(seq, fn_tmp, os.path.join(dir_out, seq_dir_name, seq_dir_name + '_double_corr.csv'), get_double_mutant_residuals.MERGE_DICT[2])
    get_double_mutant_residuals.main(seq, fn_tmp, os.path.join(dir_out, seq_dir_name, seq_dir_name + '_double_l2.csv'), get_double_mutant_residuals.MERGE_DICT[3])
  os.remove(fn_tmp)
  return loaded_models

def main(table_fn, dir_out):
  # prevent having to reload models every time
  loaded_models = None
  with open(table_fn, 'r') as fi:
    for l in fi:
      tmp = l.strip().split(',')
      if len(tmp) == 2:
        [source_file, seq_id] = tmp
        source_file = os.path.expanduser(source_file)
        seq_id = int(seq_id)
        loaded_models = process_one_sequence(seq_id, source_file, dir_out, loaded_models = loaded_models)
      if len(tmp) == 1:
        loaded_models = main_fullfile_onlysingle(tmp[0], dir_out, loaded_models = loaded_models)

def main_fullfile_onlysingle(fn_in, dir_out, loaded_models = None):
  source_file = os.path.expanduser(fn_in)
  with open(source_file, 'r') as fi:
    lines = fi.readlines()[1:]
  num_seqs = len(lines)
  for i in range(num_seqs):
    loaded_models = process_one_sequence(i, source_file, dir_out, loaded_models = loaded_models, do_double = False)
  return loaded_models

if __name__ == '__main__':
  [table_fn, dir_out] = [os.path.expanduser(q) for q in sys.argv[1:]]
  main(table_fn, dir_out)