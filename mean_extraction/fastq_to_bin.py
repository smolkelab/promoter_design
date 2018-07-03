
# NB: leaving out bio.pairwise.align alignment. Doubt it saves many reads since indels are relatively rare on Illumina platforms.
# nohup time python fastq_to_bin.py FS8_miseq_config.cfg &

import sys, os
import ConfigParser
import imp
import numpy as np
from functools import partial

def progress(count, total, status=''):
  bar_len = 60
  filled_len = int(round(bar_len * count / float(total)))
  percents = round(100.0 * count / float(total), 1)
  bar = '=' * filled_len + '-' * (bar_len - filled_len)
  sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
  sys.stdout.flush() # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)

def main(config_fn):
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.read(config_fn)
  dirs = [os.path.expanduser(config.get('Dirs',q)) for q in ['dir_in', 'dir_left', 'dir_right', 'dir_bc', 'dir_merged']]
  for d in dirs:
    if not os.path.exists(d):
      os.makedirs(d)
  dir_in, dir_left, dir_right, dir_bc, dir_merged = dirs
  fastqs = [q for q in os.listdir(dir_in) if config.get('Params','FASTQ_EXT') in q]
  fastqs.sort()
  fns_list = [tuple([os.path.join(p,q) for p in [dir_in, dir_left, dir_right, dir_bc, dir_merged]]) for q in fastqs]
  compare_mats = [partial(one_fileset, cfg = config)(f) for f in fns_list]

  final_compare_mat = np.sum(np.stack(compare_mats), axis = 0)
  np.savetxt(os.path.expanduser(config.get('Files', 'assign_summary_fn')), final_compare_mat, delimiter = ',')
  merge_all(dir_merged, os.path.expanduser(config.get('Files', 'raw_merged_fn')))
  filter_lengths(os.path.expanduser(config.get('Files', 'raw_merged_fn')), os.path.expanduser(config.get('Files', 'final_merged_fn')), int(config.get('Params','SEQ_LEN')), config.get('Params','LEFT_ONLY'))

def one_fileset(fns, cfg):
  file_in, fail_left, fail_right, fail_bc, merged = fns
  cl = -1
  lines_per_read = int(cfg.get('Params', 'LINES_PER_READ'))
  aligner = inexact_aligner(cfg, file_in)
  with open(file_in, 'r') as fi, open(fail_left, 'w') as fl, open(fail_right, 'w') as fr, open(fail_bc, 'w') as fbc, open(merged, 'w') as fm:
    line_buffer = []
    for line in fi:
      cl += 1
      line_buffer.append(line)
      if len(line_buffer) == lines_per_read:
        output, code = aligner.process_one_read(line_buffer)
        if code == 'L':
          fl.write("%s\n" % output)
        elif code == 'R':
          fr.write("%s\n" % output)
        elif code == 'BC':
          fbc.write("%s\n" % output)
        else:
          fm.write("%s\n" % output)
        line_buffer = []
  return(aligner.assignment_compare_mat)

# merge all read files into a single, alphabetically sorted file
def merge_all(dir_in, fn_out):
  fns = [os.path.join(dir_in, q) for q in os.listdir(dir_in)]
  with open(fn_out, 'w') as fo:
    for f in fns:
      with open(f, 'r') as fi:
        for l in fi:
          fo.write(l)
  # sort using system 'sort' command
  os.system('sort ' + fn_out + ' -o ' + fn_out)

# require a fixed length for all sequences
def filter_lengths(fn_in, fn_out, seq_len, left_only = False):
  if left_only == 'True':
    left_only = True
  else:
    left_only = False
  with open(fn_in, 'r') as fi, open(fn_out, 'w') as fo:
    for l in fi:
      if left_only:
        if len(l.split(',')[0]) >= seq_len:
          l = l.split(',')
          l[0] = l[0][:seq_len]
          l = ','.join(l)
          fo.write(l)
      else:
        if len(l.split(',')[0]) == seq_len:
          fo.write(l)

  # sort, just to be safe
  os.system('sort ' + fn_out + ' -o ' + fn_out)

# class for finding constant regions for one file; see static method 'inexact_align' in earlier versions.
# Implemented as a class to simplify tracking the dictionaries 'fwd_dict', 'rev_dict' caching previous alignment results.
class inexact_aligner(object):
  def __init__(self, cfg, file_in):
    self.cfg = cfg
    self.file_in = file_in.split('/')[-1].split('.')[0]
    self.constants = [self.cfg.get('Params', 'CONSTANT_FWD'), self.cfg.get('Params', 'CONSTANT_REV')]
    self.thresh = int(self.cfg.get('Params', 'CONSTANT_MISMATCH_THRESH'))
    self.depth = int(self.cfg.get('Params', 'CONSTANT_SEARCH_DEPTH'))
    self.left_only = self.cfg.get('Params','LEFT_ONLY') == 'True'
    self.fwd_dict = {}
    self.rev_dict = {}
    self.barcode_fn = cfg.get('Files','barcode_fn')
    assn_module = imp.load_source('assn_module', cfg.get('Files', 'assign_logic_fn'))
    self.assn_obj = assn_module.bin_assigner(self.barcode_fn)
    with open(self.barcode_fn, 'r') as bf:
      self.num_libs = sum(1 for l in bf)
    self.assignment_compare_mat = np.zeros(shape = [self.num_libs, self.num_libs], dtype = 'int')

  def process_one_read(self,read_in):
    seq_raw = read_in[1].strip()
    # Find the left barcode
    seq_raw_L = self.inexact_split(seq_raw, from_left = True)
    if seq_raw_L == None:
      return(seq_raw, 'L')
    [prefix, seq] = seq_raw_L
    phred = read_in[3][:len(prefix)]
    # Find the right barcode, if not LEFT_ONLY
    if self.left_only:
      postfix = None
    else:
      seq_raw_R = self.inexact_split(seq, from_left = False)
      if seq_raw_R == None:
        return(prefix + ',' + phred + ',' + seq, 'R')
      [seq, postfix] = seq_raw_R
    # Get a barcode assignment from file_code, prefix, and postfix.
    # 'assns' will be  list with either one entry (the bin assignment number)
    # or two: the 'official' assignment and a 'check' assignment derived in a different way.
    # None will be returned if an assignment can't be made.
    assns = self.assn_obj.assign_read(prefix, postfix, self.file_in)
    if not None in assns and len(assns) == 2:
      self.assignment_compare_mat[assns[0], assns[1]] += 1
    if None in assns or any([q != assns[0] for q in assns]):
      return(','.join([str(q) for q in assns]) + ',' + seq, 'BC')

    return(seq + ',' + str(assns[0]),'')

  def inexact_split(self, seq, from_left):
    if from_left:
      pattern = self.constants[0]
      search_length = self.depth + len(pattern)
      [seq_to_align, seq_rest] = [seq[:search_length], seq[search_length:]]
      idx = 0
      cache = self.fwd_dict
    else:
      pattern = self.constants[1]
      search_length = self.depth + len(pattern)
      [seq_to_align, seq_rest] = [seq[-search_length:],seq[:-search_length]]
      idx = len(seq_to_align) - len(pattern)
      cache = self.rev_dict

    # try finding a cached alignment
    if seq_to_align in cache:
      split = cache[seq_to_align]

    # try direct alignment
    else:
      split = None
      for i in range(self.depth):
        subseq = seq_to_align[idx:idx+len(pattern)]
        diffs = sum([p != q for (p,q) in zip(subseq, pattern)])
        if diffs < self.thresh:
          split = [seq_to_align[:idx], seq_to_align[idx+len(pattern):]]
          break
        if from_left:
          idx = idx + 1
        else:
          idx = idx - 1
      cache[seq_to_align] = split

    if split == None:
      return(None)

    if from_left:
      return[split[0], split[1] + seq_rest]
    else:
      return[seq_rest + split[0], split[1]]

if __name__ == '__main__':
  main(sys.argv[1])
