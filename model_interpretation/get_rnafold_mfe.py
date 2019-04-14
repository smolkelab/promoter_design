# Use RNAfold to get folding energies of sequence 5' UTRs

import sys
import os
import ConfigParser
import tempfile
import subprocess

# input file will be a CSV with header Seqs,Scores
# write a temporary corresponding .fa file
def csv_to_fasta(fn_in, dir_in, dir_tmp, utr_bases, context):
  fn_out = fn_in.split('.')[0] + '.fa'
  idx = 0
  first = True
  with open(os.path.join(dir_in, fn_in), 'r') as fi, open(os.path.join(dir_tmp, fn_out), 'w') as fo:
    for l in fi:
      # skip the header
      if first:
        first = False
      else:
        seq = l.strip().split(',')[0]
        seq = seq[utr_bases[0]:utr_bases[1]]
        seq = context[0] + seq + context[1]
        fo.write('>' + str(idx) + '\n')
        fo.write(seq + '\n')
        idx += 1
  return fn_out

def fold_file(dir_in, fn_in, dir_out):
  fn_out = fn_in.split('.')[0] + '.txt'
  with open(os.path.join(dir_out, fn_out), 'w') as fo:
    p = subprocess.Popen(['RNAfold', os.path.join(dir_in, fn_in)], stdout=fo)
    p.wait()
  return fn_out

def read_mfes(dir_in, fn_in):
  fn = os.path.join(dir_in, fn_in)
  # fn will contain three-line records: two lines of FASTA input and one of output, formatted as 'structure (mfe)'
  ans = []
  with open(fn, 'r') as fi:
    for (i,l) in enumerate(fi):
      if i % 3 == 2:
        l = l.strip().split(' ')
        l = ' '.join(l[1:])
        #try:
        #  assert(len(l) == 2)
        #except AssertionError:
        #  print(l)
        #  print(len(l))
        #  raise Exception('RNAfold output not parsed correctly')
        mfe = float(l[1:-1].strip())
        ans.append(mfe)

  return ans

def process_one_file(fn_in, dir_in, dir_tmp, dir_out, fn_out, utr_bases, context):
  # get the folding energies
  fn_tmp_fa = csv_to_fasta(fn_in, dir_in, dir_tmp, utr_bases, context)
  rnafold_out = fold_file(dir_tmp, fn_tmp_fa, dir_tmp)
  mfes = read_mfes(dir_tmp, rnafold_out)

  # remove temporary files
  os.remove(os.path.join(dir_tmp, fn_tmp_fa))
  os.remove(os.path.join(dir_tmp, rnafold_out))

  # write MFEs and original file content to output file
  first = True
  with open(os.path.join(dir_in, fn_in), 'r') as fi, open(os.path.join(dir_out, fn_out), 'w') as fo:
    for l in fi:
      if first:
        first = False
        fo.write(l.strip() + ',MFE\n')
      else:
        fo.write(l.strip() + ',' + str(mfes.pop(0)) + '\n')

def main(cfg):
  dir_tmp = tempfile.mkdtemp()
  #try:
  #dir_in = os.path.expanduser(cfg.get('Dirs', 'dir_in'))
  dir_out = os.path.expanduser(cfg.get('Dirs', 'dir_out'))
  context = cfg.get('Params', 'context').split(',')
  fns_table = cfg.get('Params', 'fns_table')
  with open(fns_table, 'r') as fi:
    for l in fi:
      [dir_in, fn_in, fn_out, utr_start, utr_end] = l.strip().split(',')
      dir_in = os.path.expanduser(dir_in)
      utr_bases = [int(utr_start), int(utr_end)]
      process_one_file(fn_in, dir_in, dir_tmp, dir_out, fn_out, utr_bases, context)

  #finally:
  os.rmdir(dir_tmp)

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  main(cfg)
