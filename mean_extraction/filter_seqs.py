import sys
import os
import ConfigParser

ALPHA = ['A','C','G','T']

# given a config file, find the read table to be filtered and the params 'seq', 
# as well as params explaining what each base in 'seq' encodes.
# For each line in the read file, copy it to the final output if it is the same length as 'seq',
# and it matches the expression.
def main_method(cfg):
  fn_in = os.path.expanduser(cfg.get('Output','file_aligned'))
  fn_out = os.path.expanduser(cfg.get('Output','file_aligned_filtered'))
  seq = cfg.get('Params','SEQ')
  # determine what bases are allowed at each position in 'seq'
  b_allowed = []
  for b in seq:
    if b in ALPHA:
      b_allowed.append([b])
    else:
      b_code = [int(q) for q in cfg.get('Params',b).split(':')]
      assert(len(b_code) == len(ALPHA))
      b_allowed.append([q for (p,q) in zip(b_code, ALPHA) if p > 0]) # if a base can possibly appear at this position, add it to the list of allowed bases
  with open(fn_in, 'r') as fi, open(fn_out,'w') as fo:
    for l in fi:
      s = l.split(',')[0]
      writeable = all([p in q for (p,q) in zip(s, b_allowed)])
      if writeable:
        fo.write(l)

if __name__ == '__main__':
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.read(sys.argv[1])
  main_method(config)
