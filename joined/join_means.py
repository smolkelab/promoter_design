import sys
import os
import ConfigParser

def process_one_dataset(t):
  mean_fn, fwd_pad, rev_pad = t
  header_done = False
  ans = []
  with open(mean_fn, 'r') as fi:
    for l in fi:
      if not header_done: # skip the first line
        header_done = True
      else:
        l = l.strip().split(',')
        assert (len(l) == 2 or len(l) == 3)
        l[0] = fwd_pad + l[0] + rev_pad
        if len(l) == 2:
          l.append(l[1])
        ans.append(l)
  return(ans)

def main_method(config):
  mean_fns = [os.path.expanduser(q) for q in config.get('Files','mean_fns').strip().split(',')]
  fwd_pads = [q for q in config.get('Params','fwd_pad').strip().split(',')]
  rev_pads = [q for q in config.get('Params','rev_pad').strip().split(',')]
  fn_out = os.path.expanduser(config.get('Files','fn_out'))
  dataset_params = zip(mean_fns, fwd_pads, rev_pads)
  with open(fn_out, 'w') as fo:
    fo.write('Seq,Strength_A,Strength_B\n')
    for i in dataset_params:
      output = process_one_dataset(i)
      for l in output:
        fo.write(','.join(l) + '\n')

if __name__ == '__main__':
  config = ConfigParser.RawConfigParser()
  config.read(sys.argv[1])
  main_method(config)
