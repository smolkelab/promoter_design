# For each line in a header-less input file with sequences in column 0,
# if acceptable, write a 'munged' form to an output file
# if not, write the line to a separate file
# If desired, take the mean of numerical values
# Output a file with a header, with sequences in 'Seq', all others outputs to be predicted,
# with names specified in configuration.
import sys
import os
import pandas
import numpy as np
import ConfigParser

def munge_one_line(l, mins, maxes, m_range, take_mean):
  l = l.strip().split(',')
  seq = l[0]
  means = np.array([float(q) for q in l[1:]])
  for i, pair in enumerate(zip(mins, maxes)):
    this_min, this_max = pair
    if means[i] < this_min or means[i] > this_max:
      return(None)
  this_range = np.max(means) - np.min(means)
  if np.abs(this_range) > m_range:
    return(None)
  # all tests are passed at this point
  if take_mean:
    final_mean = np.mean(means)
    return(seq + ',' + str(final_mean) + '\n')
  else:
    means = [str(q) for q in means]
    return(seq + ',' + ','.join(means) + '\n')

def main(fn_in, fn_out, fn_fail, mins, maxes, m_range, take_mean, outnames):
  with open(fn_in, 'r') as fi, open(fn_out, 'w') as fo, open(fn_fail, 'w') as fx:
    fo.write('Seq,' + ','.join(outnames) + '\n')
    for l in fi:
        munged_l = munge_one_line(l, mins, maxes, m_range, take_mean)
        if munged_l == None:
          fx.write(l)
        else:
          fo.write(munged_l)

if __name__ == '__main__':
  config = ConfigParser.RawConfigParser()
  config.read(sys.argv[1])
  fn_in = os.path.expanduser(config.get('Files','file_in'))
  fn_out = os.path.expanduser(config.get('Files','file_out'))
  fn_fail = os.path.expanduser(config.get('Files','file_fail'))
  take_mean = config.get('Params','mean_reps') == 'True'
  m_range = float(config.get('Params','replicate_diff'))
  mins = [float(q) for q in config.get('Params','min').strip().split(',')]
  maxes = [float(q) for q in config.get('Params','max').strip().split(',')]
  outnames = config.get('Params','outnames').strip().split(',')
  assert(len(mins) == len(maxes))
  assert(len(mins) == len(outnames))
  if take_mean:
    assert(len(outnames) == 1)

  main(fn_in, fn_out, fn_fail, mins, maxes, m_range, take_mean, outnames)
