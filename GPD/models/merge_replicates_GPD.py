# For each line in an input file,
# if acceptable, write a 'munged' form to an output file
# if not, write the line to a separate file
import sys
import pandas
import numpy as np

MUNGE_THRESH = 0.2

def munge_one_line(l, seq_id, mean_ids): # expect a header entry 'Seq'; others are numerical means
  l = l.strip().split(',')
  seq = l[seq_id]
  means = np.array( [float(q) for q in [ l[p] for p in mean_ids ] ])
  mungeable = np.abs(means[1] - means[0]) < MUNGE_THRESH
  if mungeable:
    mean = str(float(np.mean(means)))
    munged = seq + ',' + mean + '\n'
    return(munged, True)
  else:
    return(None, False)

def munge(fn_in, fn_out, fn_bad):
  seq_id = None
  unmungeables = 0
  with open(fn_in, 'r') as fi, open(fn_out, 'w') as fo, open(fn_bad, 'w') as fx:
    fo.write('Seq,Mean\n')
    for l in fi:
      if seq_id is None:
        header = l.strip().split(',')
        seq_id = [i for i,q in enumerate(header) if q == 'Seqs']
        seq_id = int(seq_id[0])
        mean_ids = [i for i,q in enumerate(header) if 'Mean' in q]
      else:
        munged, mungeable = munge_one_line(l, seq_id, mean_ids)
        if mungeable:
          fo.write(munged)
        else:
          fx.write(l)
          unmungeables += 1
  print('Lines: ' + str(num_lines))
  print('Rejected: ' + str(unmungeables))

if __name__ == '__main__':
  fn_in, fn_out, fn_bad = sys.argv[1:]
  munge(fn_in, fn_out, fn_bad)