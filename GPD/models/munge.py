
# For each line in an input file,
# if acceptable, write a 'munged' form to an output file
# if not, write the line to a separate file
import sys
import pandas
import numpy as np

MUNGE_THRESH = 0.2

# Count the lines in a file:
# https://stackoverflow.com/questions/19001402/how-to-count-the-total-number-of-lines-in-a-text-file-using-python#19001475
def count_lines(file_in):
  with open(file_in) as f:
    ans = sum(1 for _ in f)
  return(ans)

def progress(count, total, status=''):
  bar_len = 60
  filled_len = int(round(bar_len * count / float(total)))
  percents = round(100.0 * count / float(total), 1)
  bar = '=' * filled_len + '-' * (bar_len - filled_len)
  sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
  sys.stdout.flush() # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)

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
  num_lines = count_lines(fn_in)
  curr = 0
  seq_id = None
  unmungeables = 0
  with open(fn_in, 'r') as fi, open(fn_out, 'w') as fo, open(fn_bad, 'w') as fx:
    fo.write('Seq,Mean\n')
    for l in fi:
      curr += 1
      if curr % 1e5 == 0:
        progress(curr, num_lines)
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
  print('Unmungeables: ' + str(unmungeables))

if __name__ == '__main__':
  fn_in, fn_out, fn_bad = sys.argv[1:]
  munge(fn_in, fn_out, fn_bad)

# time python munge.py FS7_nextseq_means.csv FS7_means_trainable.csv FS7_nextseq_means_unmungeable.csv
