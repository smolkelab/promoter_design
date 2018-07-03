import sys
import pandas

# convert a padded de-one-hotted sequence to the original seq
def extract_seq(seq):
  seq_cut = seq[25:-24]
  try:
    assert(len(seq_cut) == 313)
  except:
    print('Seq: ' + seq)
    print(len(seq))
    raise Exception
  return(seq_cut)

def main(*args):
  [data_fn, test_seqs_fn, test_vals_fn, out_fn] = args[0]
  dat = pandas.read_csv(data_fn) # keys: 'Means_A','Means_B','Seqs'
  seqs = pandas.DataFrame({'Seqs': [extract_seq(q.strip().split(',')[1]) for q in open(test_seqs_fn, 'r') if len(q) > 10]})
  preds = pandas.read_csv(test_vals_fn) # keys: 'Means', 'Preds'
  preds = pandas.concat([preds, seqs], axis = 1) # keys: 'Means', 'Preds', 'Seqs'

  dat = pandas.merge(dat, preds, how = 'inner', on = 'Seqs')
  dat.to_csv(out_fn)

if __name__ == '__main__':
  main(sys.argv[1:])

# python join_preds_seqs.py FS7_nextseq_means.csv test_seqs.txt FS7_preds.csv FS7_preds_seqs.csv
