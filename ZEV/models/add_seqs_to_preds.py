import sys
import numpy as np
import pandas

def main(f_seqs, f_preds, fo):
  key = ['A','C','G','T']
  x = np.load(f_seqs)
  x = np.argmax(x, axis = 1)
  seqs = []
  for i in x:
    seq = ''.join([key[q] for q in i])
    seqs.append(seq)
  df = pandas.read_csv(f_preds)
  df['Seq'] = seqs
  df.to_csv(fo)

if __name__ == '__main__':
  fs, fp, fo = sys.argv[1:]
  main(fs, fp, fo)
