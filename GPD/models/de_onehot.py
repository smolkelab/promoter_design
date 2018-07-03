import sys
import numpy as np
import pandas

def main(fi, fo):
  key = ['A','C','G','T']
  x = np.load(fi)
  x = np.argmax(x, axis = 1)
  seqs = []
  for i in x:
    seq = ''.join([key[q] for q in i])
    seqs.append(seq)
  pandas.DataFrame({'Seq':seqs}).to_csv(fo)

if __name__ == '__main__':
  fi, fo = sys.argv[1:]
  main(fi, fo)
