# drop bins not in the sets of 12; for now, just take average of bins 8
# Apply thresholding by total read count after munging
import sys
import numpy as np
import pandas
from time import time

def threshold_df(df_in, thresh):
  df = df_in.copy()
  # drop all non-integer columns
  for i in df.columns.values:
    if not isint(i):
      df = df.drop(i, axis = 1)

  colnames = np.array([int(q) for q in df.columns.values])
  read_arr = np.array(df).astype('float')
  read_arr = read_arr[:,np.argsort(colnames)]
  read_cts_a = np.apply_along_axis(np.sum,1,read_arr[:,:12])
  read_cts_b = np.apply_along_axis(np.sum,1,read_arr[:,12:])
  rows_useful = np.logical_and(read_cts_a >= thresh, read_cts_b >= thresh)
  return(df_in.loc[rows_useful])

def isint(q):
  try:
    int(q)
    return True
  except ValueError:
    return False

def reduce_by_keys(dat, key_pair, frac_first = 0.5):
  assert(len(key_pair) == 2)
  dat_0 = [q*frac_first for q in dat[key_pair[0]].tolist()]
  dat_1 = [q*(1.-frac_first) for q in dat[key_pair[1]].tolist()]
  dat_to_merge = np.stack([dat_0, dat_1], axis = 0)
  dat_merged = np.sum(dat_to_merge, axis = 0)
  dat[key_pair[0]] = dat_merged
  dat = dat.drop( columns = [key_pair[1]] )
  return(dat)

def renumber(dat):
  k = list(dat.columns.values)
  k = [int(q) for q in k if isint(q)]
  k.sort()
  new_vals = [q for q in range(len(k))]
  for p, q in zip(k, new_vals):
    dat = dat.rename(index=str, columns={str(p):str(q)})
  return(dat)

def drop_excess(dat):
  k = list(dat.columns.values)
  k = [q for q in k if isint(q)]
  k.sort()
  for q in k:
    if int(q) >= 24:
      dat = dat.drop(q, axis = 1)
  return(dat)

def main(fn_in, fn_out, thresh):
  t = time()
  dat = pandas.read_csv(fn_in)
  print('Read csv: ' + str(time() - t))
  t = time()
  for q in list(dat.columns.values):
    if not isint(q) and q != 'Seq':
      dat = dat.drop(columns = [q])
  print('Drop non-int columns: ' + str(time() - t))
  t = time()
  unin_8 = ['7','8']
  in_8 = ['20','21']
  dat = reduce_by_keys(dat, unin_8)
  print('Reduce A: ' + str(time() - t))
  t = time()
  dat = reduce_by_keys(dat, in_8)
  print('Reduce B: ' + str(time() - t))
  t = time()
  dat = renumber(dat)
  print('Renumber: ' + str(time() - t))
  t = time()
  dat = drop_excess(dat)
  print('Drop excess: ' + str(time() - t))
  t = time()
  dat = threshold_df(dat, thresh)
  print('Threshold: ' + str(time() - t))
  print(dat)
  dat.to_csv(fn_out)

if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
