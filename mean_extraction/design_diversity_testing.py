import sys
import os
import pandas as pd
import hgbrian_GPU_NW
import numba

FN_TMP_SEQS = os.path.expanduser('~/facs-seq_test/seqs_tmp.txt')
FN_TMP_SEQS2 = os.path.expanduser('~/facs-seq_test/seqs_tmp2.txt')
FN_TMP_SCORES = os.path.expanduser('~/facs-seq_test/scores_tmp.txt')
MAX_GAPS = 5
CHUNK_SIZE = 500000
THREADS_PER_BLOCK = 64
RANDOM_SEED = 2017


# get all distances within a list of sequences ('reads') - 
# input is a list of strings, output is a list of (seq 1, seq 2, scores)
def get_all_dists_in_list_of_reads(reads, fn_tmp_seqs = FN_TMP_SEQS, fn_tmp_scores = FN_TMP_SCORES, max_gaps = MAX_GAPS, chunk_size = CHUNK_SIZE, threads_per_block = THREADS_PER_BLOCK):

  with open(fn_tmp_seqs, 'w') as fo:
    for i in range(len(reads)):
      for j in range(i+1, len(reads)):
        fo.write(reads[i] + '\n')
        fo.write(reads[j] + '\n')

  hgbrian_GPU_NW.score_file(fn_tmp_seqs, fn_tmp_scores, max_gaps, chunk_size, threads_per_block)
  # clean up a little
  context = numba.cuda.current_context()
  context.reset()

  ans = []
  with open(fn_tmp_scores, 'r') as fi:
    for i in range(len(reads)):
      for j in range(i+1, len(reads)):
        score = float(fi.next().strip())
        ans.append((reads[i], reads[j], score))
        # drop the next line; we don't need it
        try:
          _ = fi.next()
        except StopIteration:
          pass

  os.remove(fn_tmp_seqs)
  os.remove(fn_tmp_scores)
  return(ans)

def main(dir_in, fn_out):
  fns = os.listdir(dir_in)
  fns = [os.path.join(dir_in, q) for q in fns]
  ans = []
  for f in fns:
    print f
    p = pd.read_csv(f)
    reads = p['Seqs']
    score_tuples = get_all_dists_in_list_of_reads(reads)
    score_tuples = [q + (f,) for q in score_tuples]
    ans.extend(score_tuples)
  with open(fn_out, 'w') as fo:
    fo.write('Read_1,Read_2,Score,Filename\n')
    for line in ans:
      line = ','.join([str(q) for q in line])
      fo.write(line + '\n')

if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2])
