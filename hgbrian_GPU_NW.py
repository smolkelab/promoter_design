# Approach: for each sequence, call a function modifying hgbrian's NW code
# to get the score.

import numba
from numba import cuda, float32
import numpy as np

# for testing
import nw_align_mod
import time

import sys
import copy
from multiprocessing import Process, Pipe

# Running on a 960M GPU: CC 5.0, 32 blocks, 2048 threads per SM; not using shared memory so don't have to optimize that. 64 threads per block should give full occupany; 32 T.P.B. would be 50%.
THREADS_PER_BLOCK = 64 
A_SMALL_FLOAT = -100000.

# Constant matrices for DNA-to-integer conversion
NTLIST = ["","A","C","G","T","N","-"]
NTDICT = dict((j,i) for i,j in enumerate(NTLIST))

N_VAL = ord('N')

# For reference: it looks like Python, but it ain't.
# https://docs.anaconda.com/numbapro/CUDAPySpec
# Written as a formattable string, not just code, because cuda.local.array requires the size to be specified with raw constants.
# And why not pass in the other parameters that way while I'm at it?
KERNEL_STR =  """@cuda.jit
def hgbrian_score_kernel(seq_arr, len_arr, out_arr):
  max_gaps = {max_gaps}
  gap= {gap} #-1
  match= {match} #1
  mismatch= {mismatch} #0
  eps = 1e-8

  pos = cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x
  if pos < seq_arr.size - 1:  # Check array boundaries
      max_i = len_arr[pos]
      max_j = len_arr[pos+1]
      # Break immediately if alignment impossible
      if abs(max_i - max_j) > max_gaps:
        out_arr[pos] = A_SMALL_FLOAT
      else:
          # ensure seq_i is the longer of the two
          if max_j > max_i:
            t = max_i
            max_i = max_j
            max_j = t 
            seqi = seq_arr[pos+1,...][:max_i]
            seqj = seq_arr[pos,...][:max_j]   
          else:
            seqi = seq_arr[pos,...][:max_i]
            seqj = seq_arr[pos+1,...][:max_j]


          last_score = cuda.local.array({local_size}, dtype=float32) # 2*max_gaps + 1 # could be made faster using shared memory, but I don't think the kernel is the bottleneck
          score = cuda.local.array({local_size}, dtype=float32) # 2*max_gaps + 1

          for q in range(len(last_score)):
            score[q] = abs(q-max_gaps)*gap

          for i in range(0, max_i):
            ci = seqi[i]
            for p in range(len(last_score)):
                last_score[p] = score[p]
            for q in range(len(score)):
              # Ignore the spaces in the grid that require more than max_gaps indels
              j = i - max_gaps + q
              if j < 0 or j >= max_j:
                score[q] = A_SMALL_FLOAT
              else:
                cj = seqj[j]
                _matchscore = match if ci == cj and ci != {N} else mismatch # {N} is the 'N' character, which matches nothing

                dg_score = last_score[q] + _matchscore
                if q != len(score) - 1:
                  up_score = last_score[q+1] + gap
                else:
                  up_score = A_SMALL_FLOAT
                if q != 0:
                  lf_score = score[q-1] + gap
                else:
                  lf_score = A_SMALL_FLOAT

                if dg_score >= up_score-eps and dg_score >= lf_score-eps:
                  score[q] = dg_score
                elif lf_score >= dg_score-eps and lf_score >= up_score-eps:
                  score[q] = lf_score
                else:
                  score[q] = up_score
          # Find which position in 'score' corresponds to the final alignment ([max_i, max_j] in the original explicit matrix)
          ans = score[max_j - max_i + max_gaps]
          out_arr[pos] = ans"""

def _convert_dna_to_array(dnastr, max_len):
  ans = np.zeros(max_len)
  for i in range(len(dnastr)):
    ans[i] = NTDICT[dnastr[i]]
  return(ans)

def hgbrian_nw_score_wrapper(seqs, max_gaps):
  local_size = 2*max_gaps + 1
  ks_formatted = KERNEL_STR.format(max_gaps = str(max_gaps), local_size = str(local_size), gap = -5, match = 1, mismatch = 0, N = str(N_VAL))
  exec(ks_formatted)
  len_arr = np.array([len(q) for q in seqs])
  max_len = np.max(len_arr)
  scores = np.zeros(shape = len(seqs) - 1, dtype = 'float')
  cuda.to_device(seqs_arr)
  cuda.to_device(len_arr)
  cuda.to_device(scores)
  blockspergrid = (len(seqs) + (THREADS_PER_BLOCK - 1)) // THREADS_PER_BLOCK
  hgbrian_score_kernel[blockspergrid, THREADS_PER_BLOCK](seqs_arr, len_arr, scores)
  return(scores)

def seq_to_list(seq, max_len):
  seq = [ord(q) for q in seq]
  seq.extend([0]*(max_len - len(seq)))
  return(seq)

def seqs_to_arrays(seqs, max_len, conn):
  seq_arr = [seq_to_list(q,max_len) for q in seqs]
  seq_arr = np.array(seq_arr)
  conn.send(seq_arr)
  conn.close()

# Ugly, but seems to work. Manually spawn three processes, and split the work up among them.
def seqs_to_arrays_triplicate_wrapper(chunk, max_len):
  split_pos = len(chunk)/3
  lines_1 = copy.copy(chunk[:split_pos])
  lines_2 = copy.copy(chunk[split_pos:2*split_pos])
  lines_3 = copy.copy(chunk[2*split_pos:])

  parent_1, child_1 = Pipe()
  p1 = Process(target = seqs_to_arrays, args = (lines_1, max_len, child_1))
  p1.start()
  parent_2, child_2 = Pipe()
  p2 = Process(target = seqs_to_arrays, args = (lines_2, max_len, child_2))
  p2.start()
  parent_3, child_3 = Pipe()
  p3 = Process(target = seqs_to_arrays, args = (lines_3, max_len, child_3))
  p3.start()

  array_1 = parent_1.recv()
  array_2 = parent_2.recv()
  array_3 = parent_3.recv()
  p1.join()
  p2.join()
  p3.join()

  arrs_out = []
  for q in [array_1, array_2, array_3]:
    if q.shape[0] > 0:
      arrs_out.append(q)
  arrs_out = np.concatenate(arrs_out, axis = 0)

  return(arrs_out)


# Read lines from file_in; write scores to file_out
# Approach: iteratively read lines from file_in; pass them as conveniently-sized blocks
# to the device and invoke the kernel to align them.
# Written for clarity, not perfect I/O efficiency - CPU-bound preprocessing is bottleneck for me.
def score_file(file_in, file_out, max_gaps, chunk_size, gap = -5, match = 1, mismatch = 0):
  local_size = 2*max_gaps + 1
  ks_formatted = KERNEL_STR.format(max_gaps = str(max_gaps), local_size = str(local_size), gap = gap, match = match, mismatch = mismatch, N = str(N_VAL))
  exec(ks_formatted)
  chunk = []
  lens = []
  max_len = 0
  stream = cuda.stream() # not really necessary in this setup, but keeping in case I want to adapt this code
  with open(file_in) as fi, open(file_out, 'w') as fo:
    keep_reading = True
    while(keep_reading):
      # fill the next chunk
      while(len(chunk) < chunk_size):
        try:
          l = fi.next().split(',')[0]
          chunk.append(l)
          lens.append(len(l))
          if max_len < len(l):
            max_len = len(l)
        except StopIteration:
          keep_reading = False
          break
      if len(chunk) <= 1: # need at least 2 seqs for a comparison
        break

      # convert the chunk to arrays
      scores_local = np.zeros(shape = len(chunk) - 1, dtype = 'float')
      seqs_arr_local = seqs_to_arrays_triplicate_wrapper(chunk, max_len)
      len_arr_local = np.array(lens)
    
      # send the arrays to the device
      stream.synchronize()
      seqs_arr_shared = np.copy(seqs_arr_local)
      len_arr_shared = np.copy(len_arr_local)
      scores_shared = np.copy(scores_local)
      seqs_arr_shared, len_arr_shared, scores_shared = cuda.to_device(seqs_arr_shared, stream = stream), cuda.to_device(len_arr_shared, stream = stream), cuda.to_device(scores_shared, stream = stream)
      blockspergrid = (seqs_arr_shared.shape[0] + (THREADS_PER_BLOCK - 1)) // THREADS_PER_BLOCK
      hgbrian_score_kernel[blockspergrid, THREADS_PER_BLOCK, stream](seqs_arr_shared, len_arr_shared, scores_shared)
      stream.synchronize() # <- if the kernel was more time-intensive, we could revise this to do preprocessing while the kernel runs - makes bookkeeping to get the output back nontrivial and bug-prone though.

      # Copy the results to host, and write them
      scores_local = np.copy(scores_shared.copy_to_host())
      for q in range(scores_local.shape[0]):
        fo.write("%s\n" % str(scores_local[q]))
      write_scores = False

      # prepare for the next chunk; start the next chunk with the last sequence in this one, as we don't know the distance to its successor yet.
      chunk = [chunk[-1]]
      lens = [lens[-1]]
      max_len = len(chunk[0])


def test_nw_align_to_score(aln):
  score = 0
  for (i,j) in zip(aln[0], aln[1]):
    if i == '-' or j == '-':
        score = score - 1
    elif i == j and i != 'N':
        score = score + 1
  return(score)

if __name__ == '__main__':
    file_in, file_out, max_gaps, chunk_size = sys.argv[1:]
    max_gaps, chunk_size = int(max_gaps), int(chunk_size)
    score_file(file_in, file_out, max_gaps, chunk_size)