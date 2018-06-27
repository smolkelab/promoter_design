# Importing in the originall specified order (numba, numpy, torch, cupy) causes dlopen-related errors on my machine.
# Reordering the imports may help.

#import numba
#import numpy as np
#import torch as tr
#import cupy as cp

#import torch as tr
import numpy as np
import numba
#import cupy as cp


import sys

# enum up, left, diagonal, none
UP, LF, DG, NO = range(4)

# from string to uint8 list/array and back again
ntlist = ["","A","C","G","T","N","-"]

ntdict = dict((j,i) for i,j in enumerate(ntlist))

eps = 1e-8

def _convert_dna(dnastr):
    return [ntdict[nt] for nt in dnastr]

def _unconvert_dna(dnai):
    return ''.join(ntlist[nti] for nti in dnai)


def global_align(seqi, seqj, gap=-1, match=1, mismatch=-1, mode="torch"):
    """Needleman-Wunsch
    """

    if mode == "numba":
        return global_align_numba(seqi, seqj, gap=gap, match=match, mismatch=mismatch, mode='numpy')

    seqja = _convert_dna(seqj)
    seqia = _convert_dna(seqi)

    #if mode in ("torch", "torchcuda"):
    #    seqj = tr.from_numpy(np.array(seqja, dtype=np.uint8))
    #    seqi = tr.from_numpy(np.array(seqia, dtype=np.uint8))
    #    if mode == "torchcuda":
    #        seqi = seqi.cuda()
    #        seqj = seqj.cuda()
    #    max_j = seqj.size()[0]
    #    max_i = seqi.size()[0]
    #elif mode == "cupy":
    #    seqj = cp.array(seqja, dtype=np.uint8)
    #    seqi = cp.array(seqia, dtype=np.uint8)
    #    max_j = seqj.size
    #    max_i = seqi.size
    #elif mode == "numpy":
    if mode == 'numpy':
        seqj = np.array(seqja, dtype=np.uint8)
        seqi = np.array(seqia, dtype=np.uint8)
        max_j = seqj.shape[0]
        max_i = seqi.shape[0]
    else:
        raise ValueError

    score = np.zeros((max_i + 1, max_j + 1), dtype=np.float)
    point = np.zeros((max_i + 1, max_j + 1), dtype=np.uint8)

    score[0, 0] = 0.0
    point[0, 0] = NO
    point[0, 1:] = LF
    point[1:, 0] = UP

    score[0, :] = gap * np.arange(max_j + 1)
    score[:, 0] = gap * np.arange(max_i + 1)

    for i in range(1, max_i+1):
        ci = seqi[i-1]
        for j in range(1, max_j+1):
            cj = seqj[j-1]

            _matchscore = match if ci == cj and ci != 5 else mismatch

            dg_score = score[i-1, j-1] + _matchscore
            up_score = score[i-1, j] + gap
            lf_score = score[i, j-1] + gap

            if dg_score >= up_score-eps and dg_score >= lf_score-eps:
                score[i, j] = dg_score
                point[i, j] = DG
            elif lf_score >= dg_score-eps and lf_score >= up_score-eps:
                score[i, j] = lf_score
                point[i, j] = LF
            else:
                score[i, j] = up_score
                point[i, j] = UP

    i = max_i
    j = max_j
    #if mode in ("torch", "torchcuda"):
    #    align_i = tr.from_numpy(np.zeros(max_i + max_j, dtype=np.uint8))
    #    align_j = tr.from_numpy(np.zeros(max_i + max_j, dtype=np.uint8))
    #    if mode == "torchcuda":
    #        align_i = align_i.cuda()
    #        align_j = align_j.cuda()
    #elif mode == "cupy":
    #    align_i = cp.zeros(max_i + max_j, dtype=np.uint8)
    #    align_j = cp.zeros(max_i + max_j, dtype=np.uint8)
    #elif mode == "numpy":
    if mode == 'numpy':
        align_i = np.zeros(max_i + max_j, dtype=np.uint8)
        align_j = np.zeros(max_i + max_j, dtype=np.uint8)
    else:
        raise ValueError

    nij = 0
    while i > 0 or j > 0:
        p = point[i, j]
        s = score[i, j]

        if p == DG:
            align_j[nij] = seqj[j-1]
            align_i[nij] = seqi[i-1]
            i -= 1
            j -= 1
        elif p == LF:
            align_j[nij] = seqj[j-1]
            align_i[nij] = ntdict['-']
            j -= 1
        elif p == UP:
            align_j[nij] = ntdict['-']
            align_i[nij] = seqi[i-1]
            i -= 1
        else:
            # numba issues here with raising
            print("ERR", p) #raise ValueError('unexpected p: {}'.format(p))
            break
        nij += 1

    return _unconvert_dna(align_i.tolist()[::-1]), _unconvert_dna(align_j.tolist()[::-1])

# just jit it
global_align_numba = numba.jit(global_align)


def test_random(n=1000):
    import time
    from random import choice
    def _rseq(n): return [choice("ACGT") for _ in range(n)]

    rseq1, rseq2 = _rseq(n), _rseq(n)

    assert tr.cuda.is_available(), "cuda needed for torch and cupy"

    for mode in ["torch", "torchcuda", "numpy", "numba", "cupy"]:
        print("mode", mode)
        if mode == "cupy" and n > 500: print("{} is too slow".format(mode)); continue

        start = time.time()
        al1, al2 = global_align(rseq1, rseq2, mode=mode)
        print(al1[:100])
        print(al2[:100])
        print("{:.3g}".format((time.time()-start)))


if __name__ == "__main__":
    if sys.argv[1] in ("numpy", "numba", "torch", "torchcuda", "cupy"):
        mode = sys.argv[1]
        al1, al2 = global_align(sys.argv[2], sys.argv[3], mode=mode)
        print(al1)
        print(al2)
    elif sys.argv[1][0] in "ACGT":
        mode = "numpy"
        al1, al2 = global_align(sys.argv[1], sys.argv[2], mode=mode)
        print(al1)
        print(al2)
    else:
        print("random", sys.argv[1])
        test_random(int(sys.argv[1]))

