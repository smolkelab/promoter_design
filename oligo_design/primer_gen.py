# Generate primers of a specified length and annealing temperature,
# with a specified distance from other primers.
# Accept a random seed and starting position, to guarantee the results can be repeated.
#import primer3
import ConfigParser
import sys
import os
import math
import numpy as np
import pandas as pd
import copy
import primer3

DNA = ['A','C','G','T']
DNA_C = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
N_PROBS = [0.25, 0.25, 0.25, 0.25]
MAX_NN_LENGTH = 60 # primer3 default
PRIMER3_METHOD = 'santalucia' # primer3 default
MAX_LOOP = 30 # primer3 default

def rc(seq):
  return(''.join([DNA_C[q] for q in seq])[::-1])

class primer_set(object):

  def __init__(self, cfg):
    params = dict(cfg.items('Params'))
    self.toe_anneal = float(params['toe_anneal'])
    self.toe_tolerance = float(params['toe_tolerance'])
    self.full_anneal = float(params['full_anneal'])
    self.full_tolerance = float(params['full_tolerance'])
    self.misprime_tolerance = float(params['misprime_tolerance'])
    self.evolve_iterations = int(params['evolve_iterations']) # number of iterations allowed for one call to 'evolve_new_seq' to get one new seq
    self.seqs = {'Toehold':[], 'Full':[]}
    self.seq_base = params['SEQ']
    self.base_probs = {'N': N_PROBS} # base probabilities for the toehold: include e.g. BsaI sites this way
    for q in self.seq_base:
      if not q in DNA and not q in self.base_probs.keys():
        x = [float(q) for q in params[q].strip().split(':')]
        self.base_probs[q] = x/sum(x)
    self.len_full_add = int(params['len_full_add'])# length of sequence to add to toehold to get something that anneals at self.full_anneal
    self.add_to_left = params['add_to_left'] == 'True'
    self.mv_conc = float(params['mv_conc']) # 50
    self.dv_conc = float(params['dv_conc']) # 0
    self.dntp_conc = float(params['dntp_conc']) # 0.8
    self.dna_conc = float(params['dna_conc']) # 50
    np.random.seed(int(params['random_seed']))

  # Get the annealing temperature of two sequences
  def melting_temp(self, seq):
    return( primer3.calcTm(seq, mv_conc=self.mv_conc, dv_conc=self.dv_conc, dntp_conc=self.dntp_conc,
                      dna_conc=self.dna_conc, max_nn_length=MAX_NN_LENGTH, tm_method=PRIMER3_METHOD, salt_corrections_method=PRIMER3_METHOD) )

  # Find whether a concerning heterodimer exists between two sequences
  def get_misprime_tm(self, seq1, seq2, temp):
    tm_f1 = self.get_misprime_tm_inner(seq1, seq2, temp)
    tm_r1 = self.get_misprime_tm_inner(seq1, rc(seq2), temp)
    tm_f2 = self.get_misprime_tm_inner(seq2, seq1, temp)
    tm_r2 = self.get_misprime_tm_inner(seq2, rc(seq1), temp)
    return(max(tm_f1, tm_r1, tm_f2, tm_r2))

  def get_misprime_tm_inner(self, seq1, seq2, temp):
    heterodimer = primer3.calcHeterodimer(seq1, seq2,
                    mv_conc=self.mv_conc, dv_conc=self.dv_conc, dntp_conc=self.dntp_conc, dna_conc=self.dna_conc, 
                    temp_c=temp, max_loop=MAX_LOOP)
    if not heterodimer.structure_found:
      return -100.
    return heterodimer.tm

  def true_if_no_misprime(self, seq1, seq2, temp):
    struct_tm = self.get_misprime_tm(seq1, seq2, temp)
    #if(struct_tm >= (temp - self.misprime_tolerance)):
    #  raise Exception('oops!' + str(seq1) + ' ' + str(seq2) + ' ' + str(struct_tm))

    return struct_tm < (temp - self.misprime_tolerance) # e.g. if temp = 45, self.misprime_tolerance = 5, OK if t_anneal < 40; will be melted by 45

  # placeholder, for maybe rejecting e.g. hairpins
  def primer_passed_filter(self, seq):
    return True

  # generate a random sequence
  def get_random_seq(self, parent = None):
    if parent == None: # this is Python's fault
      parent = self.seq_base
    ans = []
    for b in parent:
      if b in DNA:
        ans.append(b)
      else:
        ans.append( np.random.choice(DNA, None, self.base_probs[b]) )
    return(''.join(ans))

  # get all variants of a random sequence, and find the one with Tm closest to some target
  # Compare to a template to see which bases aren't allowed to be changed
  def evolve_one_step(self, seq, template, target):
    assert(len(seq) == len(template))
    ans = seq
    err = abs( self.melting_temp(seq) - target )
    for i in range(len(seq)):
      if template[i] not in DNA:
        for b in DNA:
          if seq[i] != b:
            tmp = [q for q in seq]
            tmp[i] = b
            tmp = ''.join(tmp)
            new_err = abs(self.melting_temp(tmp) - target)
            if new_err < err and self.primer_passed_filter(tmp):
              err = new_err
              ans = tmp
    return(ans)

  def true_if_no_misprimes_to_seqs(self, seq):
    if len(self.seqs['Full']) == 0:
      return True
    return(all([self.true_if_no_misprime(seq, q, self.toe_anneal) for q in self.seqs['Full']]))

  # evolve a valid new toe first, then extend
  # 'self.evolve_iterations' available to get both tasks done.
  def evolve_new_seq_inner(self):
    seq = self.get_random_seq()
    is_full = False
    template = self.seq_base
    toe = None; full = None
    for i in range(self.evolve_iterations):
      if is_full: # get temperature targets based on what part of the process this is
        target = self.full_anneal
        target_low = target - self.full_tolerance; target_high = target + self.full_tolerance
      else:
        target = self.toe_anneal
        target_low = target - self.toe_tolerance; target_high = target + self.toe_tolerance
      tm = self.melting_temp(seq)
      if tm > target_low and tm < target_high:
        if not is_full:
          is_full = True; toe = seq
          new_pad = ''.join(['N']*self.len_full_add)
          if self.add_to_left:
            template = new_pad + seq
            seq = self.get_random_seq(template)
          else:
            template = seq + new_pad
            seq = self.get_random_seq(template)
        else:
          full = seq
          return(toe, full)
      seq = self.evolve_one_step(seq, template, target)
    return(None, None)

  def evolve_new_seq(self):
    toe, full = self.evolve_new_seq_inner()
    if full != None:
      if self.true_if_no_misprimes_to_seqs(full):
        self.seqs['Toehold'].append(toe)
        self.seqs['Full'].append(full)

def get_seqs(cfg):
  num_seqs = int(cfg.get('Params','num_seqs'))
  num_iters = int(cfg.get('Params','num_iters'))
  primers = primer_set(cfg)
  for i in range(num_iters):
    primers.evolve_new_seq()
    if len(primers.seqs['Full']) >= num_seqs:
      return(primers)
  if len(primers.seqs['Full']) < num_seqs:
    raise Exception(str(num_seqs) + ' primers desired, but ' + str(len(primers.seqs)) + ' found.')
  return(primers)

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True); cfg.optionxform=str
  cfg.read(sys.argv[1])
  primers = get_seqs(cfg)
  for (toe, full) in zip(primers.seqs['Toehold'], primers.seqs['Full']):
    print('Full: ' + str(full))
    print('Toe: ' + str(toe))
    print('T_melt (toe): ' + str(primers.melting_temp(toe)))
    print('T_melt (full): ' + str(primers.melting_temp(full)))
    t_ann = float(cfg.get('Params','toe_anneal'))
    worst_seq = ''
    worst_t = -1000.
    for q in primers.seqs['Full']:
      if q != full:
        new_t = primers.get_misprime_tm(full, q, t_ann)
        if new_t > worst_t:
          worst_t = new_t; worst_seq = q

    print('T_misprime (max): ' + str(worst_t))
    print('Worst misprime sequence: ' + str(worst_seq))
  print('Sequences generated: ' + str(len(primers.seqs['Full'])))
  if len(sys.argv) > 2:
    df_out = pd.DataFrame(primers.seqs)
    df_out.to_csv(sys.argv[2], index = False)
