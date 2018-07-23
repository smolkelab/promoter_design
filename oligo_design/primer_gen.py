# Generate primers of a specified length and annealing temperature,
# with a specified distance from other primers.
# Accept a random seed and starting position, to guarantee the results can be repeated.
#import primer3
import ConfigParser
import sys
import os
import math
import numpy as np
import copy

DNA = ['A','C','G','T']
N_PROBS = [0.25, 0.25, 0.25, 0.25]
MAX_NN_LENGTH = 60 # primer3 default
PRIMER3_METHOD = 'santalucia' # primer3 default
MAX_LOOP = 30 # primer3 default

class primer_set(object):

  def __init__(self, cfg)
    params = dict(cfg.items('Params'))
    self.toe_anneal = float(params['toe_anneal'])
    self.toe_tolerance = float(params['toe_tolerance'])
    self.full_anneal = float(params['full_anneal'])
    self.full_tolerance = float(params['full_tolerance'])
    self.misprime_tolerance = float(params['misprime_tolerance'])
    self.evolve_iterations = int(params['evolve_iterations']) # number of iterations allowed for one call to 'evolve_new_seq' to get one new seq
    self.seqs = []
    self.seq_base = params['SEQ']
    self.base_probs = {'N': N_PROBS} # base probabilities for the toehold: include e.g. BsaI sites this way
    for q in self.seq_base:
      if not q in DNA:
        x = [float(q) for q in params[q].strip().split(':')]
        self.base_probs[q] = x/sum(x)
    self.len_full_add = int(params['len_full_add'])# length of sequence to add to toehold to get something that anneals at self.full_anneal
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
    heterodimer = primer3.calcHeterodimer(seq1, seq2,
                    mv_conc=self.mv_conc, dv_conc=self.dv_conc, dntp_conc=self.dntp_conc, dna_conc=self.dna_conc, 
                    temp_c=temp, max_loop=MAX_LOOP])
    if not heterodimer.structure_found:
	  return 0.
	return heterodimer.tm
  
  def true_if_no_misprime(self, seq1, seq2, temp):
    struct_tm = self.get_misprime_tm(seq1, seq2, temp)
    return struct_tm < temp - self.misprime_tolerance # e.g. if temp = 45, self.misprime_tolerance = 5, OK if t_anneal < 40; will be melted by 45
    
  # placeholder, for maybe rejecting e.g. hairpins
  def primer_passed_filter(self, seq):
    return True
    
  # generate a random sequence
  def get_random_seq(self, parent = self.base_probs):
    ans = []
    for b in parent:
      if b in DNA:
        ans.append(b)
      else:
        ans.append( np.random.choice(DNA, None, self.base_probs[b]) )
    return(''.join(ans))
    
  # get all variants of a random sequence, and find the one with Tm closest to some target
  def evolve_one_step(self, seq, target):
    ans = seq
    err = abs( self.melting_temp(seq) - target )
    for i in len(seq):
      for b in DNA:
        if seq[i] != b:
          tmp = [q for q in seq]
          tmp[i] = b 
          new_err = self.melting_temp(tmp) - target
          if new_err < err and self.primer_passed_filter(tmp):
            err = new_err
            ans = tmp
    return(ans)
        
  def true_if_no_misprimes_to_seqs(self, seq):
    target_tmp = self.toe_anneal
    if len(self.seqs) == 0:
      return True
    return(all([self.true_if_no_misprime(seq, q, target_tmp) for q in self.seqs]))
    
  # evolve a valid new toe first, then extend
  # 'self.evolve_iterations' available to get both tasks done.
  def evolve_new_seq_inner(self):
    seq = self.get_random_seq()
    is_full = False
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
          is_full = True
          seq = self.get_random_seq(''.join(['N']*self.len_full_add) + seq)
        else:
          return(seq)
      seq = self.evolve_one_step(seq, target)
    return(None)
    
  def evolve_new_seq(self):
    seq = self.evolve_new_seq_inner()
    if seq != None:
      if self.true_if_no_misprimes_to_seqs(seq):
      self.seqs.append(seq)

def get_seqs(cfg):
  num_seqs = int(cfg.get('Params','num_seqs'))
  num_iters = int(cfg.get('Params','num_iters'))
  primers = primer_set(cfg)
  for i in range(num_iters):
    primers.evolve_new_seq()
	if len(primers.seqs) >= num_seqs:
	  break
  if len(primers.seqs) < num_seqs:
    return(None)
  return(primers)

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  primers = get_seqs(cfg)
  for seq in primers.seqs:
    print('Seq: ' + str(seq))
	toe = seq[-len(cfg.get('Params','SEQ'):]
	print('Toe: ' + str(toe))
	print('T_melt (toe): ' + str(primers.melting_temp(toe)))
	print('T_melt (full): ' + str(primers.melting_temp(seq)))
	t_ann = float(cfg.get('Params','toe_anneal'))
	misprimes = [(primers.get_misprime_tm(seq, q, t_ann), q) for q in primers.seqs if q != seq]
	print('T_misprime (max): ' + str(max([p for (p,q) in misprimes]))
	print('Worst misprime sequence: ' + str([q for (p,q) in misprimes if p == max([p for (p,q) in misprimes])]))