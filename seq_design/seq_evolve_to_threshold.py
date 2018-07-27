# Evolve a set of sequences; when one reaches a target threshold score, remove it (replace it with a new random sequence)
# and add the good one to an output.
# Evolve for a given number of rounds, or until enough target sequences are found.
import sys
import os
import pandas
import numpy as np
import random
from numpy.random import choice
import ConfigParser
import types
import seq_evolution
import seq_selection

class seq_evolution_thresh(seq_evolution.seq_evolution_class):
  def __init__(self, cfg):
    super(seq_evolution_class_gradient, self).__init__(cfg)
    
  # Iterate, but return and reset any sequences scoring higher than 'thresh'.
  # Reset the iteration counter for that sequence as well.
  def thresholding_iterative(self, thresh, seqs_final, num_iters):
    seqs_out = []
    while(len(seq_out) < seqs_final):
      # Round the seqs to one-hot, test them, merge the outputs from each model, and merge those to one final output per sequence.
      model_scores = self.params['merge_models'](self.params['merge_outputs'](self._test_sequences(self.round_seqs(self.seqs))))
      seq_scores = np.apply_along_axis(self.params['seq_scores'], 0, self.seqs)
      scores = model_scores + seq_scores
      done_pos = scores > thresh
      seq_out.extend([q for q in thresholded_evolver.seqs[done_pos]])
      for j, q in enumerate(done_pos):
        if q or self.curr_iters[j] == num_iters - 1: # reset this sequence if needed
          self.seqs[j] = self._populate_one_sequence()
          self.curr_iters[j] = 0
        else:
          self.curr_iters[j] += 1
      self.iterate()
    seqs_out = seqs_out[:seqs_final]
    self.seqs = np.array(seqs_out)

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  num_iters = int(cfg.get('Params','NUM_ITERS'))
  thresh = float(cfg.get('Params','THRESH'))
  seqs_final = int(cfg.get('Params','NUM_SEQS_FINAL'))
  thresholded_evolver = seq_evolution_thresh(cfg)
  thresholded_evolver.choose_best_seqs = types.MethodType(eval(cfg.get('Functions','choose_best_seqs')), thresholded_evolver)
  thresholded_evolver.thresholding_iterative(thresh, seqs_final, num_iters)
  thresholded_evolver.save_output()
  seq_selection.main(cfg)
