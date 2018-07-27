
import sys
import os
import pandas
import numpy as np
import random
from numpy.random import choice, rand
import ConfigParser
import seq_evolution
import seq_selection
import seq_gradient_evolution
import seq_evolve_to_threshold

DNA = seq_evolution.DNA

# like seq_evolution_class_gradient, but with the thresholding_iterative method from seq_evolve_to_threshold
class seq_evolution_thresh_gradient(seq_gradient_evolution.seq_evolution_class_gradient, seq_evolve_to_threshold.seq_evolution_thresh):
  def __init__(self, cfg):
    super(seq_evolution_thresh_gradient, self).__init__(cfg)

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  num_iters = int(cfg.get('Params','NUM_ITERS'))
  thresh = float(cfg.get('Params','THRESH'))
  seqs_final = int(cfg.get('Params','NUM_SEQS_FINAL'))
  thresholded_evolver = seq_evolution_thresh_gradient(cfg)
  thresholded_evolver.choose_best_seqs = types.MethodType(eval(cfg.get('Functions','choose_best_seqs')), thresholded_evolver)
  thresholded_evolver.thresholding_iterative(thresh, seqs_final, num_iters)
  thresholded_evolver.save_output()
  seq_selection.main(cfg)
