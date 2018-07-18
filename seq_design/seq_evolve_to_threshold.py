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

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  random_seed = int(cfg.get('Params','random_seed'))
  random.seed(random_seed); np.random.seed(random_seed)
  num_iters = int(cfg.get('Params','NUM_ITERS'))
  thresh = float(cfg.get('Params','THRESH'))
  seqs_final = int(cfg.get('Params','NUM_SEQS_FINAL'))
  params = seq_evolution.unpack_params(cfg)
  thresholded_evolver = seq_evolution.seq_evolution(cfg)
  thresholded_evolver.choose_best_seqs = types.MethodType(eval(cfg.get('Functions','choose_best_seqs')), thresholded_evolver)
  seqs_out = []
  for i in range(num_iters):
    # test the current set of sequences - add the good ones to seq_out
    model_scores = params['merge_models'](params['merge_outputs'](thresholded_evolver._test_sequences(thresholded_evolver.seqs)))
    seq_scores = np.apply_along_axis(params['seq_scores'], 0, thresholded_evolver.seqs)
    scores = model_scores + seq_scores
    thresholded_evolver.score_tracking.append(scores)
    done_pos = scores > thresh
    seq_out.extend([q for q in thresholded_evolver.seqs[done_pos]])
    for j, q in enumerate(done_pos):
      if q:
        thresholded_evolver.seqs[q] = thresholded_evolver._populate_one_sequence()
    if len(seq_out) >= seqs_final:
      break
    
    # iterate!
    thresholded_evolver.iterate(params, i)
  
  seqs_out = seqs_out[:seqs_final]
  thresholded_evolver.seqs = np.array(seqs_out)
  ans = thresholded_evolver.generate_report()
  ans.to_csv(os.path.expanduser(cfg.get('Files','preds_fn')), index = False)
  scores_tracked = np.array(evolver.score_tracking)
  fn_scores = os.path.expanduser(cfg.get('Files','score_fn'))
  np.savetxt(fn_scores, scores_tracked, delimiter=',')
