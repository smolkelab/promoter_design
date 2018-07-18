# Generate interesting sequences by directly generating large sets and picking those that pass
# some threshold
import sys
import os
import ConfigParser
import numpy as np
import random
import pandas
import seq_evolution
import time
import seq_selection

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  random_seed = int(cfg.get('Params','random_seed'))
  random.seed(random_seed); np.random.seed(random_seed)
  screener = seq_evolution.seq_evolution_class(cfg)
  print('screener ready')
  params = seq_evolution.unpack_params(cfg)
  thresh = float(cfg.get('Params','THRESH'))
  seqs_needed = int(cfg.get('Params','NUM_SEQS_FINAL'))
  filter_fx = eval(cfg.get('Functions','seq_scores'))
  filter_fx_bool = lambda x: float(filter_fx(x)) == 0.
  seqs_passing = []
  curr_iter = 0
  while len(seqs_passing) < seqs_needed:
    t = time.time()
    print(len(seqs_passing))
    preds = screener._test_sequences(screener.seqs)
    preds = preds[:,np.newaxis,...] # dummy axis for "variants," which we aren't using here
    print(preds.shape)
    preds = np.apply_along_axis(params['merge_outputs'],2,preds)
    print(preds.shape)
    scores = np.apply_along_axis(params['merge_models'],2,preds).squeeze() # 1D array of scores
    print(np.max(scores))
    seqs_passing.extend([p for (p,q) in zip(screener.seqs, scores) if q > thresh and filter_fx_bool(p)])
    screener._populate_sequences()
    print('Iter ' + str(curr_iter) + ': ' + str(time.time() - t))
  seqs_passing = seqs_passing[:seqs_needed]
  screener.seqs = np.array(seqs_passing)
  ans = screener.generate_report()
  ans.to_csv(os.path.expanduser(cfg.get('Files','preds_fn')), index = False)
  
  seq_selection.main(cfg)
