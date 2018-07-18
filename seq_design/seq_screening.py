# Generate interesting sequences by directly generating large sets and picking those that pass
# some threshold
import sys
import os
import ConfigParser
import numpy as np
import pandas
from seq_evolution import seq_evolution

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  screener = seq_evolution(cfg)
  print('screener ready')
  params = {'num_mutations': int(cfg.get('Params','NUM_MUTATIONS')),
            'keep_parent': cfg.get('Params','KEEP_PARENT') == 'True',
            'merge_outputs': eval(cfg.get('Functions','merge_outputs')),
            'merge_models': eval(cfg.get('Functions','merge_models'))}
  thresh = float(cfg.get('Params','THRESH'))
  seqs_needed = int(cfg.get('Params','NUM_SEQS_FINAL'))
  seqs_passing = []
  while len(seqs_passing) < seqs_needed:
    print(len(seqs_passing))
    preds = screener._test_sequences(screener.seqs)
    preds = preds[:,np.newaxis,...] # dummy axis for "variants," which we aren't using here
    preds = np.apply_along_axis(params['merge_outputs'],2,preds)
    scores = np.apply_along_axis(params['merge_models'],2,preds).squeeze() # 1D array of scores
    print(np.max(scores))
    seqs_passing.extend([p for (p,q) in zip(screener.seqs, scores) if q > thresh])
    screener._populate_sequences()
  seqs_passing = seqs_passing[:seqs_needed]
  screener.seqs = np.array(seqs_passing)
  ans = screener.generate_report()
  ans.to_csv(os.path.expanduser(cfg.get('Files','preds_fn')), index = False)
