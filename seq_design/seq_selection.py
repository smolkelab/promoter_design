# Given a CSV with a column 'Seqs' for sequences, columns of form 'output_name' + _ + 'model_#'
# for predictions, select N sequences matching some criteria
# exception: if < N sequences available.
import sys
import os
import imp
import pandas
import random
import numpy as np
import ConfigParser

def load_from_csv(cfg):
  fn_in = os.path.expanduser(cfg.get('Files','preds_fn'))
  df_in = pandas.read_csv(fn_in)
  seqs = df_in['Seqs'].tolist()
  pred_names = [q for q in list(df_in) if q != 'Seqs']
  output_names = cfg.get('Params','OUTPUT_NAMES').strip().split(',')
  ans = []
  for o in output_names:
    o_names = [q for q in pred_names if q.split('_')[0] == o]
    o_names.sort()
    o_dat = [np.array(df_in[q]) for q in o_names]
    o_dat = np.stack(o_dat, axis = 1) # shape is n_seqs, n_models
    ans.append(o_dat)
  ans = np.stack(ans, axis = 1) # shape is n_seqs, n_outputs, n_models
  return(ans, seqs)

def get_scores(dat, cfg):
  scores = np.apply_along_axis(eval(cfg.get('Functions','merge_outputs')),1,dat) # n_seqs, m_models
  scores = np.apply_along_axis(eval(cfg.get('Functions','merge_models')),1,scores) # n_seqs
  return(scores)

# identify and reject sequences containing undesired motifs (e.g. BsaI sites)
# return boolean list of whether or not each sequence passes
def filter_sequences_by_motif(seqs, cfg):
  r_motifs = cfg.get('Params','REJECT_MOTIFS').strip().split(',')
  seq_passes = [all(r not in s for r in r_motifs) for s in seqs]
  return(seq_passes)

def filter_sequences_by_score(seqs, scores, cfg):
  thresh = float(cfg.get('Params','THRESH'))
  seqs_scores_out = [(p,q) for (p,q) in zip(seqs, scores) if q > thresh]
  return([seq for seq, _ in seqs_scores_out], np.array([score for _, score in seqs_scores_out]))

def pick_top_N(seqs, scores, N):
  # 'sort', 'argsort' put smallest first
  sorted_order = np.argsort(scores)
  assert(len(sorted_order) > 0)
  scores = scores[sorted_order]
  seqs_sorted = np.array(seqs)
  seqs_sorted = seqs_sorted[sorted_order]
  outpairs = zip(seqs_sorted.tolist(), scores.tolist())
  ans = []
  for i in range(N):
    try:
      ans.append(outpairs.pop())
    except IndexError:
      break
  return([seq for seq, _ in ans], [score for _, score in ans])

def pick_random(seqs, scores, N):
  seqs_scores = zip(seqs, scores)
  seqs_scores_out = np.random.choice(np.array(seqs_scores), size = N, replace = False).tolist()
  return([seq for seq, _ in seqs_scores_out], [score for _, score in seqs_scores_out])

def main(cfg):
  dat, seqs = load_from_csv(cfg)
  random_seed = int(cfg.get('Params','random_seed'))
  random.seed(random_seed); np.random.seed(random_seed)
  scores = get_scores(dat, cfg)

  # filter out sequences by motif
  seq_passes = filter_sequences_by_motif(seqs, cfg)
  with open(os.path.expanduser(cfg.get('Files','rejected_fn')),'w') as ff:
    ff.write('Seqs,Scores\n')
    for i,p in enumerate(seq_passes):
      if not p:
        ff.write(seqs[i] + ',' + str(scores[i]) + '\n')
  # apply the filter
  seqs = [s for (p,s) in zip(seq_passes,seqs) if p]
  scores = np.array([s for (p,s) in zip(seq_passes, scores) if p])

  # filter by strength threshold
  seqs, scores = filter_sequences_by_score(seqs, scores, cfg)
  # Select final sequences
  n_seqs = int(cfg.get('Params','NUM_SEQS_FINAL'))
  pick_top = cfg.get('Params','PICK_TOP') == 'True'
  if pick_top:
    seqs, scores = pick_top_N(seqs, scores, n_seqs)
  else:
    seqs, scores = pick_random(seqs, scores, n_seqs)
  with open(os.path.expanduser(cfg.get('Files','selected_fn')),'w') as fo:
    fo.write('Seqs,Scores\n')
    for seq, score in zip(seqs,scores):
      fo.write(seq + ',' + str(score) + '\n')

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  main(cfg)
