# Overall goal: Given a starting pattern and a set of models, generate new sequences
# In this file - generally useful functions, in the 'seq_evolution' class. Implementation-specific stuff elsewhere.
import sys
import os
import pandas
import numpy as np
import random
from numpy.random import choice
import ConfigParser
import imp
import types
import math
import seq_evolution
import seq_selection

DNA = ['A','C','G','T']

class seq_evolution_class(object):

  def __init__(self, cfg):
    print('Building evolver...')
    self.cfg = cfg
    self.dna_dict = {}
    for (i,q) in enumerate(DNA):
      p_vec = np.zeros(shape = (len(DNA),))
      p_vec[i] = 1.
      self.dna_dict[q] = p_vec
    self.score_tracking = []
    self.num_seqs = int(self.cfg.get('Params', 'NUM_SEQS'))
    self.num_variants = int(self.cfg.get('Params','NUM_VARIANTS'))
    self.random_seed = int(self.cfg.get('Params', 'RANDOM_SEED'))
    self._get_base_probs()
    random.seed(self.random_seed); np.random.seed(self.random_seed)
    self._populate_sequences()
    self._prepare_models()
    print('Evolver built.')

  def _get_base_probs(self):
    seq = self.cfg.get('Params','SEQ')
    self.mutable = np.where([q not in DNA for q in seq])[0] # positions of bases that can be mutated
    self.base_probs = np.zeros(shape = (len(DNA), len(seq)))
    for i,X in enumerate(seq):
      if X not in self.dna_dict:
        wts = [float(q) for q in self.cfg.get('Params',X).strip().split(':')]
        wts = [q/float(sum(wts)) for q in wts]
        assert(len(wts) == len(DNA))
        self.dna_dict[X] = np.array(wts)
      self.base_probs[:,i] = self.dna_dict[X]

  def _populate_sequences(self):
    seqs = np.zeros((self.num_seqs,) + self.base_probs.shape)
    idx0 = np.arange(self.num_seqs)
    for i in range(self.base_probs.shape[1]):
      if(np.max(self.base_probs[:,i]) == 1.):
        seqs[:,:,i] = self.base_probs[:,i]
      else:
        idx1 = choice(self.base_probs.shape[0], size = self.num_seqs, p = self.base_probs[:,i])
        seqs[idx0,idx1,i] = 1.
    self.seqs = seqs

  def _populate_one_sequence(self):
    s = np.zeros(shape = self.base_probs.shape)
    for i in range(self.base_probs.shape[1]):
      s[choice(self.base_probs.shape[0], size = None, p = self.base_probs[:,i]), i] = 1.
    return(s)

  def _prepare_models(self):
    # get the model-building code and the filenames of weights to use
    do_model = imp.load_source('do_model', os.path.expanduser(self.cfg.get('Files','model_fn')))
    self.shift = do_model.SHIFT
    wts_dir = os.path.expanduser(self.cfg.get('Dirs','weights_dir'))
    self.weight_fns = [os.path.join(wts_dir, q) for q in os.listdir(wts_dir) if self.cfg.get('Params','WTS_EXT') in q]
    self.weight_fns.sort()
    # instantiate the models; create 'dat_to_use' to pass in shape information in a format that 'do_model' expects
    self.models = []
    dat_to_use = [ [np.zeros(1,), np.zeros(1,)],[np.zeros(1,), np.zeros(1,)],[np.zeros(1,), np.zeros(1,)], self.base_probs.shape[1] ] # last entry is length of sequence
    self.output_names = self.cfg.get('Params','OUTPUT_NAMES').strip().split(',')
    self.num_outputs = len(self.output_names)
    for i in self.weight_fns:
      this_model = do_model.do_model(dat_to_use, self.num_outputs, train = False)
      this_model.load_weights(i)
      self.models.append(this_model)

  # input is one-hot array with shape (num_seqs, len(DNA), len(seq); output is float array with shape (num_seqs, num_outputs, len(self.models))
  # NB in practice, an array with multiple new variants of each sequence will need to be carefully reshaped;
  # see '_reshaping_test_sequences', below
  def _test_sequences(self, seq_arr):
    seq_arr = np.swapaxes(seq_arr, 1,2)
    preds = [q.predict(seq_arr[:,0:seq_arr.shape[1]-self.shift+1,:]) for q in self.models]
    preds = np.stack(preds, axis = -1)
    return(preds)

  # seq_arr has arbitrary shape, but the last two dimensions are to be left alone,
  # and the rest are to be merged to one. Run _test_sequences, and undo the reshaping.
  def _reshaping_test_sequences(self, seq_arr):
    shape_in = seq_arr.shape
    seq_arr_reshaped = np.reshape(seq_arr, (np.prod(shape_in[:-2]), shape_in[-2], shape_in[-1]))
    preds_reshaped = self._test_sequences(seq_arr_reshaped)
    new_shape = shape_in[:-2] + preds_reshaped.shape[-2:]
    preds = np.reshape(preds_reshaped, new_shape)
    return(preds)

  def generate_report(self):
    def de_onehot(seq_arr):
      print(seq_arr.shape)
      x = np.argmax(seq_arr, axis = 1)
      print(x.shape)
      seqs = []
      for i in x:
        seq = ''.join([DNA[q] for q in i])
        seqs.append(seq)
      return(seqs)
    seqs_out = de_onehot(self.seqs)
    ans = {}
    ans['Seqs'] = seqs_out
    final_preds = self._test_sequences(self.seqs)
    for i,p in enumerate(self.output_names):
      for j in range(final_preds.shape[-1]):
        ans[p + '_' + str(j)] = final_preds[:,i,j]
    return(pandas.DataFrame(ans))

  def mutate_one_seq(self, seq, num_mutations):
    seq_orig = np.copy(seq)
    seq = np.copy(seq)
    bases_to_change = np.random.choice(self.mutable, size = num_mutations, replace = False)
    for b in bases_to_change:
      new_base_idxes = np.where(seq[:,b] == 0.)[0]
      new_base_probs = self.base_probs[new_base_idxes,b]
      new_base_probs = new_base_probs/np.sum(new_base_probs)
      new_base_idx = np.random.choice(new_base_idxes,1,p = new_base_probs)
      # old approach, doesn't take base_probs into account
      #new_base_idx = np.random.choice(np.where(seq[:,b] == 0.)[0],1)
      new_base = np.zeros(shape = (len(DNA),))
      new_base[new_base_idx] = 1.
      seq[:,b] = new_base
    assert(np.any(seq_orig != seq))
    return(seq)

  # num_mutations is number of mutations each new sequence should have.
  # keep_parent is boolean; if 'True', one of the variants will be the original sequence.
  def mutate_seqs(self, num_mutations, keep_parent):
    nv = self.num_variants
    if keep_parent:
      nv -= 1
    ans = []
    for s in self.seqs:
      vars = [s] if keep_parent else []
      for i in range(nv):
        new_seq = self.mutate_one_seq(s, num_mutations)
        assert(np.any(s != new_seq))
        vars.append(new_seq)
      ans.append(np.stack(vars, axis = 0))
    ans = np.stack(ans, axis = 0)

    return(ans) # has shape ( n_seqs, num_mutations, len(DNA), len(seq) )

# Functions relying on things to be implemented in daughter classes

  def iterate(self, params, iter_idx):
    print('Iteration: ' + str(iter_idx))
    # generate mutated variants of current seqs - daughter classes to determine
    vars = self.mutate_seqs(params['num_mutations'][iter_idx], params['keep_parent'][iter_idx])
    preds = self._reshaping_test_sequences(vars)
    self.seqs = self.choose_best_seqs(vars, preds, params)

  # Placeholder, to be overridden
  def choose_best_seqs(self, vars, preds, params):
    raise Exception('choose_best_seqs needs to be overridden')

  # support the most basic loop; do some number of iterations and save the results.
  def basic_iterative(self, choose_best_seqs, num_iters, params):
    # override choose_best_seqs: cf. https://tryolabs.com/blog/2013/07/05/run-time-method-patching-python/
    self.choose_best_seqs = types.MethodType(choose_best_seqs, self)
    # sequences were populated and models prepared by __init__
    for i in range(num_iters):
      self.iterate(params, i)
    return( self.generate_report() )

# Some helper functions for optimization
def merge_outputs_AR_useful(axis_in):
  # there should be two outputs: A (uninduced) and B (induced)
  # This gives a sigmoid mapping 1.0 -> 0.1; 1.1 -> 0.5; 1.2 -> 0.9; goal is to ensure Pred_B is useful
  # while otherwise optimizing activation ratio (B - A)
  k1 = math.log(1./9.); k2 = math.log(9.)
  m = (k2-k1)/0.2; b = k1 - m
  b_t = m*axis_in[1] + b
  b_sig = math.exp(b_t)/(math.exp(b_t)+1.)
  return( b_sig + axis_in[1] - axis_in[0])

# just what it looks like - marginally more readable than a lambda, maybe
def mean_minus_sd(axis_in):
  return( np.mean(axis_in) - np.std(axis_in) )

# as above.
def get_induced(axis_in):
  return(axis_in[1])
 
# Is there a region of the sequence of length 'filter_len' with GC content < gc_min or > gc_max?
# if yes, apply a penalty 'present_penalty' to the score.
# Idea is that any such regions will be strongly selected against.
def gc_filter(seq_arr, filter_len = 20, gc_min = 0.2, gc_max = 0.8, present_penalty = 10.): # seq_arr has shape (len(DNA), len(seq))
  is_gc = np.apply_along_axis(np.sum, 0, seq_arr[[q in 'GC' for q in DNA],:])
  is_at = np.logical_not(is_gc).astype('float')
  most_gc = _max_window_score(is_gc, filter_len)
  most_at = _max_window_score(is_at, filter_len)
  if most_gc/float(filter_len) > gc_max or most_at/float(filter_len) > (1. - gc_min):
    return(-1.*present_penalty)
  return(0.)

# find the stretch with the most '1.' in an array, given a window size
def _max_window_score(arr_in, w_s):
  cs = np.pad(np.cumsum(arr_in), (1,0), mode = 'constant')
  scores = cs[w_s:] - cs[:-w_s] # cumsum of windows - i.e., number of '1.' in each window assuming all 0. or 1.
  return(np.max(scores))

def greedy_choose_best_seqs(self, vars, preds, params):
  # join multiple outputs from models into a single score; 'preds' is shape (n_seqs, n_mutations, n_outputs, n_models)
  preds = np.apply_along_axis(params['merge_outputs'], 2, preds) # preds is now shape (n_seqs, n_mutations, n_models)
  model_scores = np.apply_along_axis(params['merge_models'], 2, preds) # preds is now shape(n_seqs, n_mutations)
  seq_scores = np.zeros(shape = model_scores.shape)
  for i in range(preds.shape[0]):
    for j in range(preds.shape[1]):
      seq_scores[i,j] = params['seq_scores'](vars[i,j,...])
  scores = model_scores + seq_scores
  best_ids = np.apply_along_axis(np.argmax, 1, scores) # shape (n_seqs,); ID of best variant of each sequence
  # vars has shape ( n_seqs, num_mutations, len(DNA), len(seq) )
  new_seqs = []
  best_scores = []
  for i in range(best_ids.shape[0]):
    new_seqs.append(vars[i,best_ids[i],...])
    best_scores.append(scores[i,best_ids[i]])
  self.score_tracking.append(best_scores)
  return( np.stack(new_seqs, axis = 0) )

# apply rules for extracting parameters from a config
def unpack_params(cfg):
  params = {'merge_outputs': eval(cfg.get('Functions','merge_outputs')),
            'merge_models': eval(cfg.get('Functions','merge_models')),
            'seq_scores': eval(cfg.get('Functions','seq_scores'))}
  # format for these: e.g. 50:5,30:3,20:1 for 50 cycles with 5 mutations each, 30 with 3, 20 with 1
  def decompress_pm(cfg, key):
    pm = cfg.get('Params',key)
    pm = pm.strip().split(',')
    pm = [(p.split(':')[0], p.split(':')[1]) for p in pm]
    pm = [[p]*int(q) for (p,q) in pm]
    return(flatten(pm))
  params['num_mutations'] = [int(q) for q in decompress_pm(cfg, 'NUM_MUTATIONS')]
  params['keep_parent'] = [q == 'True' for q in decompress_pm(cfg, 'KEEP_PARENT')]
  assert(len(params['num_mutations'])) == int(cfg.get('Params','NUM_ITERS'))
  assert(len(params['keep_parent'])) == int(cfg.get('Params','NUM_ITERS'))
  return(params)

# cf. http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
def flatten(l):
  out = []
  for item in l:
    if isinstance(item, (list, tuple)):
      out.extend(flatten(item))
    else:
      out.append(item)
  return out
  
# Assume this is a GPD-like task, and use 'greedy_choose_best_seqs'
if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  random_seed = int(cfg.get('Params','random_seed'))
  random.seed(random_seed); np.random.seed(random_seed)
  num_iters = int(cfg.get('Params','NUM_ITERS'))
  evolver = seq_evolution_class(cfg)
  params = unpack_params(cfg)
  ans = evolver.basic_iterative(greedy_choose_best_seqs, num_iters, params)
  ans.to_csv(os.path.expanduser(cfg.get('Files','preds_fn')), index = False)
  scores_tracked = np.array(evolver.score_tracking)
  fn_scores = os.path.expanduser(cfg.get('Files','score_fn'))
  np.savetxt(fn_scores, scores_tracked, delimiter=',')
  
  seq_selection.main(cfg)
