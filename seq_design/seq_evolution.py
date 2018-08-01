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
    self.params = self.unpack_params()
    self.dna_dict = {}
    for (i,q) in enumerate(DNA):
      p_vec = np.zeros(shape = (len(DNA),))
      p_vec[i] = 1.
      self.dna_dict[q] = p_vec
    self.score_tracking = []
    self.num_seqs = int(self.cfg.get('Params', 'NUM_SEQS'))
    self.num_variants = int(self.cfg.get('Params','NUM_VARIANTS'))
    self.random_seed = int(self.cfg.get('Params', 'RANDOM_SEED'))
    self.curr_iters = np.zeros((self.num_seqs,), dtype = 'int')
    
    seq = self.cfg.get('Params','SEQ')
    self.mutable = np.where([q not in DNA for q in seq])[0] # positions of bases that can be mutated
    self.base_probs = self.encode_seq(seq)
    
    random.seed(self.random_seed); np.random.seed(self.random_seed)
    self._populate_sequences()
    self._prepare_models()
    print('Evolver built.')

  # apply rules for extracting parameters from a config
  def unpack_params(self):
    params = {'merge_outputs': eval(self.cfg.get('Functions','merge_outputs')),
            'merge_models': eval(self.cfg.get('Functions','merge_models'))}

    # Not all of the sequence needs to be filtered (i.e. GC content-scored): the 5' homology is pretty GC-poor.
    # When calling params['seq_scores'], automatically strip out parts of the input we don't want to look at.
    seq_fx = eval(self.cfg.get('Functions','seq_scores'))
    params['seq_indices'] = [int(q) for q in self.cfg.get('Params','SEQ_INDICES').strip().split(',')]
    params['seq_scores'] = lambda x: seq_fx(x[:,params['seq_indices'][0]:params['seq_indices'][1]])
    # format for these: e.g. 50:5,30:3,20:1 for 50 cycles with 5 mutations each, 30 with 3, 20 with 1
    
    # cf. http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
    def flatten(l):
      out = []
      for item in l:
        if isinstance(item, (list, tuple)):
          out.extend(flatten(item))
        else:
          out.append(item)
      return out
    
    def decompress_pm(cfg, key):
      pm = cfg.get('Params',key)
      pm = pm.strip().split(',')
      pm = [(p.split(':')[0], p.split(':')[1]) for p in pm]
      pm = [[p]*int(q) for (p,q) in pm]
      return(flatten(pm))

    params['num_mutations'] = [int(q) for q in decompress_pm(self.cfg, 'NUM_MUTATIONS')]
    params['keep_parent'] = [q == 'True' for q in decompress_pm(self.cfg, 'KEEP_PARENT')]
    params['gradient_step'] = [float(q) for q in decompress_pm(self.cfg, 'GRADIENT_STEP')]
    params['normalize_power'] = [float(q) for q in decompress_pm(self.cfg, 'NORMALIZE_POWER')]
    assert(len(params['num_mutations'])) == int(self.cfg.get('Params','NUM_ITERS'))
    assert(len(params['keep_parent'])) == int(self.cfg.get('Params','NUM_ITERS'))
    assert(len(params['gradient_step'])) == int(self.cfg.get('Params','NUM_ITERS'))
    assert(len(params['normalize_power'])) == int(self.cfg.get('Params','NUM_ITERS'))
    return(params)

  def encode_seq(self, seq):
    ans = np.zeros(shape = (len(DNA), len(seq)))
    for i,X in enumerate(seq):
      if X not in self.dna_dict:
        wts = [float(q) for q in self.cfg.get('Params',X).strip().split(':')]
        wts = [q/float(sum(wts)) for q in wts]
        assert(len(wts) == len(DNA))
        self.dna_dict[X] = np.array(wts)
      ans[:,i] = self.dna_dict[X]
    return(ans)

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

  def de_onehot(self, seq_arr):
    x = np.argmax(seq_arr, axis = 1)
    seqs = []
    for i in x:
      seq = ''.join([DNA[q] for q in i])
      seqs.append(seq)
    return(seqs)

  def round_seqs(self, seq_arr):
    seqs = self.de_onehot(seq_arr)
    seqs_enc = [self.encode_seq(q) for q in seqs]
    return( np.stack(seqs_enc, axis = 0) )
    
  def generate_report(self):
    seqs_out = self.de_onehot(self.seqs)
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
  def mutate_seqs(self, num_mutations_arr, keep_parent_arr):
    ans = []
    for (i,s) in enumerate(self.seqs):
      keep_parent = keep_parent_arr[i]
      nv = self.num_variants
      if keep_parent:
        nv -= 1
      vars = [s] if keep_parent else []
      for i in range(nv):
        new_seq = self.mutate_one_seq(s, num_mutations_arr[i])
        assert(np.any(s != new_seq))
        vars.append(new_seq)
      ans.append(np.stack(vars, axis = 0))
    ans = np.stack(ans, axis = 0)

    return(ans) # has shape ( n_seqs, num_mutations, len(DNA), len(seq) )

# Functions relying on things to be implemented in daughter classes

  def map_to_key(self, key):
    key = self.params[key]
    ans = np.zeros(self.curr_iters.shape[0], dtype = type(key[0]))
    for (i,q) in enumerate(self.curr_iters):
      ans[i] = key[q]
    return(ans)

  def iterate(self):
    print('Iteration (0): ' + str(self.curr_iters[0]))
    # generate mutated variants of current seqs - daughter classes to determine
    vars = self.mutate_seqs(self.map_to_key('num_mutations'), self.map_to_key('keep_parent'))
    preds = self._reshaping_test_sequences(vars)
    self.seqs = self.choose_best_seqs(vars, preds)

  # Placeholder, to be overridden
  def choose_best_seqs(self, vars, preds):
    raise Exception('choose_best_seqs needs to be overridden')

  # support the most basic loop; do some number of iterations and save the results.
  def basic_iterative(self, choose_best_seqs, num_iters):
    # override choose_best_seqs: cf. https://tryolabs.com/blog/2013/07/05/run-time-method-patching-python/
    self.choose_best_seqs = types.MethodType(choose_best_seqs, self)
    # sequences were populated and models prepared by __init__
    for i in range(num_iters):
      self.iterate()
      self.curr_iters += 1

  def save_output(self):
    ans = self.generate_report()
    ans.to_csv(os.path.expanduser(self.cfg.get('Files','preds_fn')), index = False)
    scores_tracked = np.array(self.score_tracking)
    fn_scores = os.path.expanduser(self.cfg.get('Files','score_fn'))
    np.savetxt(fn_scores, scores_tracked, delimiter=',')

# Some helper functions for optimization
'''def merge_outputs_AR_useful(axis_in):
  # there should be two outputs: A (uninduced) and B (induced)
  # This gives a sigmoid mapping 1.0 -> 0.1; 1.1 -> 0.5; 1.2 -> 0.9; goal is to ensure Pred_B is useful
  # while otherwise optimizing activation ratio (B - A)
  k1 = math.log(1./9.); k2 = math.log(9.)
  m = (k2-k1)/0.2; b = k1 - m
  b_t = m*axis_in[1] + b
  b_sig = math.exp(b_t)/(math.exp(b_t)+1.)
  return( b_sig + axis_in[1] - axis_in[0])'''
def merge_outputs_AR(axis_in):
  return(axis_in[1] - axis_in[0])

# just what it looks like - marginally more readable than a lambda, maybe
def mean_minus_sd(axis_in):
  return( np.mean(axis_in) - np.std(axis_in) )

# as above.
def get_induced(axis_in):
  return(axis_in[1])
 
# Is there a region of the sequence of length 'filter_len' with GC content < gc_min or > gc_max?
# if yes, apply a penalty 'present_penalty' to the score.
# Idea is that any such regions will be strongly selected against.
def gc_filter(seq_arr, filter_len = 20, gc_min = 0.3, gc_max = 0.7, present_penalty = 10.): # seq_arr has shape (len(DNA), len(seq))
  is_gc = np.apply_along_axis(np.sum, 0, seq_arr[[i for (i,q) in enumerate(DNA) if q in 'GC'],:])
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

def greedy_choose_best_seqs(self, vars, preds):
  # join multiple outputs from models into a single score; 'preds' is shape (n_seqs, n_mutations, n_outputs, n_models)
  preds = np.apply_along_axis(self.params['merge_outputs'], 2, preds) # preds is now shape (n_seqs, n_mutations, n_models)
  model_scores = np.apply_along_axis(self.params['merge_models'], 2, preds) # preds is now shape(n_seqs, n_mutations)
  seq_scores = np.zeros(shape = model_scores.shape)
  for i in range(preds.shape[0]):
    for j in range(preds.shape[1]):
      seq_scores[i,j] = self.params['seq_scores'](vars[i,j,...])
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

# Assume this is a GPD-like task, and use 'greedy_choose_best_seqs'
if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  num_iters = int(cfg.get('Params','NUM_ITERS'))
  evolver = seq_evolution_class(cfg)
  evolver.basic_iterative(greedy_choose_best_seqs, num_iters)
  evolver.save_output()
  seq_selection.main(cfg)
