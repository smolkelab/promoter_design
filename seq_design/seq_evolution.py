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

DNA = ['A','C','G','T']

class seq_evolution(object):

  def __init__(self, cfg):
    self.cfg = cfg
    self.dna_dict = {}
    for (i,q) in enumerate(DNA):
      p_vec = np.zeros(shape = (len(DNA),))
      p_vec[i] = 1.
      self.dna_dict[q] = p_vec
    self.num_seqs = int(self.cfg.get('Params', 'NUM_SEQS'))
    self.num_variants = int(self.cfg.get('Params','NUM_VARIANTS'))
    self.random_seed = int(self.cfg.get('Params', 'RANDOM_SEED'))
    self._get_base_probs()
    random.seed(self.random_seed); np.random.seed(self.random_seed)
    self._populate_sequences()
    self._prepare_models()
  
  def _get_base_probs(self):
    seq = cfg.get('Params','SEQ')
    self.mutable = np.where([q not in DNA for q in seq])[0] # positions of bases that can be mutated
    self.base_probs = np.zeros(shape = (len(DNA), len(seq)))
    for i,X in enumerate(seq):
      if X not in self.dna_dict:
        wts = [float(q) for q in cfg.get('Params',X).strip().split(':')]
        wts = [q/float(sum(wts)) for q in wts]
        assert(len(wts) == len(DNA))
        self.dna_dict[X] = np.array(wts)
      self.base_probs[i,:] = self.dna_dict[X]
  
  def _populate_sequences(self):
    seqs = []
    for q in range(self.num_seqs):
      s = np.zeros(shape = self.base_probs.shape)
      for i in range(self.base_probs.shape[1]):
        s[choice(self.base_probs.shape[0], size = None, p = self.base_probs[i,:]), i] = 1.
      seqs.append(s)
    self.seqs = np.stack(seqs, axis = 0)

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
      this_model = do_model(dat_to_use, num_outputs, train = False)
      this_model.load_weights(i)
      self.models.append(this_model)

  # input is one-hot array with shape (num_seqs, len(DNA), len(seq); output is float array with shape (num_seqs, num_outputs, len(self.models))
  # NB in practice, an array with multiple new variants of each sequence will need to be carefully reshaped;
  # see '_reshaping_test_sequences', below
  def _test_sequences(self, seq_arr):
    preds = [q.predict(seq_arr[:,0:seq_arr.shape[1]-self.shift+1,:]) for q in self.models]
    preds = np.stack(preds, axis = -1)
    return(preds)

  # seq_arr has arbitrary shape, but the last two dimensions are to be left alone,
  # and the rest are to be merged to one. Run _test_sequences, and undo the reshaping.
  def _reshaping_test_sequences(self, seq_arr):
    shape_in = seq_arr.shape
    seq_arr_reshaped = np.reshape(seq_arr, (np.prod(shape_in[:-2]), shape_in[-2], shape_in[-1]))
    preds_reshaped = self._test_sequences(seq_arr_reshaped)
    preds = np.reshape(preds_reshaped, (shape_in[:-2], preds_reshaped.shape[-2], preds_reshaped.shape[-1]))
    return(preds) 
    
  def generate_report(self):
    def de_onehot(seq_arr):
      x = np.argmax(seq_arr, axis = 2)
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

  def mutate_one_seq(self, seq, num_mutations)
    bases_to_change = np.random.choice(self.mutable, size = num_mutations, replace = False)
    for b in bases_to_change:
      new_base_idx = np.random.choice(np.where(seq[b,:] == 0.)[0],1)
      new_base = np.zeros(shape = (len(DNA),))
      new_base[new_base_idx] = 1.
      seq[b,:] = new_base
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
        vars.append(self.mutate_one_seq(s, num_mutations))
      ans.append(np.stack(vars, axis = 0))
    return( np.stack(ans, axis = 0) ) # has shape ( n_seqs, num_mutations, len(DNA), len(seq) )

# Functions relying on things to be implemented in daughter classes
    
  def iterate(self, **kwargs):
    # generate mutated variants of current seqs - daughter classes to determine 
    vars = self.mutate_seqs(kwargs['num_mutations'], kwargs['keep_parent'])
    preds = self._reshaping_test_sequences(vars)
    self.seqs = self.choose_best_seqs(vars, preds, kwargs)

  # Placeholder, to be overridden
  def choose_best_seqs(self, vars, preds, **kwargs):
    raise Exception('choose_best_seqs needs to be overridden')

  # support the most basic loop; do some number of iterations and save the results.
  def basic_iterative(self, choose_best_seqs, num_iters, **kwargs):
    # override choose_best_seqs: cf. https://tryolabs.com/blog/2013/07/05/run-time-method-patching-python/
    self.choose_best_seqs = types.MethodType(choose_best_seqs, self)
    # sequences were populated and models prepared by __init__
    for i in range(num_iters):
      self.iterate(kwargs)
    return( self.generate_report() )

def greedy_choose_best_seqs(self, vars, preds, **kwargs):
  # join multiple outputs from models into a single score; 'preds' is shape (n_seqs, n_mutations, n_outputs, n_models)
  preds = np.apply_along_axis(kwargs['merge_outputs'], 2, preds) # preds is now shape (n_seqs, n_models)
  scores = np.apply_along_axis(kwargs['merge_models'], 2, preds) # just take the average of model outputs
  best_ids = np.apply_along_axis(np.argmax, 1, scores) # shape (n_seqs,); ID of best variant of each sequence
  # vars has shape ( n_seqs, num_mutations, len(DNA), len(seq) )
  new_seqs = []
  for i in range(best_ids.shape[0]):
    new_seqs.append(vars[i,best_ids[i],...])
  return( np.stack(new_seqs, axis = 0) )

  
# Assume this is a GPD-like task, and use 'greedy_choose_best_seqs'
if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  num_iters = int(cfg.get('Params','NUM_ITERS'))
  evolver = seq_evolution(cfg)
  kwargs = {'num_mutations': 5, 'keep_parent': True, 'merge_outputs': np.mean, 'merge_models': np.mean}
  ans = evolver.basic_iterative(greedy_choose_best_seqs, num_iters, kwargs)
  ans.to_csv(os.path.expanduser(cfg.get('Files','output_fn')), index = False)