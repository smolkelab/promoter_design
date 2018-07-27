# Given an ensemble of models, evolve a random sequence to fulfill an objective.
# cf. https://blog.keras.io/how-convolutional-neural-networks-see-the-world.html

from keras import backend as K
from functools import partial

import sys
import os
import pandas
import numpy as np
import random
from numpy.random import choice, rand
import ConfigParser
import seq_evolution
import seq_selection

def debug_dists(seq1, seq2):
  return np.sqrt(np.sum((seq1 - seq2)**2))

DNA = seq_evolution.DNA

class seq_evolution_class_gradient(seq_evolution.seq_evolution_class):

  def __init__(self, cfg):
    # prepare the models, populate the sequences, et cetera
    self.init_noise = float(cfg.get('Params','INIT_NOISE'))
    super(seq_evolution_class_gradient, self).__init__(cfg)
    self.loss_tensor_fx = self._get_loss_tensor_fx(cfg) # given a model, get an output tensor
    iterates = [self._get_iterate_fx_from_model(q) for q in self.models]
    def _get_mean_grad(iterates, input):
      iterate_outputs = [q(input) for q in iterates]
      losses = [p for (p,q) in iterate_outputs]
      grads = [q for (p,q) in iterate_outputs]
      #final_grad = K.mean(K.stack(grads, axis = 0), axis = 0)
      final_grad = np.mean(np.stack(grads, axis = 0), axis = 0)
      return(losses, final_grad)
    # given an input, get a list of losses, and the mean gradient of the input, as Numpy arrays
    self.losses_and_grads = partial(_get_mean_grad, iterates)

    # For gradient updates - keep from "updating" immutable bases
    self.mutable_mask = np.zeros(self.base_probs.shape[1])
    self.mutable_mask[self.mutable] = 1.
    self.mutable_mask = self.mutable_mask[np.newaxis,:self.mutable_mask.shape[0]-self.shift+1, np.newaxis]

  def _get_loss_tensor_fx(self, cfg):
    def loss_wrapper(merge_outputs_keras, model_in):
      model_output = model_in.output
      loss = merge_outputs_keras(model_output)
      return(loss)
    merge_outputs_keras = eval(cfg.get('Functions','merge_outputs_keras')) # needs to be written in terms of the backend
    loss_tensor_fx = partial(loss_wrapper, merge_outputs_keras)
    return(loss_tensor_fx) # takes a model, returns a loss tensor

  # given a model, get a function taking an input tensor and returning loss and gradient tensors.
  def _get_iterate_fx_from_model(self, model_in):
    input_seq = model_in.input
    loss = self.loss_tensor_fx(model_in)
    grads = K.gradients(loss, input_seq)[0] # gradient at the input
    # following example, normalize gradient
    # cf. https://github.com/keras-team/keras/blob/master/examples/conv_filter_visualization.py
    normalize = lambda x: x / (K.sqrt(K.mean(K.square(x))) + K.epsilon())
    grads = normalize(grads)
    iterate = K.function([input_seq], [loss, grads])
    return(iterate)

  # raise all values to some power, and renormalize -
  # this encourages one value to dominate, approximating one-hot encoding.
  def _norm_bias(self, seqs, norm_power):
    print(np.min(seqs))
    print(np.max(seqs))
    seqs = np.power(seqs, norm_power)
    seq_normalize = np.apply_along_axis(np.sum, 2, seqs)[...,np.newaxis]
    seqs = seqs/seq_normalize
    return(seqs)

  # given an ensemble of models and a list of sequences, update with the mean of gradients
  def _update_seq_ensemble(self, step, norm_power):
    seqs = self.seqs_iter
    losses, gradient = self.losses_and_grads([seqs])
    losses = np.apply_along_axis(np.mean,0,np.stack(losses))
    print(losses)
    gradient = np.multiply(gradient, self.mutable_mask)
    seqs += gradient*step
    seqs = np.clip(seqs,0.,1.) # enforce between 0 and 1
    seq_normalize = np.apply_along_axis(np.sum, 2, seqs)[...,np.newaxis]
    seqs = seqs/seq_normalize
    seqs = self._norm_bias(seqs, norm_power)
	self.seqs_iter = np.array(seqs)
    return(losses)

  def _generate_n_sequences(self, n):
    seqs = np.zeros((n,) + self.base_probs.shape)
    for i in range(self.base_probs.shape[1]):
      if(np.max(self.base_probs[:,i]) == 1.):
        seqs[:,:,i] = self.base_probs[:,i]
      else:
        probs = self.base_probs[:,i][np.newaxis,...]
        probs = np.repeat(probs, self.num_seqs, axis = 0)
        noise = rand(probs.shape[0], probs.shape[1])*self.init_noise
        seqs[:,:,i] = probs + noise
    normalize_arr = np.apply_along_axis(np.sum,1, seqs)[:,np.newaxis,:]
    seqs = seqs/normalize_arr
    return(seqs)

  def _populate_sequences(self):
    seqs = self._generate_n_sequences(self.num_seqs)
    self.seqs = seqs

  # In seq_evolve_to_threshold, we want to generate replacement sequences one at a time.
  def _populate_one_sequence(self):
    return( self._generate_n_sequences(1) )

  # override method in parent class; update via gradient
  def iterate(self): # , params, iter_idx
    print('Iteration (0): ' + str(self.curr_iters[0]))
    losses = self._update_seq_ensemble(self.map_to_key('gradient_step'), self.map_to_key('normalize_power')
    self.score_tracking.append(losses)

  def basic_iterative(self, num_iters):
    # sequences were populated and models prepared by __init__
    seqs_iter = np.swapaxes(self.seqs, 1, 2)
    removed_pad = seqs_iter[:,seqs_iter.shape[1]-self.shift+1:,:]
    self.seqs_iter = seqs_iter[:,0:seqs_iter.shape[1]-self.shift+1,:]
	
    for i in range(num_iters):
      self.iterate()
	  self.curr_iters += 1

    seqs_iter = np.concatenate([self.seqs_iter, removed_pad], axis = 1)
    self.seqs = np.swapaxes(seqs_iter, 1, 2)

  def generate_report(self):
    seqs_out = self.de_onehot(self.seqs)
    ans = {'Seqs': seqs_out}
    orig_preds = self._test_sequences(self.seqs)
	self.seqs = self.round_seqs(self.seqs)
    final_preds = self._test_sequences(self.seqs)
    print(np.mean((final_preds - orig_preds)[:,1,:],axis=1))
    print(final_preds.shape)
    for i,p in enumerate(self.output_names):
      for j in range(final_preds.shape[-1]):
        this_pred = final_preds[:,i,j]
        this_name = p + '_' + str(j)
        ans[this_name] = this_pred
    return(pandas.DataFrame(ans))

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  num_iters = int(cfg.get('Params','NUM_ITERS'))
  evolver = seq_evolution_class_gradient(cfg)
  evolver.basic_iterative(num_iters)
  evolver.save_output()
  seq_selection.main(cfg)