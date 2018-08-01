# Given an ensemble of models, evolve a random sequence to fulfill an objective.
# cf. https://blog.keras.io/how-convolutional-neural-networks-see-the-world.html

import keras
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
    scores_fx_name = cfg.get('Functions','seq_scores_keras')
    if scores_fx_name == 'None':
      self.seq_scores_fx = None
    else:
      get_seq_scores_fx = eval(scores_fx_name)
      [seq_start, seq_end] = [int(q) for q in cfg.get('Params', 'SEQ_INDICES').strip().split(',')]
      self.seq_scores_fx = get_seq_scores_fx(self.models[0], seq_start, seq_end)
    iterates = [self._get_iterate_fx_from_model(q) for q in self.models]
    def _get_mean_grad(iterates, seq_scores_fx, input):
      iterate_outputs = [q(input) for q in iterates]
      losses = [p for (p,q) in iterate_outputs]
      grads = [q for (p,q) in iterate_outputs]
      #final_grad = K.mean(K.stack(grads, axis = 0), axis = 0)
      final_grad = np.mean(np.stack(grads, axis = 0), axis = 0)
      if seq_scores_fx != None:
        seq_losses, seq_grad = seq_scores_fx(input)
        losses += seq_losses
        final_grad += seq_grad
      return(losses, final_grad)
    # given an input, get a list of losses, and the mean gradient of the input, as Numpy arrays
    self.losses_and_grads = partial(_get_mean_grad, iterates, self.seq_scores_fx)

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
    seqs = np.power(seqs, norm_power[:,np.newaxis, np.newaxis])
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
    seqs += gradient*step[:,np.newaxis, np.newaxis]
    seqs = np.clip(seqs,0.,1.) # enforce between 0 and 1
    seq_normalize = np.apply_along_axis(np.sum, 2, seqs)[...,np.newaxis]
    seqs = seqs/seq_normalize
    seqs = self._norm_bias(seqs, norm_power)
    self.seqs_iter = np.array(seqs)
    return(losses)

  def rand_range(self, shape_0, shape_1, min, max):
    r = rand(shape_0, shape_1)*(max - min)
    r = r + min
    return(r)

  def _generate_n_sequences(self, n):
    seqs = np.zeros((n,) + self.base_probs.shape)
    for i in range(self.base_probs.shape[1]):
      if(np.max(self.base_probs[:,i]) == 1.):
        seqs[:,:,i] = self.base_probs[:,i]
      else:
        probs = self.base_probs[:,i][np.newaxis,...]
        probs = np.repeat(probs, n, axis = 0)
        noise = self.rand_range(probs.shape[0], probs.shape[1], 1. - self.init_noise, 1. + self.init_noise)
        seqs[:,:,i] = probs*noise
        #noise = rand(probs.shape[0], probs.shape[1])*self.init_noise
        #seqs[:,:,i] = probs + noise
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
    seqs_iter = np.swapaxes(self.seqs, 1, 2)
    removed_pad = seqs_iter[:,seqs_iter.shape[1]-self.shift+1:,:]
    self.seqs_iter = seqs_iter[:,:seqs_iter.shape[1]-self.shift+1,:]

    print('Iteration (0): ' + str(self.curr_iters[0]))
    losses = self._update_seq_ensemble(self.map_to_key('gradient_step'), self.map_to_key('normalize_power'))

    seqs_iter = np.concatenate([self.seqs_iter, removed_pad], axis = 1)
    self.seqs = np.swapaxes(seqs_iter, 1, 2)

    self.score_tracking.append(losses)

  def basic_iterative(self, num_iters):
    for i in range(num_iters):
      self.iterate()
      self.curr_iters += 1

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

# GC filter function: use non-trainable convolution layer to get GC content and add to gradient
def get_gc_fx_from_model(model_in, seq_start, seq_end, filter_len = 20, gc_min = 0.3, gc_max = 0.7, penalty_wt = 100.):
  input_seq = model_in.input
  def get_os(start, end, input_shape):
    input_shape = list(input_shape)
    input_shape[1] = input_shape[1] - seq_start + seq_end # assume seq_start is positive, seq_end negative
    return(tuple(input_shape))
  get_output_shape = partial(get_os, seq_start, seq_end)
  input_seq_sliced = keras.layers.Lambda(lambda x: x[:,seq_start:seq_end,:], output_shape = get_output_shape)(input_seq)
  w_weights = np.zeros((filter_len,len(DNA),1))
  for (i,q) in enumerate(DNA):
    if q in 'GC':
      w_weights[:,i,:] = 1./filter_len
  b_weights_min = np.array([gc_min])
  b_weights_max = np.array([-1.*gc_max])
  gc_min = keras.layers.Conv1D(filters = 1, kernel_size = filter_len, strides=1, padding='valid', weights = [-1.*w_weights, b_weights_min], activation = 'relu')(input_seq_sliced)
  gc_min = keras.layers.Lambda(lambda x: K.min(K.abs(x), axis = 1))(gc_min)
  gc_max = keras.layers.Conv1D(filters = 1, kernel_size = filter_len, strides=1, padding='valid', weights = [w_weights, b_weights_max], activation = 'relu')(input_seq_sliced)
  gc_max = keras.layers.Lambda(lambda x: K.max(K.abs(x), axis = 1))(gc_max)
  gc_penalty = keras.layers.Add()([gc_min, gc_max])
  gc_penalty = keras.layers.Lambda(lambda x: -1.*x*penalty_wt)(gc_penalty)
  
  # ensure shape matches 'output' below
  gc_penalty = keras.layers.Lambda(lambda x: x[:,0])(gc_penalty)
  grads = K.gradients(gc_penalty, input_seq)[0]
  gc_fx = K.function([input_seq], [gc_penalty, grads])
  return(gc_fx)

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  num_iters = int(cfg.get('Params','NUM_ITERS'))
  evolver = seq_evolution_class_gradient(cfg)
  evolver.basic_iterative(num_iters)
  evolver.save_output()
  seq_selection.main(cfg)
