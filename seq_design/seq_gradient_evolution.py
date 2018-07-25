# Given an ensemble of models, evolve a random sequence to fulfill an objective.
# cf. https://blog.keras.io/how-convolutional-neural-networks-see-the-world.html

from keras import backend as K
from functools import partial

import sys
import os
import pandas
import numpy as np
import random
from numpy.random import choice, random
import ConfigParser
import seq_evolution
import seq_selection

def debug_dists(seq1, seq2):
  return np.sqrt(np.sum((seq1 - seq2)**2))

DNA = seq_evolution.DNA

class seq_evolution_class_gradient(seq_evolution.seq_evolution_class):

  def __init__(self, cfg):
    # prepare the models, populate the sequences, et cetera
    super(seq_evolution_class_gradient, self).__init__(cfg)
    self.init_noise = float(cfg.get('Params','INIT_NOISE'))
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

  # given an ensemble of models and a list of sequences, update with the mean of gradients
  def _update_seq_ensemble(self, seqs, step):
    losses, gradient = self.losses_and_grads([seqs])
    print(np.apply_along_axis(np.mean,0,np.stack(losses)))
    gradient = np.multiply(gradient, self.mutable_mask)
    seqs += gradient*step
    seq_normalize = np.apply_along_axis(np.sum, 2, seqs)[...,np.newaxis]
    seqs = seqs/seq_normalize
    return(np.array(seqs))

  def _populate_sequences(self):
    seqs = np.zeros((self.num_seqs,) + self.base_probs.shape)
    idx0 = np.arange(self.num_seqs)
    for i in range(self.base_probs.shape[1]):
      if(np.max(self.base_probs[:,i]) == 1.):
        seqs[:,:,i] = self.base_probs[:,i]
      else:
        #idx1 = choice(self.base_probs.shape[0], size = self.num_seqs, p = self.base_probs[:,i])
        #seqs[idx0,idx1,i] = 1.
        seqs[idx0, idx1,:] = self.base_probs[:,idx1] + rand(self.base_probs.shape[0])*self.init_noise
        seqs[idx0, idx1,:] /= np.sum(seqs[idx0, idx1,:])
    self.seqs = seqs
    
  # override method in parent class; update via gradient
  def iterate(self, params, iter_idx):
    self.seqs_iter = self._update_seq_ensemble(self.seqs_iter, float(params['gradient_step'][iter_idx]))

  def basic_iterative(self, num_iters, params):
    # sequences were populated and models prepared by __init__
    # unlike in base class, self.iterate will operate on Keras tensors,
    # not Numpy arrays
    #for i in range(self.seqs.shape[0]):
    #  print(debug_dists(self.seqs[0], self.seqs[i]))

    seqs_iter = np.swapaxes(self.seqs, 1, 2)
    removed_pad = seqs_iter[:,seqs_iter.shape[1]-self.shift+1:,:]
    self.seqs_iter = seqs_iter[:,0:seqs_iter.shape[1]-self.shift+1,:]
    for i in range(num_iters):
      self.iterate(params, i)

    seqs_iter = np.concatenate([self.seqs_iter, removed_pad], axis = 1)
    self.seqs = np.swapaxes(seqs_iter, 1, 2)
    return( self.generate_report() )

  def generate_report(self):
    def de_onehot(seq_arr):
      #print(seq_arr.shape)
      x = np.argmax(seq_arr, axis = 1)
      #print(x.shape)
      seqs = []
      for i in x:
        seq = ''.join([DNA[q] for q in i])
        seqs.append(seq)
      return(seqs)
    seqs_out = de_onehot(self.seqs)
    #print(seqs_out)
    ans = {'Seqs': seqs_out}

    #print(self._test_sequences(self.seqs))
    #for i in range(self.seqs.shape[0]):
    #  print(debug_dists(self.seqs[0], self.seqs[i]))

    for (i,s) in enumerate(self.seqs):
      x = np.argmax(s, axis = 0)
      seq_onehot = np.zeros(self.base_probs.shape)
      for j in range(seq_onehot.shape[1]):
        seq_onehot[x[j],j] = 1.
      self.seqs[i] = seq_onehot
    final_preds = self._test_sequences(self.seqs)
    print(final_preds)
    #print(final_preds.shape)

    #for i in range(self.seqs.shape[0]):
    #  print(debug_dists(self.seqs[0], self.seqs[i]))

    for i,p in enumerate(self.output_names):
      for j in range(final_preds.shape[-1]):
        this_pred = final_preds[:,i,j]
        this_name = p + '_' + str(j)
        ans[this_name] = this_pred
    return(pandas.DataFrame(ans))

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  random_seed = int(cfg.get('Params','random_seed'))
  random.seed(random_seed); np.random.seed(random_seed)
  num_iters = int(cfg.get('Params','NUM_ITERS'))
  evolver = seq_evolution_class_gradient(cfg)
  params = seq_evolution.unpack_params(cfg)
  ans = evolver.basic_iterative(num_iters, params)
  ans.to_csv(os.path.expanduser(cfg.get('Files','preds_fn')), index = False)
  scores_tracked = np.array(evolver.score_tracking)
  fn_scores = os.path.expanduser(cfg.get('Files','score_fn'))
  np.savetxt(fn_scores, scores_tracked, delimiter=',')
  seq_selection.main(cfg)
