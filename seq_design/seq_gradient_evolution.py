# Given an ensemble of models, evolve a random sequence to fulfill an objective.
# cf. https://blog.keras.io/how-convolutional-neural-networks-see-the-world.html

from keras import backend as K
from functools import partial

import seq_evolution

class seq_evolution_class_gradient(seq_evolution.seq_evolution_class):

  def __init__(self, cfg):
    # prepare the models, populate the sequences, et cetera
    super(seq_evolution_class_gradient, self).__init__()
    self.loss = self._get_loss(cfg)

  def _get_loss(self, cfg):
    def loss_wrapper(merge_outputs_keras, model_in):
      layer_output = model.layers[-1].output
      loss = merge_outputs_keras(layer_output)
      loss = K.mean(loss)
      return(loss)
    merge_outputs_keras = eval(cfg.get('Functions','merge_outputs')) # needs to be written in terms of the backend
    loss_fx = partial(loss_wrapper, merge_outputs_keras)
    return(loss_fx)

  # get the normalized gradient for a sequence and one model
  def _gradient_seq_model(self, seq, model):
    grads = K.gradients(self.loss, seq)[0] # gradient at the input
    # following example, normalize gradient
    grads = grads/(K.sqrt(K.mean(K.square(grads))) + 1e-5)
    return(grads)

  # given an ensemble of models and a sequence, update with the mean of gradients  
  def _update_seq_ensemble(self, seq, step)
    gradients = np.stack([self._gradient_seq_model(seq, q) for q in self.models], axis = 0)
    gradients = np.apply_along_axis(np.mean, 0, gradients)
    seq += gradients*step
    return(seq)
    
  # override method in parent class; update via gradient
  def iterate(self, params, iter_idx):
    self.seqs = self._update_seq_ensemble(self.seqs, float(params['gradient_step'][iter_idx]))
    #new_seqs = []
    #for s in self.seqs:
    #  s_new = self._update_seq_ensemble(s, float(params['gradient_step'][iter_idx]))
    #  new_seqs.append(s_new)
    #self.seqs = np.stack(new_seqs, axis = 0)

  def basic_iterative(self, num_iters, params):
    # sequences were populated and models prepared by __init__
    for i in range(num_iters):
      self.iterate(params, i)
    return( self.generate_report() )

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
    for (i,s) in enumerate(self.seqs):
      x = np.argmax(s, axis = 1)
      ans = np.zeros(s.shape)
      ans[x] = 1.
      self.seqs[i] = ans
    final_preds = self._test_sequences(self.seqs)
    for i,p in enumerate(self.output_names):
      for j in range(final_preds.shape[-1]):
        ans[p + '_' + str(j)] = final_preds[:,i,j]
    return(pandas.DataFrame(ans))    

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True)
  cfg.read(sys.argv[1])
  random_seed = int(cfg.get('Params','random_seed'))
  random.seed(random_seed); np.random.seed(random_seed)
  num_iters = int(cfg.get('Params','NUM_ITERS'))
  evolver = seq_evolution_class_gradient(cfg)
  params = unpack_params(cfg)
  ans = evolver.basic_iterative(num_iters, params)
  ans.to_csv(os.path.expanduser(cfg.get('Files','preds_fn')), index = False)
  scores_tracked = np.array(evolver.score_tracking)
  fn_scores = os.path.expanduser(cfg.get('Files','score_fn'))
  np.savetxt(fn_scores, scores_tracked, delimiter=',')
  
  seq_selection.main(cfg)