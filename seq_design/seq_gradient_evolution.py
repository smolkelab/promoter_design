# Given an ensemble of models, evolve a random sequence to fulfill an objective.
# cf. https://blog.keras.io/how-convolutional-neural-networks-see-the-world.html

from keras import backend as K
from functools import partial

import seq_evolution

class seq_evolution_class_gradient(seq_evolution.seq_evolution_class):

  def __init__(self, cfg):
    # prepare the models, populate the sequences, et cetera
    super(seq_evolution_class_gradient, self).__init__()
    self.loss_tensor_fx = self._get_loss_tensor_fx(cfg) # given a model, get an output tensor
    iterates = [self._get_iterate_fx_from_model(q) for q in self.models]
    def _get_mean_grad(iterates, input):
      iterate_outputs = [q(input) for q in iterates]
      losses = [p for (p,q) in iterate_outputs]
      grads = [q for (p,q) in iterate_outputs]
      final_grad = K.mean(K.stack(grads, axis = 0), axis = 0)
      return(losses, final_grad)
    # given an input, get a list of losses, and the mean gradient of the input
    self.losses_and_grads = partial(_get_mean_grad, iterates)
    

  def _get_loss_tensor_fx(self, cfg):
    def loss_wrapper(merge_outputs_keras, model_in):
      #model_output = model_in.layers[-1].output
      model_output = model_in.output
      loss = merge_outputs_keras(model_output)
      return(loss)
    merge_outputs_keras = eval(cfg.get('Functions','merge_outputs')) # needs to be written in terms of the backend
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
    losses, gradient = self.losses_and_grads(seqs)
    print(losses)
    seqs += gradient*step
    return(seqs)
    
  # override method in parent class; update via gradient
  def iterate(self, params, iter_idx):
    self.seqs = self._update_seq_ensemble(self.seqs, float(params['gradient_step'][iter_idx]))

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