# Getting the configs and shell scripts right for sequence design is important.
# Given a CSV table containing the info for each experiment, build the config files for each experiment,
# as well as matching shell scripts.
# The shell scripts have numbers to help keep track of each experiment.

import sys
import os
import pandas
import ConfigParser

# which promoter is this?
PROMOTERS = {'GPD':['TACGTAAATAATTAATAGTAGTGACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTCTGGGTGNNNNNNNNNNNGGCATCCANNNNNNNNNNNNNNNNNNNNNNNNNGGCATCCANNNNNNNNATCCCAGCCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTATATAAAGMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMCACCAAGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHATGTCTAAAGGTGAAGAATTATTCAC'],
             'ZEV':['TTTATCATTATCAATACTCGCCATTTCAAAGAATACGTAAATAATTAATAGTAGTGACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCGTGGGCGNNNNNNNGCGTGGGCGNNNNNNNNNGCGTGGGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATAAGTATATAAAGACGGMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMCACCAAGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHATGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGG']} 
OBJECTIVES = {'Strength':['strong','np.mean','lambda x: K.mean(x, axis = 1)'],
              'Induced strength':['induced','seq_evolution.get_induced','lambda x: x[:,1]'],
              'AR':['AR','seq_evolution.merge_outputs_AR','lambda x: x[:,1] - x[:,0]']} # what are we optimizing for? (controls merge_outputs)
FILTERS = {True:['gcfilter','seq_evolution.gc_filter'],False:['nofilter','lambda x: 0']} # Are we applying the GC content filter? (controls seq_scores)
FUNCTIONS = {'Mean':['mean','np.mean'],'Mean-sd':['1sd','seq_evolution.mean_minus_sd']} # How are we combining the outputs from each model?
STRATEGIES = {'Screening':['screen','10000','seq_screening.py'],
              'Evolution to threshold':['evolve-thresh','100','seq_evolve_to_threshold.py'],
              'Evolution: cycle-limited':['evolve-cycle','100','seq_evolution.py'],
              'Gradient to threshold':['gradient-thresh','500','seq_gradient_evolve_to_threshold.py'],
              'Gradient: cycle-limited':['gradient-cycle','500','seq_gradient_evolution.py']} # What evolution/screening strategy are we using?

if __name__ == '__main__':
  exps_p = pandas.read_csv(sys.argv[1])
  exps = {}
  for p in exps_p:
    exps[p] = exps_p[p].tolist()
  assert(all([q in PROMOTERS for q in exps['Promoter']]))
  assert(all([q in OBJECTIVES for q in exps['Objective']]))
  assert(all([q in FILTERS for q in exps['Filter']]))
  assert(all([q in FUNCTIONS for q in exps['Function']]))
  assert(all([q in STRATEGIES for q in exps['Strategy']]))

  for i, (promoter, objective, filter, function, strategy, threshold) in enumerate(zip(exps['Promoter'], exps['Objective'], exps['Filter'], 
                                                                        exps['Function'], exps['Strategy'], exps['Threshold'])):
    fn_stem = '_'.join([str(i), promoter, OBJECTIVES[objective][0], FILTERS[filter][0], FUNCTIONS[function][0], STRATEGIES[strategy][0]])
    cfg = ConfigParser.RawConfigParser()
    # Force case sensitivity, cf. https://stackoverflow.com/questions/1611799/preserve-case-in-configparser
    cfg.optionxform=str

    cfg.add_section('Dirs')
    cfg.set('Dirs','weights_dir','~/facs-seq_test/joined/final_weights')
    
    cfg.add_section('Files')
    cfg.set('Files','model_fn','~/facs-seq/models/do_model.py')
    cfg.set('Files','preds_fn','~/facs-seq_test/seq_designs/preds/' + fn_stem + '.csv')
    cfg.set('Files','rejected_fn', '~/facs-seq_test/seq_designs/rejected/' + fn_stem + '.txt')
    cfg.set('Files','selected_fn', '~/facs-seq_test/seq_designs/selected/' + fn_stem + '_selected.txt')
    cfg.set('Files','score_fn', '~/facs-seq_test/seq_designs/scores/' + fn_stem + '.csv')
    
    cfg.add_section('Functions')
    cfg.set('Functions','merge_outputs_keras',OBJECTIVES[objective][2])
    cfg.set('Functions','merge_outputs',OBJECTIVES[objective][1])
    cfg.set('Functions','merge_models',FUNCTIONS[function][1])
    cfg.set('Functions','seq_scores',FILTERS[filter][1])
    cfg.set('Functions','choose_best_seqs','seq_evolution.greedy_choose_best_seqs')
    
    cfg.add_section('Params')
    cfg.set('Params','SEQ',PROMOTERS[promoter][0])
    cfg.set('Params','N','25:25:25:25')
    cfg.set('Params','M','28:09:09:54')
    cfg.set('Params','H','33:33:0:33')
    cfg.set('Params','INIT_NOISE','2e-1')
    cfg.set('Params','NUM_SEQS',STRATEGIES[strategy][1])
    cfg.set('Params','NUM_VARIANTS','20')
    cfg.set('Params','WTS_EXT','.h5')
    cfg.set('Params','OUTPUT_NAMES','A,B')
    cfg.set('Params','RANDOM_SEED','2017')
    cfg.set('Params','NUM_ITERS','250')
    cfg.set('Params','REJECT_MOTIFS','GGTCTC,GAGACC')
    cfg.set('Params','THRESH',threshold)
    cfg.set('Params','NUM_SEQS_FINAL','120')
    cfg.set('Params','PICK_TOP','True')
    cfg.set('Params','NUM_MUTATIONS','20:50,10:50,3:50,1:100')
    cfg.set('Params','KEEP_PARENT','False:200,True:50')
    cfg.set('Params','GRADIENT_STEP','5e-2:100,1e-2:100,2e-3:50')
    cfg.set('Params','NORMALIZE_POWER','1.:100,1.01:50,1.05:50,1.1:50')
    cfg.write(open(os.path.join('designs',fn_stem + '.cfg'), 'w'))
    script = 'nohup time python ../' + STRATEGIES[strategy][2] + ' ' + fn_stem + '.cfg > ~/facs-seq_test/seq_designs/logs/' + fn_stem + '.log &\n'
    script_fn = os.path.join('designs',fn_stem + '.sh')
    with open(script_fn, 'w') as sf:
      sf.write(script)
    os.system('chmod +x ' + script_fn)
