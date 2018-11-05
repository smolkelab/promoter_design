# Given a FASTA file containing sequences for each oligo ordered and a means table,
# get the actual designed sequences, remove adapters,
# and add columns to the means table 'Experiment' for the ID of this subpool and 'Design' for the ID of this design within the subpool.
# Save the table with added columns as a new file.

import sys
import pandas as pd

def simple_gg(fwd, rev):
  fwd = fwd.split('GAGACC')[0][:-1] # keep the GG site in 'fwd'
  rev = rev.split('GGTCTC')[1][5:] # remove the GG site from 'rev'
  return fwd + rev

def strip_adapters(s, fwd_adapter, rev_adapter):
  s = s.split(fwd_adapter)[1]
  s = s.split(rev_adapter)[0]
  return s

def process_one_design(lines, fwd_adapter, rev_adapter):
  lines = [q.strip() for q in lines]
  [name, fwd, _, rev] = lines
  seq = strip_adapters(simple_gg(fwd, rev), fwd_adapter, rev_adapter)
  [experiment,_,design,_] = name.split('|')
  experiment = experiment.split('_')[1]
  design = int(design)
  return(seq, experiment, design)


def main(fn_designs, fn_means, fn_out, fwd_adapter, rev_adapter):
  design_dict = {}
  # build a dictionary linking designed sequences to the 'Experiment' and 'Design' values
  with open(fn_designs, 'r') as fd:
    lines = []; first = True
    for l in fd:
      lines.append(l)
      if len(lines) == 4:
        seq, experiment, design = process_one_design(lines, fwd_adapter, rev_adapter)
        design_dict[seq] = (experiment, design)
        lines = []

  # read in the means
  means = pd.read_csv(fn_means)
  
  def lookup_error_tolerant(dict, key):
    try:
      ans = dict[key]
    except KeyError:
      ans = ('None',-1)
    return ans
  
  dict_output = [lookup_error_tolerant(design_dict, q) for q in means['Seqs']]
  experiments = [q[0] for q in dict_output]
  designs = [q[1] for q in dict_output]
  means['Experiment'] = experiments
  means['Design'] = designs
  means.to_csv(fn_out, index = False)

if __name__ == '__main__':
  [fn_designs, fn_means, fn_out, fwd_adapter, rev_adapter] = [q.strip() for q in sys.argv[1:]]
  main(fn_designs, fn_means, fn_out, fwd_adapter, rev_adapter)
