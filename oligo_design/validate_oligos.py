# TO DO: keep oligo_design.py from generating BsaI sites in the random pads!!!

# Given a table of oligos and a reference set of input sequences, do the following:
# Confirm that the design works:
# all oligo pairs generate their intended sequence
# all oligo pools have the same pool primers and different GG sites for each sequence
# all oligo pairs have the same constant primer, and PCR with constant primers generates intended product
# - just use exact matching as a proxy for amplification; hopefully the Tm calculations in primer_gen.py
# took care of any major issues there
# Then, design primers:
# constant primers for the entire table
# pool primers for each pool
# Name them appropriately based on extracted pool names.

import sys
import os
import ConfigParser
import re
import pandas as pd

DNA = ['A','C','G','T']
DNA_C = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def rc(seq):
  return(''.join([DNA_C[q] for q in seq])[::-1])

class simulate_gg(object):
  def __init__(self, fwd_primer_cfg, rev_primer_cfg, oligo_design_cfg):
    self.fwd_gg_regex = re.compile(oligo_design_cfg.get('Params','fwd_gg_cut').replace('N','.'))
    self.rev_gg_regex = re.compile(oligo_design_cfg.get('Params','rev_gg_cut').replace('N','.'))
    self.consts_f = oligo_design_cfg.get('Params','fwd_homol') + oligo_design_cfg.get('Params','fwd_toehold')
    self.consts_r = oligo_design_cfg.get('Params','rev_toehold') + oligo_design_cfg.get('Params','rev_homol')
    df_f = pd.read_csv(os.path.expanduser( fwd_primer_cfg.get('Files','out_fn') ))
    df_r = pd.read_csv(os.path.expanduser( rev_primer_cfg.get('Files','out_fn') ))
    self.pools_f = df_f['Full'].tolist()
    self.pools_r = df_r['Full'].tolist()
    self._toedict = {self.consts_f: oligo_design_cfg.get('Params','fwd_toehold'),
                     self.consts_r: oligo_design_cfg.get('Params','rev_toehold')}
    for (p,q) in zip(df_f['Toehold'], df_f['Full']):
      self._toedict[q] = p
    for (p,q) in zip(df_r['Toehold'], df_r['Full']):
      self._toedict[q] = p

  def _get_toe(self, full):
    return(self._toedict[full])

  def _digest_fwd(self, seq):
    seq_split = self.fwd_gg_regex.split(seq)
    assert(len(seq_split) == 2)
    seq_split = seq_split[0]
    seq = seq_split[:-4]
    gg_site = seq_split[-4:]
    return(seq, gg_site)

  def _digest_rev(self, seq):
    seq_split = self.rev_gg_regex.split(seq)
    assert(len(seq_split) == 2)
    seq_split = seq_split[1]
    seq = seq_split[4:]
    gg_site = seq_split[:4]
    return(seq, gg_site)

  def digest(self, seq, is_fwd):
    if is_fwd:
      return(self._digest_fwd(seq))
    else:
      return(self._digest_rev(seq))

  # Return the part of a sequence amplified by PCR primers (accounting for added bases),
  # or None if no exact matches to template for both sequences
  def simple_pcr(self, template, fwd, rev):
    f_toe = self._get_toe(fwd); r_toe = self._get_toe(rc(rev))
    print(f_toe)
    print(r_toe)
    print(template)
    if f_toe not in template or r_toe not in template:
      return(None)
    template = ''.join(template.split(f_toe)[1:])
    template = template.split(r_toe)[0]
    return( fwd + template + rev )

  def gg_two_piece(self, fwd, rev):
    fwd_seq, fwd_gg = self.digest(fwd, is_fwd = True)
    rev_seq, rev_gg = self.digest(rev, is_fwd = False)
    assert(fwd_gg == rev_gg)
    return(fwd_seq + fwd_gg + rev_seq, fwd_gg)

  '''def primer_search_pcr(self, template, is_fwd):
    if is_fwd: 
      pcrs = [(self.simple_pcr(template, self.consts_f, q), q) for q in self.pools_r]
      pcrs = [(p,q) for (p,q) in pcrs if p != None]
    else:
      pcrs = [(self.simple_pcr(template, q, self.consts_r), q) for q in self.pools_f]
    pcrs = [(p,q) for (p,q) in pcrs if p != None]
    print(pcrs)
    assert(len(pcrs) == 1)
    return(pcrs[0]) # product, primer tuple'''

  '''def simulate_pair(self, oligo_f, oligo_r):
    pcr_prod_f, pool_primer_f = self.primer_search_pcr(oligo_f, True)
    pcr_prod_r, pool_primer_r = self.primer_search_pcr(oligo_r, False)
    gg_prod, gg_site = self.gg_two_piece(pcr_prod_f, pcr_prod_r)
    return(pool_primer_f, pool_primer_r, gg_site, gg_prod)'''
	
  def simulate_line(self, fwd_oligo, rev_oligo, pool_fwd, pool_rev, const_fwd, const_rev):
    fwd_prod = self.simple_pcr(fwd_oligo, const_fwd, pool_rev)
	rev_prod = self.simple_pcr(rev_oligo, pool_fwd, const_rev)
	gg_prod = self.gg_two_piece(fwd_prod, rev_prod)
	assert(gg_prod == const_fwd + line_dict['Design'] + rc(const_rev))
	return(gg_prod)

class primers_from_oligos_by_name(object): # map oligos to primers by their name
  def __init__(self, oligo_design_cfg):
	primer_fn = oligo_design_cfg.get('Files','primer_fn')
	self.name_to_primer = {}
	lines = []
	with open(primer_fn, 'r') as fi:
	  for l in fi:
	    lines.append(l.strip())
		if len(lines) == 2:
		  [name, primer] = lines; lines = []
		  name = name.split('>')[1] # drop that '>'
		  if name == 'Const|F':
		    self.const_f = primer
		  elif name == 'Const_R':
		    self.const_r = primer
		  else:
		    self.name_to_primer[name] = primer
  def get_consts(self):
    return(self.const_f, self.const_r)
  def match_name(self, oligo_name):
    oligo_name = name.split('>')[1] # drop that '>'
	oligo_name = oligo_name.split('|')
	assert(len(oligo_name) == 4) # experiment, pool, sequence, F/RawConfigParser
	# drop the 'sequence' part
	oligo_name = oligo_name[:2] + oligo_name[3]
	oligo_name = '|'.join(oligo_name)
	return(self.name_to_primer(oligo_name))
		    
	  
'''def validate_df(df):
  pool_ids = df['pool_id'].unique()
  # oligos-wide tests:
  # "stem" name same for all sequences
  stem = pd.Series({'stem': [q.split('|')[0] for q in df['pool_id']]})
  print(stem)
  print(stem.unique())
  assert(len(stem.unique()) == 1)
  # all fwd constants same
  assert(len(df['const_primer_f'].unique()) == 1)
  # all rev constants same
  assert(len(df['const_primer_r'].unique()) == 1)
  # break df into pools
  pools = [df.loc[df['pool_id'] == q,:] for q in pool_ids]
  # within-pool tests:
  fwds = {'fwds':[]}
  revs = {'revs':[]}
  for pool in pools:
    # all fwd primers same
    pool_f = pool['pool_primer_f'].unique()
    assert(len(pool_f) == 1); fwds['fwds'].append(pool_f[0])
    # all rev primers same
    pool_r = pool['pool_primer_r'].unique()
    assert(len(pool_r) == 1); revs['revs'].append(pool_r[0])
    # all gg sites different
    assert(len(pool['gg_site'].unique()) == len(pool['gg_site']))
  # between-pool tests:
  # all fwd distinct
  fwds = pd.Series(fwds); assert(len(fwds.unique()) == len(fwds))
  # all rev distinct
  revs = pd.Series(revs); assert(len(revs.unique()) == len(revs))'''

# not exhaustive - see validate_df for some other things that should be checked
def simulate_oligo_file(fn_in, simulator, matcher):
  with open(fn_in 'r') as fi:
    lines = []
	for l in fi:
	  lines.append(l.strip())
	  if len(lines) == 4:
	    [fwd_name, fwd_oligo, rev_name, rev_oligo] = lines; lines = []
		pool_fwd = matcher.match_name(fwd_name)
		pool_fwd = matcher.match_name(rev_name)
		const_fwd, const_fwd = matcher.get_consts()
		gg_prod = simulator.simulate_line(fwd_oligo, rev_oligo, pool_fwd, pool_rev, const_fwd, const_rev)

def cfg_from_key(cfg_in, key):
  ans = ConfigParser.RawConfigParser(allow_no_value=True); ans.optionxform=str
  fn = os.path.expanduser(cfg_in.get('Configs',key))
  ans.read(fn)
  return(ans)

def get_simulator(cfg):
  fwd_primer_cfg, rev_primer_cfg, oligo_design_cfg = [cfg_from_key(cfg, q) for q in ['fwd_primer_cfg', 'rev_primer_cfg', 'oligo_design_cfg']]
  return( simulate_gg(fwd_primer_cfg, rev_primer_cfg, oligo_design_cfg) )
  
def main(cfg):
  simulator = get_simulator(cfg)
  oligo_design_cfg = cfg_from_key(cfg, 'oligo_design_cfg')
  matcher = primers_from_oligos_by_name(oligo_design_cfg)
  fn_in = os.path.expanduser(oligo_design_cfg.get('Files','oligo_fn'))
  fn_out = os.path.expanduser(cfg.get('Files','fn_out'))
  sim_df = simulate_oligo_file(fn_in, simulator, matcher)
  

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True); cfg.optionxform=str
  cfg.read(sys.argv[1])
  main(cfg)
