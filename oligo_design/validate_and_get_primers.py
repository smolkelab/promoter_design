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
    for (p,q) in zip(df_f['Toehold'], self.pools_f):
      self._toedict[q] = p
    for (p,q) in zip(df_r['Toehold'], self.pools_r):
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
    f_toe = self._get_toe(fwd); r_toe = self._get_toe(rev)
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

  def primer_search_pcr(self, template, is_fwd):
    if is_fwd: 
      pcrs = [(self.simple_pcr(template, self.consts_f, q), q) for q in self.pools_r]
      pcrs = [(p,q) for (p,q) in pcrs if p != None]
    else:
      pcrs = [(self.simple_pcr(template, q, self.consts_r), q) for q in self.pools_f]
    pcrs = [(p,q) for (p,q) in pcrs if p != None]
    print(pcrs)
    assert(len(pcrs) == 1)
    return(pcrs[0]) # product, primer tuple

  def simulate_pair(self, oligo_f, oligo_r):
    pcr_prod_f, pool_primer_f = self.primer_search_pcr(oligo_f, True)
    pcr_prod_r, pool_primer_r = self.primer_search_pcr(oligo_r, False)
    gg_prod, gg_site = self.gg_two_piece(pcr_prod_f, pcr_prod_r)
    return(pool_primer_f, pool_primer_r, gg_site, gg_prod)

def extract_from_oligo_pair(lines):
  [name_f, oligo_f, name_r, oligo_r] = lines
  assert(name_f[-2:] == '|F'); assert(name_r[-2:] == '|R')
  name_f = name_f[:-2]; name_r = name_r[:-2]
  assert(name_f == name_r)
  name = name_f.split('|'); assert(len(name) == 3)
  pool_id = '|'.join(name[:2])
  return(pool_id, name[2], oligo_f, oligo_r) # pool ID, sequence ID within pool, oligos

def simulate_oligo_file(fn_in, simulator):
  lines = []
  pool_ids = []
  seq_ids = []
  oligos_f = []
  oligos_r = []
  pool_primers_f = []
  pool_primers_r = []
  gg_sites = []
  gg_prods = []
  print(fn_in)
  with open(fn_in, 'r') as fi:
    for l in fi:
      l = l.strip()
      print(l)
      lines.append(l)
      if len(lines) == 4:
        pool_id, seq_id, oligo_f, oligo_r = extract_from_oligo_pair(lines); lines = []
        pool_primer_f, pool_primer_r, gg_site, gg_prod = simulator.simulate_pair(oligo_f, oligo_r)
        pool_ids.append(pool_id)
        seq_ids.append(seq_id)
        oligos_f.append(oligo_f)
        oligos_r.append(oligo_r)
        pool_primers_f.append(pool_primer_f)
        pool_primers_r.append(pool_primer_r)
        gg_sites.append(gg_site)
        gg_prods.append(gg_prod)
  consts_f = [simulator.consts_f]*len(pool_ids)
  consts_r = [simulator.consts_r]*len(pool_ids)
  ans = {'pool_id': pool_ids, 'seq_id': seq_ids, 'oligo_f': oligos_f, 'oligo_r': oligos_r, 'pool_primer_f': pool_primers_f,
           'pool_primer_r': pool_primers_r, 'gg_site': gg_sites, 'gg_prod': gg_prods, 'const_primer_f': consts_f, 'const_primer_r': consts_r}
  return(pd.DataFrame(ans))

def validate_df(df):
  print(df)
  pool_ids = df['pool_id'].unique()
  print(pool_ids)
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
  revs = pd.Series(revs); assert(len(revs.unique()) == len(revs))
  
def validate_and_generate_oligo_fasta(df, fn_out):
  validate_df(df)
  stem = df['pool_ids'][0].split('|')[0]
  # constant oligos
  const_f = df['const_primer_f'][0]
  const_r = df['const_primer_r'][0]
  # aggregate pools
  pool_ids = df['pool_id'].unique()
  pool_ids.sort()
  pools_collected = {}
  for q in pool_ids:
    primer_f = df.loc[df['pool_id'] == q,'pool_primer_f'][0]
    primer_r = df.loc[df['pool_id'] == q,'pool_primer_r'][0]
    pools_collected[q] = {'primer_f':primer_f, 'primer_r':primer_r}
  with open(fn_out, 'w') as fo:
    # constant oligos
    fo.write('>' + stem + '|const_f' + '\n'); fo.write(const_f + '\n')
    fo.write('>' + stem + '|const_r' + '\n'); fo.write(rc(const_r) + '\n')
    # pools
    for q in pool_ids:
      fo.write('>' + q + '|pool_f\n'); fo.write(pools_collected[q]['primer_f'] +'\n' )
      fo.write('>' + q + '|pool_r\n'); fo.write(rc(pools_collected[q]['primer_r']) +'\n' )

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
  fn_in = os.path.expanduser(oligo_design_cfg.get('Files','oligo_fn'))
  fn_out = os.path.expanduser(cfg.get('Files','fn_out'))
  sim_df = simulate_oligo_file(fn_in, simulator)
  validate_and_generate_oligo_fasta(sim_df, fn_out)

def debug_sim():
  cfg = ConfigParser.RawConfigParser(allow_no_value=True); cfg.optionxform=str
  cfg.read('validate_and_get_primers.cfg')
  return(get_simulator(cfg))

if __name__ == '__main__':
  cfg = ConfigParser.RawConfigParser(allow_no_value=True); cfg.optionxform=str
  cfg.read(sys.argv[1])
  main(cfg)
