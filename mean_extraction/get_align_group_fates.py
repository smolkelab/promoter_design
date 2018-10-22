# Determine how successful the read clustering process is.
# For each read group ID, determine how many reads were collected, and what its fate was.
# Output is CSV with columns for group ID, # reads in 'cluster_output' (before getting consensus sequences),
# reads in 'ambig', reads in 'N', reads in 'raw_read_table'.

import sys
import os
import ConfigParser

class fate_summary(object):
  def __init__(self, config):
    self.fn_cluster_output = os.path.expanduser(config.get('Files_Intermediate', 'cluster_output'))
    self.fn_ambig = os.path.expanduser(config.get('Output', 'file_ambig'))
    self.fn_N = os.path.expanduser(config.get('Output', 'file_N'))
    self.fn_aln = os.path.expanduser(config.get('Output', 'file_aligned'))
    self.fn_out = os.path.expanduser(config.get('Output', 'file_fates'))
    
    # each line in 'N' and 'ambig' is a read (last comma-separated field is read ID)
    # each line in 'aln' is a successful group.
    # create records of the format {group_id: {reads_aln:, reads_N:, reads_ambig:}}
    self.keyset = set()
    self.records = {}
    self.lines_to_reads(self.fn_cluster_output, 'reads_pre')
    self.lines_to_reads(self.fn_ambig, 'reads_ambig')
    self.lines_to_reads(self.fn_N, 'reads_N')
    self.lines_to_groups(self.fn_aln, 'reads_aln')
    
  def _new_record(self, idx):
    self.keyset.add(idx)
    self.records[idx] = {'reads_pre':0, 'reads_aln':0, 'reads_N':0, 'reads_ambig':0}
    
  def lines_to_reads(self, fn, key_id):
    with open(fn, 'r') as fi:
      for l in fi:
        gp_id = int(l.strip().split(',')[-1])
        if gp_id not in self.keyset:
          self._new_record(gp_id)
        self.records[gp_id][key_id] = self.records[gp_id][key_id] + 1
    
  def lines_to_groups(self, fn, key_id): # group ID is field 1, 2 on up are reads (sum these)
    with open(fn, 'r') as fi:
      for l in fi:
        l = l.strip().split(',')[1:]
        gp_id = int(l[0])
        num_reads = sum([int(q) for q in l[1:]])
        if gp_id not in self.keyset:
          self._new_record(gp_id)
        self.records[gp_id][key_id] = self.records[gp_id][key_id] + num_reads
    
  def write_to_file(self):
    with open(self.fn_out, 'w') as fo:
      fo.write('Group,Pre,Ambig,N,Aln\n')
      for k in self.keyset:
        rec = self.records[k]
        fo.write(','.join([str(q) for q in [k, rec['reads_pre'], rec['reads_ambig'], rec['reads_N'], rec['reads_aln']]]) + '\n')
      

def main_method(cfg):
  fates = fate_summary(cfg)
  fates.write_to_file()
