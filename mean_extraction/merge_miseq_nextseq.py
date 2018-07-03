
# Merge sequences from MiSeq with mean estimates from NextSeq
# approach: read in NextSeq file as a dictionary; map sequence beginning to other values
# eliminate ambiguous sequences: all inputs in NextSeq will be unique,
# but multiples in MiSeq need to be dropped as they are found
import sys
import os
import ConfigParser

class LineJoiner(object):

  def __init__(self, short_seqs_fn, req_len):
    self.ss_dict = {}
    self.ss_set = set()
    self.req_len = req_len
    self.misses = 0
    self.drops = 0
    self.hits = 0
    self.needs_header = True
    with open(short_seqs_fn, 'r') as fi:
      for l in fi:
        if self.needs_header: # ignore the first line
          l = l.strip().split(',')
          self.seq_pos = [i for (i,q) in enumerate(l) if q == 'Seqs']
          assert(len(self.seq_pos) == 1)
          self.seq_pos = self.seq_pos[0]
          self.needs_header = False
        else:
          l = l.strip().split(',')
      	  seq = l.pop(self.seq_pos)
      	  assert(len(seq) == self.req_len)
      	  val_str = ','.join(l)
      	  # val_str is the 'payload'; 'None' corresponds to the matching long sequence.
      	  # If there are multiple long sequences, we can drop it entirely.
      	  self.ss_dict[seq] = [val_str, None]
      	  self.ss_set.add(seq)

  def merge_line(self, long_seq):
    s_seq = long_seq[:self.req_len]
    if s_seq not in self.ss_set:
      self.misses += 1
      return
    # drop the sequence if it's ambiguous
    if self.ss_dict[s_seq][1] != None:
      self.ss_set.remove(s_seq)
      del self.ss_dict[s_seq]
      self.drops += 1
      return
    # add the sequence
    self.ss_dict[s_seq][1] = long_seq
    self.hits += 1

  def dump_lines(self):
    short_seqs = self.ss_dict.keys()
    short_seqs.sort()
    lines_out = []
    for s in short_seqs:
      ans = self.ss_dict[s]
      if ans[1] != None:
      	ans = ','.join([ans[1],ans[0]]) + '\n'
      	lines_out.append(ans)
    return(lines_out)

def main(fn_miseq, fn_nextseq, fn_out, req_len):
  joiner = LineJoiner(fn_nextseq, req_len)
  with open(fn_miseq, 'r') as fm:
    for l in fm:
      joiner.merge_line(l.split(',')[0])
  lines_out = joiner.dump_lines()
  with open(fn_out, 'w') as fo:
    for l in lines_out:
      fo.write(l)
  return(joiner.misses, joiner.drops, joiner.hits)

if __name__ == '__main__':
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.read(sys.argv[1])

  fn_miseq = os.path.expanduser(config.get('Input','file_miseq'))
  fn_nextseq = os.path.expanduser(config.get('Input','file_nextseq'))
  fn_out = os.path.expanduser(config.get('Output','file_out'))
  req_len = int(config.get('Params','REQ_LEN'))
  misses, drops, hits = main(fn_miseq, fn_nextseq, fn_out, req_len)
  print('Misses: ' + str(misses))
  print('Drops: ' + str(drops))
  print('Hits: ' + str(hits))
