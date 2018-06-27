# Input: CSV table with fields
# Seq, read_id_F, group_F, read_id_R, group_R
# Output: Python dictionary with key:read_id_F, value:reduced group ID
# Approach: group_F and group_R are nodes in a graph;
# there are edges between nodes that are mapped to the same read.
# (so this is a bipartite graph).
# Find the connected components of the graph; specifically, get a connected component ID
# for each read.
# Also, take a file of format (seq, bin, read_id) and create a file of format (seq, bin, final_group_id);
# sort this output file by final_group_id.

import sys
import os

from time import time
def timer(text, mark):
  print(text + ': ' + "{:.3g}".format(time() - mark))

class group_node(object):

  def __init__(self):
    self.reads = []
    self.groups = []
    self.touched = False

  def add_read(self, read, group):
    if read not in self.reads:
      self.reads.append(int(read))
    if group not in self.groups:
      self.groups.append(int(group))

# h/t Lukasz Kidzinski
# Function draws a progress bar in consol. It's especially useful for CPU verison which takes hours
# and an empty screen is frightening to people
def progress(count, total, status=''):
  pass
#    bar_len = 60
#    filled_len = int(round(bar_len * count / float(total)))

#    percents = round(100.0 * count / float(total), 1)
#    bar = '=' * filled_len + '-' * (bar_len - filled_len)

#    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
#    sys.stdout.flush()

# input file is CSV; each row has:
# read_id, fwd_group, rev_group
def process_input_file(filename_in):
  ans = []
  with open(filename_in, 'r') as fi:
    for line in fi:
      line_split = line.strip().split(',')
      read_id = line_split[0]
      ans.append([int(read_id), int(line_split[1]), int(line_split[2])])
  return(ans)

# Python won't let me do this recursively - stack depth limitations :(
'''
def recurse_over_lists(list_fwd, list_rev, curr_group, use_fwd, idx):
  if use_fwd:
    this_node = list_fwd[idx]
  else:
    this_node = list_rev[idx]

  if this_node.touched: # go back!
    return(set())
  this_node.touched = True
  curr_group = curr_group.union(this_node.reads)
  for group_id in this_node.groups:
    curr_group = curr_group.union( recurse_over_lists(list_fwd, list_rev, curr_group, not use_fwd, group_id) )
  
  return(curr_group)
'''

# instead of calling the same function recursively, try keeping a "to-do list" of sorts within one function call, using a while loop.
def dont_recurse_over_lists(list_all):
  curr_pos = 0
  groups_out = []
  
  while curr_pos != None:
    curr_gp_members_visited = []
    curr_gp_members_to_visit = [curr_pos]
    while len(curr_gp_members_to_visit) > 0:
      next_member = curr_gp_members_to_visit.pop()
      curr_gp_members_visited.append(next_member)
      list_all[next_member].touched = True
      for i in list_all[next_member].groups:
        if not list_all[i].touched:
          curr_gp_members_to_visit.append(i)
    curr_pos = update_list_pointers(curr_pos, list_all)
    groups_out.append(curr_gp_members_visited)
    if curr_pos != None:
      progress(curr_pos, len(list_all))
    else:
      progress(len(list_all), len(list_all))
  print('Out of group reduction')
  print('Number of groups: ' + str(len(groups_out)))
  reads_out = []
  for (i, group_list) in enumerate(groups_out):
    progress(i, len(groups_out))
    next_reads = []
    for gp in group_list:
      next_reads.extend(list_all[gp].reads)
    reads_out.append(set(next_reads))

  return(reads_out)

def update_list_pointers(idx, this_list):
  while(True):
    if idx == len(this_list) or idx == None:
      idx = None
      break
    if not this_list[idx].touched:
      break
    idx += 1
  return(idx)

def link_read_groups(raw_dat):
  num_fwd = max([int(q[1]) for q in raw_dat]) + 1
  num_rev = max([int(q[2]) for q in raw_dat]) + 1
  list_fwd = [group_node() for q in range(num_fwd)]
  list_rev = [group_node() for q in range(num_rev)]
  # for each read, add it to the appropriate entry in list_fwd and list_rev,
  # along with links to the other lists.
  # Then, join the two lists, for convenience.
  for q in raw_dat:
    [read_num, fwd_gp, rev_gp] = q
    list_fwd[int(fwd_gp)].add_read(read_num, rev_gp + len(list_fwd))
    list_rev[int(rev_gp)].add_read(read_num, fwd_gp)

  list_all = list_fwd
  list_all.extend(list_rev)

  # get rid of some stuff
  list_fwd = []
  list_rev = []
  raw_dat = []

  final_groups = dont_recurse_over_lists(list_all)
  return(final_groups)

def get_read_group_dict(groups):
  dict_out = {}
  for gp_id, gp in enumerate(groups):
    progress(gp_id, len(groups))
    for read_id in gp:
      if read_id in dict_out:
        raise Exception('read_id ' + str(read_id) + ' appeared more than once in final_groups')
      dict_out[read_id] = gp_id
  return(dict_out)

def main_method(file_idd_reads, filename_orig_groups, file_group_list, filename_out):
  t = time()
  lists = process_input_file(filename_orig_groups)
  final_groups = link_read_groups(lists)
   #write the final_groups, as an intermediate output
  with open(file_group_list, 'w') as fgl:
    for group in final_groups:
      line_out = ','.join([str(q) for q in group])
      fgl.write("%s\n" % line_out)
  timer('Link read groups', t)

  # convert the final groups into a dictionary for assigning groups to read IDs
  #t = time()
  group_dict = get_read_group_dict(final_groups)
  with open(file_idd_reads) as fi, open(filename_out, 'w') as fo:
    for line in fi:
      seq, bin_id, read = line.split(',')
      gp = group_dict[int(read.strip())]
      line_out = ','.join([str(seq), str(bin_id), str(gp)])
      fo.write("%s\n" % line_out)

  # Sort the final output, by group:
  os.system('sort -s -t , -k3,3 ' + filename_out + ' -o ' + filename_out)
  timer('Output and sort final group assignments',t)

#if __name__ == '__main__':
#  [filename_in, filename_out] = sys.argv[1:3]
#  main_method(filename_in, filename_out)
