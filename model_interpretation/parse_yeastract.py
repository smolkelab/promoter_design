# Output from the YEASTRACT database looks like this:
#Target Sequence: 48 (size 247)
#Back to top top
#Abf1p TNNCGTNNNNNNTGAT  -193  F
#Ash1p YTGAT -182  F
#Ash1p YTGAT -70 F
#Ash1p YTGAT -154  R
#Gcn4p TGACTMT -233  R
#Gcr1p CTTCC -103  F
#Gcr1p CTTCC -52 F
#Gcr1p CTTCC -14 F
#Gcr1p CWTCC -103  F
#Gcr1p CWTCC -52 F
#Gcr1p CWTCC -14 F
#Mot3p AAGAGG  -74 R
#Mot3p AAGAGG  -94 R
#Mot3p TMGGAA  -7  R
#Nrg1p CCCTC -81 F
#Rtg1p, Rtg3p  GTCAC -140  R
#Stb5p CGGNS -220  F
#Stb5p CGGNS -135  F
#Stb5p CGGNS -129  R
#Stb5p CGGNS -133  R
#Haa1p SMGGSG  -132  R
#
# Create a table like this:
#Seq_ID Abf1p Ash1p [etc...]
#48     1     3
# etc. - row for each sequence analyzed

import sys
import os

def get_raw_results(fn, new_result_key = 'Target Sequence: ', skip_keys = ['Back to top top']):
  with open(fn, 'r') as f:
    lines = f.readlines() 
  ans = []; curr_res = None
  for l in lines:
    if new_result_key in l:
      seq_name = l.split(new_result_key)[1].split(' ')[0]
      if curr_res is not None:
        ans.append(curr_res)
      curr_res = [seq_name,[]]
    else:
      if not any([q in l for q in skip_keys]):
        curr_res[1].append(l)
  return ans

def count_tfs_one_result(res):
  [seq_name, lines] = res
  ans = {}
  for l in lines:
    l = l.strip().split('\t')[0]
    try:
      ans[l] = ans[l] + 1
    except KeyError:
      ans[l] = 1
  return [seq_name, ans]

def get_all_keys(dict_list):
  ans = set()
  for d in dict_list:
    ks = set(d.keys())
    ans = ans.union(ks)
  ans = list(ans)
  ans.sort()
  return ans

def pad_result(res, all_keys):
  [seq_name, seq_dict] = res
  for k in all_keys:
    if k not in seq_dict.keys():
      seq_dict[k] = 0
  return [seq_name, seq_dict]

def write_file(fn_out, res_list, all_keys, rownames_name = 'Seq_ID'):
  with open(fn_out, 'w') as fo:
    x = [rownames_name]
    x.extend(all_keys)
    header = '\t'.join(x)
    fo.write(header + '\n')
    for r in res_list:
      [seq_name, ans_dict] = r
      cts = [str(ans_dict[q]) for q in all_keys]
      x = [seq_name]
      x.extend(cts)
      fo.write('\t'.join(x) + '\n')

def main(fn_in, fn_out):
  result_list = get_raw_results(fn_in)
  result_list = [count_tfs_one_result(q) for q in result_list]
  all_keys = get_all_keys([q[1] for q in result_list])
  result_list = [pad_result(q, all_keys) for q in result_list]
  write_file(fn_out, result_list, all_keys)

if __name__ == '__main__':
  [fn_in, fn_out] = [os.path.expanduser(q) for q in sys.argv[1:]]
  main(fn_in, fn_out)
