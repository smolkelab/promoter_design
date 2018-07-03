
# Given a "barcode filename" containing a CSV of prefixes, postfixes, and filename codes,
# and further given a prefix, a postfix, and a code,
# return either a single integer (the 0-delimited row in CSV matching the inputs)
# or two integers (for two ways of calculating, whether the first is 'official' and the second is for checking)
# Return None if a legal assignment cannot be made.

# Barcode file should have no header, and contain prefixes, postfixes, and filename codes in that order.
import os

class bin_assigner(object):

  def __init__(self, barcode_fn):

    def add_to_list(dict, k, v):
      if k in dict:
        dict[k].append(v)
      else:
        dict[k] = [v]
      return(dict)

    self.fix_dict = {}
    self.code_dict = {}
    self.prefix_only_dict = {}
    i = 0
    with open(barcode_fn, 'r') as fb:
      for l in fb:
        l = [q.strip() for q in l.split(',')]
        self.code_dict = add_to_list(self.code_dict, l[-1], i)
        fix = '|'.join([str(len(q)) for q in l[:-1]])
        self.fix_dict[fix] = i

        # collect the prefixes separately - even when the postfix is unavailable,
        # we can use the prefix to crosscheck
        prefix = fix.split('|')[0]
        self.prefix_only_dict = add_to_list(self.prefix_only_dict, prefix, i)
        i += 1

  def assign_read(self, prefix, postfix, file_code):
    file_codes = self._possible_fixes_codeonly(file_code)
    fix_codes = self._possible_fixes_prefixonly(prefix)
    if fix_codes == None or file_codes == None:
      return([None])
    inter_code = list(set(file_codes) & set(fix_codes))
    assert len(inter_code) < 2
    if len(inter_code) == 0:
      return([None])
    return(inter_code)

  def _possible_fixes_prefixonly(self, prefix):
    prefix = str(len(prefix))
    if prefix in self.prefix_only_dict:
      return(self.prefix_only_dict[prefix])
    return(None)

  def _possible_fixes_codeonly(self, code):
    if code in self.code_dict:
      return(self.code_dict[code])
    return(None)
