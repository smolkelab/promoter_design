import sys
import os
import pandas as pd

SCRIPTS_DIR = 'designs'
SCRIPTS_TABLE = 'build_design_cfgs_key.csv'

def main(this_instance, num_instances):
  df = pd.read_csv(SCRIPTS_TABLE)
  os.chdir(SCRIPTS_DIR)
  scripts = [q for q in os.listdir('.') if q[-2:] == 'sh']
  tags = [int(q.split('_')[0]) for q in scripts]
  tag_to_script = dict(zip(tags, scripts))
  for i in df['ID']:
    if i % num_instances == this_instance:
      os.system('./' + tag_to_script[i])

if __name__ == '__main__':
  [this_instance, num_instances] = [int(q) for q in sys.argv[1:]]
  main(this_instance, num_instances)