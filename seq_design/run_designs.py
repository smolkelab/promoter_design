import sys
import os
import pandas as pd

SCRIPTS_DIR = 'designs'
SCRIPTS_TABLE = 'build_design_cfgs_key.csv'

def main(this_instance, num_instances):
  scripts = [q for q in os.listdir(SCRIPTS_DIR) if q[-2:] == 'sh']
  tags = [int(q.split('_')[0]) for q in scripts]
  tag_to_script = dict(zip(tags, scripts))
  df = pd.read_csv(SCRIPTS_TABLE)
  for i in df['ID']:
    if i % num_instances == this_instance:
      os.system('./' + os.path.join(SCRIPTS_DIR, tag_to_script[i]))

if __name__ == '__main__':
  [this_instance, num_instances] = [int(q) for q in sys.argv[1:]]
  main(this_instance, num_instances)