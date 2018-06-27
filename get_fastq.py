import sys
import os
import ConfigParser
config = ConfigParser.RawConfigParser(allow_no_value=True)
config.read(sys.argv[1])

dir_in = config.get('Dirs','dir_in')
dir_out = config.get('Dirs','dir_out')
for f in os.listdir(dir_in):
  f = os.path.join(dir_in,f)
  os.system('gunzip -k ' + f)
  os.system('mv ' + f[:-3] + ' ' + dir_out) # last 3 chars are '.gz'