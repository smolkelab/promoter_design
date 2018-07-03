import sys
import os
import ConfigParser
config = ConfigParser.RawConfigParser(allow_no_value=True)
config.read(sys.argv[1])

dir_in = config.get('Dirs','dir_in')
dir_out = config.get('Dirs','dir_out')
f_r1 = [q for q in os.listdir(dir_in) if '_R1_' in q]
f_r2 = [q for q in os.listdir(dir_in) if '_R2_' in q]
f_r1.sort()
f_r2.sort()
f_p = zip(f_r1, f_r2)
for f in f_p:
  out_name = f[0].split('L001')[0]
  os.system('pear -f ' + os.path.join(dir_in,f[0]) + ' -r ' + os.path.join(dir_in,f[1]) +
                ' -o ' + os.path.join(dir_out,out_name) + '  -j 8 -v 1 -g 2 -n 20')
