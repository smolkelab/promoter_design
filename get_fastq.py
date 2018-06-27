import os
f_in = os.listdir('rawdata')
for f in f_in:
  f = os.path.join('rawdata',f)
  os.system('gunzip -k ' + f)
  os.system('mv ' + f[:-3] + ' fastq') # last 3 chars are '.gz'
