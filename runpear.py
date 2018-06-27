import os
INPUT = '/home/benkotopka/FS8/FASTQ'
OUTPUT = '/home/benkotopka/FS8/PEAR_out'
f_r1 = [q for q in os.listdir(INPUT) if '_R1_' in q]
f_r2 = [q for q in os.listdir(INPUT) if '_R2_' in q]
f_r1.sort()
f_r2.sort()
f_p = zip(f_r1, f_r2)
for f in f_p:
  out_name = f[0].split('L001')[0]
  os.system('pear -f ' + os.path.join(INPUT,f[0]) + ' -r ' + os.path.join(INPUT,f[1]) +
                ' -o ' + os.path.join(OUTPUT,out_name) + '  -j 8 -v 1 -g 2 -n 20')
