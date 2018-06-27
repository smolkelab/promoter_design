# patch for quick-and-dirty means extraction for prelim analysis of MiSeq data
# don't bother with the multiple 8 bins for now; drop 0_8-1 and 1_8-2 entirely
# also drop all other non-bin columns
# want 24 numerical columns total
# header from read table 'fullseqs_read_table.txt' in FS7:
#,1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,24,3,4,5,6,7,8,9,Seq

#AGGATAAGAT...AAA,100007*,0,0,0,0,0,2,0,0*,0,0,0,0,0,0,0,0,0,0,0,1,0,0*,0,0,0,0,0*,1*,3*,0*

def munge_line(l):
  l = l.strip().split(',')
  l = [l[0]] + l[2:9] + l[10:23] + l[24:28]
  il = [int(q) for q in l[1:]]
  if sum(il[:12]) < 10 or sum(il[12:]) < 10:
    return(None)
  return(','.join(l))

def main(fn_in, fn_out):
  with open(fn_in, 'r') as fi, open(fn_out, 'w') as fo:
    header = ',' + ','.join([str(x) for x in range(1,25)]) + ',Seq\n'
    fo.write(header)
    for l in fi:
      ml = munge_line(l)
      if ml != None:
        fo.write(ml + '\n')

if __name__ == '__main__':
  main('FS8_miseq_raw_read_table.txt', 'FS8_miseq_munged_reads.txt')
