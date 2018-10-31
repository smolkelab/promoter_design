GPD_READ_TABLE_FN = ~/facs-seq_test/design_testing/miseq_GPD/final_merged_reads_GPD.csv
ZEV_READ_TABLE_FN = ~/facs-seq_test/design_testing/miseq_ZEV/final_merged_reads_ZEV.csv
OUT_FN = ~/facs-seq_test/design_testing/merged_read_table.csv
with open(GPD_READ_TABLE_FN, 'r') as fg, open(ZEV_READ_TABLE_FN, 'r') as fz, open(OUT_FN, 'w') as fo:
  for l in fg:
    fo.write(l)
  for l in fz:
    fo.write(l)
