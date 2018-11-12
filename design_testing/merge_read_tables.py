import os
import pandas as pd

GPD_READ_TABLE_FN = os.path.expanduser('~/facs-seq_test/design_testing/miseq_GPD/final_merged_reads_GPD.csv')
ZEV_READ_TABLE_FN = os.path.expanduser('~/facs-seq_test/design_testing/miseq_ZEV/final_merged_reads_ZEV.csv')
OUT_FN = os.path.expanduser('~/facs-seq_test/design_testing/merged_read_table.csv')

pd_gpd = pd.read_csv(GPD_READ_TABLE_FN)
pd_zev = pd.read_csv(ZEV_READ_TABLE_FN)
pd_out = pd.concat([pd_gpd, pd_zev])
pd_out.to_csv(OUT_FN, index = False)