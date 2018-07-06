mkdir -p ~/facs-seq_test/ZEV
nohup ./ZEV/miseq/assign_bins_miseq_ZEV.sh >> ~/facs-seq_test/ZEV/master.log &
nohup ./ZEV/miseq/build_table_align_ZEV.sh >> ~/facs-seq_test/ZEV/master.log &
nohup ./ZEV/nextseq/assign_bins_nextseq_ZEV.sh >> ~/facs-seq_test/ZEV/master.log &
nohup ./ZEV/nextseq/build_table_exact_ZEV.sh >> ~/facs-seq_test/ZEV/master.log &
nohup ./ZEV/nextseq/fit_means_ZEV.sh >> ~/facs-seq_test/ZEV/master.log &
nohup ./ZEV/merge_miseq_nextseq_ZEV.sh >> ~/facs-seq_test/ZEV/master.log &
nohup ./ZEV/filter_means_ZEV.sh >> ~/facs-seq_test/ZEV/master.log &
nohup ./ZEV/models/train_model_ZEV.sh >> ~/facs-seq_test/ZEV/master.log &
