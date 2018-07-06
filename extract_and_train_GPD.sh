nohup ./GPD/miseq/assign_bins_miseq_GPD.sh >> ~/facs-seq_test/GPD/master.log &
nohup ./GPD/miseq/build_table_align_GPD.sh >> ~/facs-seq_test/GPD/master.log &
nohup ./GPD/nextseq/assign_bins_nextseq_GPD.sh >> ~/facs-seq_test/GPD/master.log &
nohup ./GPD/nextseq/build_table_exact_GPD.sh >> ~/facs-seq_test/GPD/master.log &
nohup ./GPD/nextseq/fit_means_GPD.sh >> ~/facs-seq_test/GPD/master.log &
nohup ./GPD/merge_miseq_nextseq_GPD.sh >> ~/facs-seq_test/GPD/master.log &
nohup ./GPD/filter_means_GPD.sh >> ~/facs-seq_test/GPD/master.log &
nohup ./GPD/models/train_model_GPD.sh >> ~/facs-seq_test/GPD/master.log &