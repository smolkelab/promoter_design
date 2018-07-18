nohup time python ../seq_evolve_to_threshold.py 5_GPD_strong_nofilter_1sd_evolve-thresh.cfg > ~/facs-seq_test/seq_designs/logs/5_GPD_strong_nofilter_1sd_evolve-thresh.log &
nohup python ../seq_selection.py 5_GPD_strong_nofilter_1sd_evolve-thresh.cfg >> ~/facs-seq_test/seq_designs/logs/5_GPD_strong_nofilter_1sd_evolve-thresh.log &
