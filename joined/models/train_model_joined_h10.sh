mkdir -p ~/facs-seq_test/joined/models
nohup python ~/facs-seq/models/model_trainer_h10.py train_model_joined_h10.cfg > ~/facs-seq_test/joined/models/train_h10.log &
