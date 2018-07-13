mkdir -p ~/facs-seq_test/joined/models
nohup python ~/facs-seq/models/model_trainer_h20.py train_model_joined_h20.cfg > ~/facs-seq_test/joined/models/train_h20.log &
