mkdir -p ~/facs-seq_test/joined/models
nohup python ~/facs-seq/models/model_trainer_dropout.py train_model_joined_dropout.cfg > ~/facs-seq_test/joined/models/train_dropout.log &
