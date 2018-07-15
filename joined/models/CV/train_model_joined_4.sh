mkdir -p ~/facs-seq_test/joined/models
nohup python ~/facs-seq/models/model_trainer.py train_model_joined_CV.cfg 4 > ~/facs-seq_test/joined/models/train.log &
