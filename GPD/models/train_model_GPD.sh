mkdir -p ~/facs-seq_test/GPD/models
nohup python ~/facs-seq/models/model_trainer.py train_model_GPD.cfg > ~/facs-seq_test/GPD/models/train.log &
