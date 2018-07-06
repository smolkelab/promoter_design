mkdir -p ~/facs-seq_test/ZEV/models
nohup python ~/facs-seq/models/model_trainer.py train_model_ZEV.cfg > ~/facs-seq_test/ZEV/models/train.log &
