#!/bin/bash
# My first script

theano-cache clear
theano-cache purge.
export PATH=$PATH:/usr/local/cuda/bin/:/usr/bin/python2.7
cp data_train/* extract_gmm_10cv/
cd extract_gmm_10cv/
perl gmm_all.pl
cd ../GCCA_toolbox/
matlab -r runallbnrc_sp1
cd ../extract_gmm_10cv/
matlab -r highml_sp
perl data_all2.pl
cd ../cnn_data/
perl cnn_all2.pl
cd ../cnn
python2 convolutional_mlp.py
perl cnn_features.pl
matlab -r format2
matlab -r mpqa_rnn
cd ../fuzzy
matlab -r load_test3;
matlab -r Classification_ANFIS2;


