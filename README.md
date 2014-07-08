signal-peptide
==============

* signal_hsmm_nonparametricMB.R - working algorithm
* pub_pos_train.txt - training data
* pub_pos_test.txt - test data (signal peptides)
* pub_neg_test.txt - test data (without signal peptides)


------------------------------------
program for exploring properties of signal peptides. 

functions.R - all important functions
procedures.R - current training and test procedure

Requirements:
* dhmm package https://github.com/sobek44/dhmm (load by library(hsmm))

Other stuff:
* aagui - package for aminoacid clusterization https://github.com/michbur/aagui