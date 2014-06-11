#skrypt
#libraries -------
library(hsmm)
library(pROC)
library(pbapply)
#skrypt

#algorytm viterbiego
source("functions.R")
source("mydurationViterbi.R")

measure_region <- function(region, max.length = 32) {
  lengths <- table(region)
  res <- rep(0, max.length)
  lengths <- lengths[as.numeric(names(lengths))>0] #removing lengths smaller than 1
  
  start_l <- min(as.numeric(names(lengths)))
  end_l <- max(as.numeric(names(lengths)))
  if(prod(start_l:end_l %in% as.numeric(names(lengths)))){
    max_length <- length(lengths) #if all lengths are present in training set
  } else{
    max_length <- 1
    sl <- sum(lengths)
    while(sum(lengths[1:max_length])/sl <= 0.51) {
      max_length <- which.min(start_l:end_l %in% as.numeric(names(lengths))) - 1
      start_l <- lengths[max_length + 2]  #to assure that max_length wouldn't be too small
    }
  }
  max_length <- min(max_length, max.length)
  
  prop_lengths <- lengths[1:max_length]/sum(lengths[1:max_length])
  res[as.numeric(names(prop_lengths))] <- prop_lengths
  res
}


signal_hsmm_train <- function(train_data, test_data, aa_group, max_length = 32) {
  ts <- calc_t(train_data, aa_group)
  
  t1 <- ts[["t1"]]
  t2 <- ts[["t2"]]
  t3 <- ts[["t3"]]
  t4 <- ts[["t4"]]
  
  overall <- t4 #table(degenerate(unlist(analized_sequences), aa5))
  overall.probs <- overall/sum(overall)          
  overall.probs.log = log(overall.probs) #for viterbi
  
  lengths <- ts[["lengths"]]
  params <- apply(lengths, 2, measure_region, max.length = max.length)
  params <- cbind(params, rep(1/max.length, max.length))
  
  #setting params for hmm -------
  additional_margin = 10
  pipar <- c(1,0,0,0)
  tpmpar <- matrix(c(0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1,
                     0, 0, 0, 0), 4, byrow = TRUE)
  od <- matrix(c((t1/sum(t1))[1:4],
                 (t2/sum(t2))[1:4],
                 (t3/sum(t3))[1:4],
                 (t4/sum(t4))[1:4]), 4, byrow = TRUE)
  signal_hsmm(train_data, aa_group, pipar = pipar, tpmpar = tpmpar, od = od, 
                   overall.probs.log = overall.probs.log, params = params)
}

signal_hsmm <- function(list_prot, aa_group, pipar, tpmpar, 
                             od, overall.probs.log, params, max.length = 50) { 
  t(vapply(list_prot, function(prot) {
    probka <- as.numeric(degenerate(prot, aa_group)[1:max.length])
    probka <- na.omit(probka)
    viterbi.res <- duration.viterbi(probka, pipar, tpmpar, od, params)
    viterbi_path <- viterbi.res[["path"]]
    c.site <- ifelse(any(viterbi_path==3), max(which(viterbi_path==3)), length(probka))
    cuts <- c(cuts, c.site)
    prob.signal <- viterbi.res[["viterbi"]][c.site, viterbi_path[c.site]]
    prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka[1:c.site], 0)
    c(prob.signal = prob.signal, prob.non = prob.non, c.site = c.site)
  }, c(0, 0, 0)))
}

#building training set ----
train_size <- 800
test_size <- 100 #real number of test size is two times larger
train_pos <- read_uniprot("euk.txt", euk = TRUE)
train_neg <- read.fasta("euk_not.fasta", seqtype = "AA")

ind_pos <- sample(1:length(train_pos))
ind_neg <- sample(1:length(train_neg))

save(ind_pos, ind_neg, file = "session_dump.Rdata")

train_dat <- train_pos[ind_pos[1:train_size]] 
test_dat <- c(train_pos[ind_pos[(train_size + 1):(train_size + test_size)]],
              train_neg[ind_neg[1:test_size]])

ets <- c(rep(0, test_size), rep(1, test_size))
real_cs <- sapply(train_pos[ind_pos[(train_size + 1):(train_size + test_size)]], 
                  function(protein) attr(protein, "sig")[2])
#training
res <- signal_hsmm_train(train_dat, test_dat, aa5)
characteristics <- analyze_bihmm(res, ets)
cs <- data.frame(real.cs = real_cs, pred.cs = res[1:test_size, 3])

