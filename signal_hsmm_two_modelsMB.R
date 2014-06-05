#skrypt
#libraries -------
library(hsmm)
library(pROC)
library(pbapply)
#skrypt

#algorytm viterbiego
source("myViterbi.R")
source("functions.R")

bihmm <- function(list_prot, aa_group, pipar, tp, od, overall.probs.log, max.length = 50) {
  t(vapply(list_prot, function(prot) {
    probka <- as.numeric(degenerate(prot, aa_group)[1:max.length])
    probka <- na.omit(probka)
    viterbi.res <- viterbi(probka, pipar, tp, od)
    viterbi_path <- viterbi.res[["path"]]
    c.site <- ifelse(any(viterbi_path==3), max(which(viterbi_path==3)), length(probka))
    prob.signal <- viterbi.res[["viterbi"]][c.site, viterbi_path[c.site]]
    prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka[1:c.site], 0)
    c(prob.signal = prob.signal, prob.non = prob.non, c.site = c.site)
  }, c(0, 0, 0)))
}

test_bihmm <- function(train_data, test_data, aa_group) {
  ts <- calc_t(train_data, aa_group)
  
  t1 <- ts[["t1"]]
  t2 <- ts[["t2"]]
  t3 <- ts[["t3"]]
  t4 <- ts[["t4"]]
  
  overall <- t4 #table(degenerate(unlist(analized_sequences), aa5))
  overall.probs <- overall/sum(overall)          
  overall.probs.log = log(overall.probs) #for viterbi
  
  #setting params for hmm -------
  additional.aminoacids = 10 #aminoacids choosen after cleavage site
  pipar <- c(1,0,0,0)
  tp <- matrix(c(0.814, 0.186, 0, 0,
                 0, 0.91, 0.09, 0,
                 0, 0, 0.78, 0.22,
                 0, 0, 0, 1-1/additional.aminoacids), 4, byrow = TRUE)
  
  od <- matrix(c((t1/sum(t1))[1:4],
                 (t2/sum(t2))[1:4],
                 (t3/sum(t3))[1:4],
                 (t4/sum(t4))[1:4]), 4, byrow = TRUE)
  
  bihmm(test_data, aa_group, pipar, tp, od, overall.probs.log)
}


signal_hsmm_train <- function(train_data, aa_group) {
  ts <- calc_t(train_data, aa_group)
  
  t1 <- ts[["t1"]]
  t2 <- ts[["t2"]]
  t3 <- ts[["t3"]]
  t4 <- ts[["t4"]]
  
  overall <- t4 #table(degenerate(unlist(analized_sequences), aa5))
  overall.probs <- overall/sum(overall)          
  overall.probs.log = log(overall.probs) #for viterbi
  
  #setting params for hmm -------
  additional.aminoacids = 10 #aminoacids choosen after cleavage site
  pipar <- c(1,0,0,0)
  tp <- matrix(c(0.814, 0.186, 0, 0,
                 0, 0.91, 0.09, 0,
                 0, 0, 0.78, 0.22,
                 0, 0, 0, 1-1/additional.aminoacids), 4, byrow = TRUE)
  
  od <- matrix(c((t1/sum(t1))[1:4],
                 (t2/sum(t2))[1:4],
                 (t3/sum(t3))[1:4],
                 (t4/sum(t4))[1:4]), 4, byrow = TRUE)
  
  list(pipar = pipar, tp = tp, od = od, overall.probs.log = overall.probs.log)
}

signal_hsmm_test <- function(prot, aa_group, signal_hsmm_train_res, max.length = 50) {
  pipar <- signal_hsmm_train_res$pipar
  tp <- signal_hsmm_train_res$tp
  od <- signal_hsmm_train_res$od
  overall.probs.log <- signal_hsmm_train_res$overall.probs.log
  probka <- as.numeric(degenerate(prot, aa_group)[1:max.length])
  probka <- na.omit(probka)
  viterbi.res <- viterbi(probka, pipar, tp, od)
  viterbi_path <- viterbi.res[["path"]]
  c.site <- ifelse(any(viterbi_path==3), max(which(viterbi_path==3)), length(probka))
  prob.signal <- viterbi.res[["viterbi"]][c.site, viterbi_path[c.site]]
  prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka[1:c.site], 0)
  list(prob.signal = prob.signal, prob.non = prob.non, c.site = c.site, path = viterbi_path)  
}

#ets - vector of etiquettes (0 if signal)
#cs - vector of real cleave sites
analyze_bihmm <- function(bihmm_res, ets) {
  pred <- bihmm_res[, 1] - bihmm_res[, 2]
  pred_exp <- exp(pred - max(pred))
  list(auc = auc(response = ets, predictor = pred_exp),
       roc = roc(response = ets, predictor = pred_exp, plot=FALSE))
}

aa1 = list('1' = c("G", "A", "P", "V", "L", "I", "M"), 
           '2' = c("F", "W", "Y"), 
           '3' = c("K", "R", "H"), 
           '4' = c("D", "E"), 
           '5' = c("S", "T", "C", "N", "Q"))

aa2 = list('1' = c("G", "A", "P", "V", "L", "I", "M", "F"), 
           '2' =c("K", "R", "H"), 
           '3' =c("D", "E"), 
           '4' =c("S", "T", "C", "N", "Q", "Y", "W"))

aa3 = list('1' = c("G", "A", "P", "V", "L", "I", "M", "F", "W", "S", "T", "C", "N", "Q", "Y"), 
           '2' = c("D", "E"), 
           '3' = c("K", "R", "H"))

aa4 = list('1' = c("G", "A", "V", "L", "I"), 
           '2' = c("M", "C"), 
           '3' = c("F", "W", "Y"), 
           '4' = c("S", "T", "N", "Q"), 
           '5' = c("K", "R", "H"), 
           '6' = c("D", "E"), 
           '7' ="P")

aa5 = list('1' = c("K", "R", "H"),
           '2' = c("V","I","L","M","F","W","C"),
           '3' = c("S", "T", "N", "Q"),
           '4' = c("D","E","A","P","Y","G"))

aa6 <- list('1' = c("A", "C", "Q", "E"), 
            '2' = c("R", "H",  "K"), 
            '3' = c("N", "D", "G", "P", "S", "T"), 
            '4' = c("I", "L",  "M", "F", "W", "Y", "V"))

#building training set ----
train_size <- 800
test_size <- 100 #real number of test size is two times larger
train_pos <- read_uniprot("euk.txt", euk = TRUE)
train_neg <- read.fasta("euk_not.fasta", seqtype = "AA")

#one hundred repeats
hundred_reps <- pblapply(1:1000, function(unnecessary_argument) {
  ind_pos <- sample(1:length(train_pos))
  ind_neg <- sample(1:length(train_neg))
  
  train_dat <- train_pos[ind_pos[1:train_size]] 
  test_dat <- c(train_pos[ind_pos[(train_size + 1):(train_size + test_size)]],
                train_neg[ind_neg[1:test_size]])
  
  ets <- c(rep(0, test_size), rep(1, test_size))
  real_cs <- sapply(train_pos[ind_pos[(train_size + 1):(train_size + test_size)]], 
                    function(protein) attr(protein, "sig")[2])
  #training
  res <- test_bihmm(train_dat, test_dat, aa5)
  characteristics <- analyze_bihmm(res, ets)
  cs <- data.frame(real.cs = real_cs, pred.cs = res[1:test_size, 3])
  list(chars = characteristics, cs = cs)
})

all_cs <- do.call("rbind", lapply(hundred_reps, function(i) i[["cs"]]))
plot(density((all_cs[[1]] - all_cs[[2]])), main = "Density of cleave site misplacement")


plot(density(sapply(hundred_reps, function(i) i$chars$auc)))

#building final testing set ----
test_neg <- lapply(read.fasta("test_neg.fasta", seqtype = "AA"), toupper)
# test_neg <- test_neg[sapply(test_neg[1:100], length) > 100]
test_pos <- read_uniprot("test_pos.txt", euk = TRUE)
final_res <- test_bihmm(train_pos, append(test_pos, test_neg), aa5)

#special ets for signal.hsmm
ets <- c(rep(0, length(test_pos)), rep(1, length(test_neg)))
#normal ets
ets2 <- c(rep(1, length(test_pos)), rep(0, length(test_neg)))


signal_hsmm <- analyze_bihmm(final_res, ets)
signalp4 <- read_signalp4("signalp41_res.txt")
phobius <- read_phobius("phobius_res.txt")
predsi <- read.table("predsi_res.txt", sep = "\t", header = TRUE)
predsi[,4] <- as.numeric(predsi[,4]) - 1



save(ets, ets2, signal_hsmm, signalp4, phobius, predsi, hundred_reps, file = "prezent5_data.Rdata")
