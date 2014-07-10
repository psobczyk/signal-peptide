library(seqinr)
library(XML)
library(hmeasure)
library(pROC)

#READ UNIPROT DATA -----------------------------
#helper function to get seqs from  .txt files
preliminary_seqs <- function(all_lines, signal) {
  prot_ids <- grep("\\<ID   ", all_lines)
  seqs_start <- grep("\\<SQ   ", all_lines) + 1
  seqs_end <-  grep("^//", all_lines) - 1
  #test if protein has information regarding signal
  
  prot_sig <- rep(NA, length(prot_ids))
  
  if (signal) {
    signals <- grep("FT   SIGNAL", all_lines)
    all_ids <- sort(c(prot_ids, signals), method = "quick")
    prot_sig[which(all_ids %in% signals) - 1L:length(signals)] <- signals
  }
  
  rbind(prot_ids, seqs_start, seqs_end, prot_sig)
}

#removes proteins with probable or potential signal peptides
#removes proteins without cleavage site for signalase
#return indices of proteins which surely have signal peptide
remove_unsure <- function(all_lines, all_seqs) {
  signals <- all_seqs[4, ]
  #remove with unknown signal peptide end
  (1L:ncol(all_seqs))[-unique(unlist(lapply(c(">", "<1", "?", "Or "), function(pattern) 
    grep(pattern, all_lines[signals], fixed = TRUE))))]
}

#proteins not encoded in nucleus
remove_nonnuclear <- function(all_lines, all_seqs) {
  nucl <- grep("OG   ", all_lines)
  all_ids <- sort(c(all_seqs[1, ], nucl), method = "quick")
  setdiff(1L:ncol(all_seqs), which(all_ids %in% nucl) - 1L:length(nucl))
}

#removes noncleavable seqs
remove_notcleaved <- function(all_lines, all_seqs) {
  not_cleaved <- grep("Not cleaved", all_lines)
  all_ids <- sort(c(all_seqs[1, ], not_cleaved), method = "quick")
  setdiff(1L:ncol(all_seqs), which(all_ids %in% not_cleaved) - 1L:length(not_cleaved))
}


#remove seqs with atypical and/or not identified aas
get_atyp <- function(list_prots) {
  which(colSums(sapply(list_prots, function(prot) 
    sapply(c("B", "U", "X", "Z"), function(atyp_aa)
      atyp_aa %in% prot))) > 0)
}


read_uniprot <- function(connection, euk) {
  
  all_lines <- readLines(connection)
  
  all_seqs <- preliminary_seqs(all_lines, signal = TRUE) 
  
  #remove unsure signal peptides
  only_sure <- remove_unsure(all_lines, all_seqs)
  #remove notcleaved signal peptides
  cleaved <- remove_notcleaved(all_lines, all_seqs)
  sure_cleaved <- intersect(only_sure, cleaved) 
  
  if (euk) {
    #remove proteins directed to nucleus
    only_nonnuclear <- remove_nonnuclear(all_lines, all_seqs)
    sure_cleaved <- intersect(only_nonnuclear, sure_cleaved)
  } 
  
  sure_seqs <- all_seqs[, sure_cleaved]
  
  
  list_prots <- lapply(1L:ncol(sure_seqs), function(i) {
    start_seq <- sure_seqs[2,i]
    end_seq <- sure_seqs[3,i]
    
    ith_seq <- strsplit(gsub(" ", "", paste0(all_lines[start_seq:end_seq], collapse = "")), "")[[1]]
    
    class(ith_seq) <- "SeqFastaAA"
    aa_name <- strsplit(all_lines[sure_seqs[1,i]], "   ")[[1]][2]
    attr(ith_seq, "name") <- aa_name
    attr(ith_seq, "Annot") <- paste0(">", aa_name)
    attr(ith_seq, "class") <- "SeqFastaAA"
    
    #to do - think about something smarter than suppressWarnings
    sig <- suppressWarnings(as.numeric(strsplit(strsplit(all_lines[sure_seqs[4,i]], "SIGNAL       ")[[1]][2], " ")[[1]]))
    sig <- as.numeric(na.omit(sig))
    attr(ith_seq, "sig") <- sig
    #line is preserved just to have an additional source of information
    attr(ith_seq, "line") <- all_lines[sure_seqs[4,i]]
    
    ith_seq
  })
  
  names(list_prots) <- sapply(list_prots, function(i) attr(i, "name"))
  atypical <- get_atyp(list_prots)
  if (length(atypical) > 0){
    list_prots[-atypical]
  } else {
    list_prots
  }
}

# READ OUTPUT FROM VARIOUS PREDICTORS --------------------------
read_predtat <- function(connection) {
  all_lines <- readLines(connection)
  #get decisions
  sig_pep <- grepl("Sec signal peptide predicted", all_lines)
  #get cleavage sites
  cleave_sites <- sapply(sapply(all_lines[sig_pep], function(i) 
    strsplit(i, "cleavage site: ", fixed = TRUE)[[1]][2], USE.NAMES = FALSE), function(j)
      strsplit(j, " ")[[1]][c(1, 3)], USE.NAMES = FALSE)
  sig_start <- rep(NA, length(all_lines))
  sig_end <- sig_start
  sig_start[sig_pep] <- as.numeric(cleave_sites[1, ])
  sig_end[sig_pep] <- as.numeric(cleave_sites[2, ])
  
  res <- data.frame(signal.peptide = sig_pep, 
                    sig.start = sig_start,
                    sig.end = sig_end) 
  rownames(res) <- sapply(strsplit(all_lines, " - "), function(i) i[1])
  res
}


read_signalp41 <- function(connection) {
  all_lines <- readLines(connection)
  do.call(rbind, lapply(all_lines, function(i) {
    line <- strsplit(i, " ")[[1]]
    line <- line[!line == ""]
    res <- data.frame(signal.peptide = line[10] == "Y",
                      sig.start = ifelse(line[10] == "Y", 1, NA),
                      sig.end = ifelse(line[10] == "Y", as.numeric(line[5]) - 1, NA))
    rownames(res) <- line[1]
    res
  }))
}

read_signalp3_nn <- function(connection) {
  all_lines <- readLines(connection)
  do.call(rbind, lapply(all_lines, function(i) {
    line <- strsplit(i, " ")[[1]]
    line <- line[!line == ""]
    res <- data.frame(signal.peptide = line[14] == "Y",
                      sig.start = ifelse(line[14] == "Y", 1, NA),
                      sig.end = ifelse(line[14] == "Y", as.numeric(line[6]) - 1, NA))
    rownames(res) <- line[1]
    res
  }))
}

read_signalp3_hmm <- function(connection) {
  all_lines <- readLines(connection)
  do.call(rbind, lapply(all_lines, function(i) {
    line <- strsplit(i, " ")[[1]]
    line <- line[!line == ""]
    res <- data.frame(signal.peptide = line[7] == "Y",
                      sig.start = ifelse(line[7] == "Y", 1, NA),
                      sig.end = ifelse(line[7] == "Y", as.numeric(line[4]) - 1, NA))
    rownames(res) <- line[1]
    res
  }))
}


read_predsi <- function(connection) {
  dat <- read.table(connection, sep = "\t")
  data.frame(signal.peptide = dat[, 4] == "Y",
             sig.start = ifelse(dat[, 4] == "Y", 1, NA),
             sig.end = ifelse(dat[, 4] == "Y", as.numeric(dat[, 3]), NA),
             row.names = dat[, 1])
}

read_phobius <- function(connection) {
  all_lines <- readLines(connection)
  all_lines <- all_lines[-1]
  splited <- strsplit(all_lines, " ")
  #remove "" characters
  purged <- t(sapply(splited, function(i) i[i != ""]))
  cl_sites <- sapply(which(purged[, 3] == "Y"), function(i)
    as.numeric(strsplit(strsplit(purged[i,4], "/")[[1]][1], "c")[[1]][[2]]))
  res <- data.frame(signal.peptide = purged[, 3] == "Y",
                    sig.start = ifelse(purged[, 3] == "Y", 1, NA),
                    sig.end = rep(NA, nrow(purged)), 
                    row.names = purged[, 1])
  res[purged[, 3] == "Y", "sig.end"] <- cl_sites
  res
}

read_philius <- function(connection) {
  require(XML)
  all_dat <- xmlToList(xmlTreeParse(connection, asTree = TRUE))
  seq_dat_id <- 1L:(length(all_dat)/2)*2
  #data for table
  table_dat <- sapply(seq_dat_id, function(i) 
    unlist(all_dat[i][[1]][[1]][c(24, 22)]))
  cleaved <- sapply(table_dat, function(i)
    !(is.null(i[1]) || is.na(i[1])))
  res <- data.frame(signal.peptide = cleaved,
                    sig.start = ifelse(cleaved, 1, NA),
                    sig.end = rep(NA, length(seq_dat_id)),
                    row.names = unlist(all_dat[1L:(length(all_dat)/2)*2 - 1]))
  res[cleaved, "sig.end"] <- as.numeric(sapply(table_dat[cleaved], function(i) i[2]))
  res
}

# HEURISTIC ALGORITHM

#function to find n, h and c regions in signal peptide
find_nhc <- function(protein, signal = NULL) {
  if (is.null(signal)) 
    signal <- attr(protein, "sig")
  
  sig <- protein[signal[1]:signal[2]]
  start_c <- length(sig) - 2
  
  #noh number of hydrophobic residues
  noh <- 0
  while(noh < 2 && start_c > 1) {
    start_c <- start_c - 1
    noh <- noh + ifelse(sig[start_c] %in% c("A", "I", "L", "M", "F", "W", "V"), 1, -1)
    noh <- ifelse(noh < 0, 0, noh)
  }
  start_c <- start_c + 2
  
  start_h <- start_c - 6
  #if statement to prevent negative values
  if (start_h > 1) {
    #nonh number of nonhydrophobic residues
    nonh <- 0
    #noc number of charged
    noc <- 0
    while(nonh < 3 && noc == 0 && start_h > 1) {
      start_h <- start_h - 1
      nonh <- nonh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), -1, 1)
      nonh <- ifelse(nonh < 0, 0, nonh)
      noc <- ifelse(sig[start_h] %in% c("R", "H", "K", "D", "E"), 1, 0)
    }
  } else {
    start_h <- 1
  }
  
  #prestart_c = start_c - 1 - prevents start_h > start_c
  prestart_c <- start_c - 1
  noh <- 0
  while(noh == 0 && start_h < prestart_c) {
    start_h <- start_h + 1
    noh <- noh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), 1, 0)
  }
  #c(start_n = signal[1], start_h = start_h, start_c = start_c, cs = signal[2])
  c(start_n = 1, start_h = start_h, start_c = start_c, cs = signal[2])
}


# SIGNAL-HSMM ------------------------------------
source("mydurationViterbi.R")

aa5 = list('1' = c("K", "R", "H"),
           '2' = c("V","I","L","M","F","W","C"),
           '3' = c("S", "T", "N", "Q"),
           '4' = c("D","E","A","P","Y","G"))

degenerate <- function(seq, aa_group) {
  for (i in 1L:length(aa_group)) {
    seq[seq %in% aa_group[[i]]] <- i
  }
  seq
}



calc_t <- function(list_prots, aa_list) {
  nhc <- t(vapply(list_prots, find_nhc, rep(0, 4)))
  
  n_region <- NULL
  h_region <- NULL
  c_region <- NULL
  rest <- NULL
  
  for(i in 1L:length(list_prots)){
    region_starts <- nhc[i, ]
    n_region <- c(n_region, list_prots[[i]][1:(region_starts[2] - 1)])
    h_region <- c(h_region, list_prots[[i]][region_starts[2]:(region_starts[3] - 1)])
    c_region <- c(c_region, list_prots[[i]][region_starts[3]:(region_starts[4] - 1)])
    rest <- c(rest, list_prots[[i]][region_starts[4]:length(list_prots[[i]])])
  }
  
  t1 <- rep(0, length(aa_list))
  temp <- table(degenerate(n_region, aa_list))
  t1[as.numeric(names(temp))] <- temp
  names(t1) <- 1:length(aa_list)
  
  t2 <- rep(0, length(aa_list))
  temp <- table(degenerate(h_region, aa_list))
  t2[as.numeric(names(temp))] <- temp
  names(t2) <- 1:length(aa_list)
  
  t3 <- rep(0, length(aa_list))
  temp <- table(degenerate(c_region, aa_list))
  t3[as.numeric(names(temp))] <- temp
  names(t3) <- 1:length(aa_list)
  
  t4 <- rep(0, length(aa_list))
  temp <- table(degenerate(rest, aa_list))
  t4[as.numeric(names(temp))] <- temp
  names(t4) <- 1:length(aa_list)
  
  len_c <- nhc[, "cs"] - nhc[, "start_c"]
  len_h <- nhc[, "start_c"] - nhc[, "start_h"]
  len_n <- nhc[, "start_h"] - nhc[, "start_n"]
  lengths <- matrix(c(len_n, len_h, len_c), ncol = 3)
  colnames(lengths) <- c("n", "h", "c")
  
  list(mean_cs = mean(nhc[, 4]), sd_cs = sd(nhc[, 4]), t1 = t1, t2 = t2, t3 = t3, t4 = t4, 
       lengths = lengths)
}

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
      start_l <- start_l + 1  #to assure that max_length wouldn't be too small
      max_length <- ifelse(max_length == 0, length(lengths), max_length)
    }
  }
  max_length <- min(max_length, max.length)
  
  prop_lengths <- lengths[1:max_length]/sum(lengths[1:max_length])
  res[as.numeric(names(prop_lengths))] <- prop_lengths
  res
}


signal_hsmm_train <- function(train_data, test_data, aa_group, max.length = 32) {
  train_data <- lapply(train_data, toupper)
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
  
  decisions <- signal_hsmm(test_data, aa_group, pipar = pipar, tpmpar = tpmpar, od = od, 
                           overall.probs.log = overall.probs.log, params = params)
  #change output to normal decision
  cbind(prob.sig = exp(decisions[,1] - decisions[,2]), 
        sig.end = decisions[, 3])
  #debug version of output
#   list(prob.sig = exp(decisions[,1] - decisions[,2]), 
#        sig.end = decisions[, 3], pipar = pipar, tpmpar = tpmpar, od = od, 
#        overall.probs.log = overall.probs.log, params = params, lengths = lengths)
}

signal_hsmm <- function(list_prot, aa_group, pipar, tpmpar, 
                        od, overall.probs.log, params, max.length = 50) { 
  t(vapply(list_prot, function(prot) {
    probka <- as.numeric(degenerate(toupper(prot), aa_group)[1:max.length])
    probka <- na.omit(probka)
    viterbi.res <- duration.viterbi(probka, pipar, tpmpar, od, params)
    viterbi_path <- viterbi.res[["path"]]
    c.site <- ifelse(any(viterbi_path==3), max(which(viterbi_path==3)), length(probka))
    prob.signal <- viterbi.res[["viterbi"]][c.site, viterbi_path[c.site]]
    prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka[1:c.site], 0)
    c(prob.signal = prob.signal, prob.non = prob.non, c.site = c.site)   
  }, c(0, 0, 0)))
}

#READ AND PREPROCESS DATA ------------------------------------
#select:(keyword:signal) AND reviewed:yes AND created:[1950 TO 2010] AND taxonomy:"Eukaryota [2759]" AND annotation:(type:signal confidence:experimental)
#3676 proteins
pos_train <- read_uniprot("pub_pos_train.txt", euk = TRUE)
#as above but year 2004
pos_train_hard <- read_uniprot("pos_hard_data.txt", euk = TRUE)
#as above but year 1997
pos_train_hardest <- read_uniprot("pos_hardest_data.txt", euk = TRUE)
#as above but year 1990
pos_train_ultrahard <- read_uniprot("pos_ultrahard_data.txt", euk = TRUE)

#code below commented - uncomment if there is a need for changing data set
# #(keyword:signal) AND reviewed:yes AND created:[2011 TO 2014] AND taxonomy:"Eukaryota [2759]" AND annotation:(type:signal confidence:experimental)
# #140 proteins
# pos_test <- read_uniprot("pub_pos_test.txt", euk = TRUE)
# 
# #NOT (keyword:signal) AND reviewed:yes AND created:[2011 TO 2014] AND taxonomy:"Eukaryota [2759]"
# #10634 seqsbut see below
# neg_test <- read.fasta("pub_neg_test.fasta")
# 
# #eliminate atypical aa and too short seqs
# atyp_aa <- which(sapply(neg_test, function(i) any(i %in% c("x", "j", "z", "b"))))
# too_short <- which(sapply(neg_test, length) < 50)
# neg_test <- neg_test[-unique(c(atyp_aa, too_short))]
# #9795 seqs
# 
# chosen_neg <- sample(1:length(neg_test))[1:280]
# 
# #save test data
# test_dat <- c(pos_test, neg_test[chosen_neg])
# write.fasta(test_dat, names = names(test_dat), file.out = "pub_benchmark.fasta")

#BENCHMARK - PRED-TAT ------------------------------------
#http://www.compgen.org/tools/PRED-TAT/submit
eval_predtat <- read_predtat("eval_predtat.txt")

#BENCHMARK - SignalP4.1 no tm ------------------------------------
#http://www.cbs.dtu.dk/services/SignalP/
eval_signalp41notm <- read_signalp41("eval_signalp41notm.txt")

#BENCHMARK - SignalP4.1 tm ------------------------------------
#http://www.cbs.dtu.dk/services/SignalP/
eval_signalp41tm <- read_signalp41("eval_signalp41tm.txt")

#BENCHMARK - SignalP3.0 neural nets ------------------------------------
#http://www.cbs.dtu.dk/services/SignalP-3.0/
eval_signalp3nn <- read_signalp3_nn("eval_signalp3nn.txt")

#BENCHMARK - SignalP3.0 hidden markov models ------------------------------------
#http://www.cbs.dtu.dk/services/SignalP-3.0/
eval_signalp3hmm <- read_signalp3_hmm("eval_signalp3hmm.txt")

#BENCHMARK - PredSi ------------------------------------
#http://www.predisi.de/
eval_predsi <- read_predsi("eval_predsi.txt")

#BENCHMARK - Phobius ------------------------------------
#http://phobius.sbc.su.se/
eval_phobius <- read_phobius("eval_phobius.txt")

#BENCHMARK - Philius ------------------------------------
#http://www.yeastrc.org/philius/pages/philius/runPhilius.jsp
eval_philius <- read_philius("eval_philius.xml")

#BENCHMARK - signal-hsmm ------------------------------------
eval_signalhsmm <- signal_hsmm_train(pos_train, read.fasta("pub_benchmark.fasta"), aa5)
eval_signalhsmm2 <- signal_hsmm_train(pos_train_hard, read.fasta("pub_benchmark.fasta"), aa5)
eval_signalhsmm3 <- signal_hsmm_train(pos_train_hardest, read.fasta("pub_benchmark.fasta"), aa5)
eval_signalhsmm4 <- signal_hsmm_train(pos_train_ultrahard, read.fasta("pub_benchmark.fasta"), aa5)

all_preds <- data.frame(c(rep(TRUE, 140), rep(FALSE, 280)),
                        eval_predtat[, "signal.peptide"],
                        eval_signalp41notm[, "signal.peptide"],
                        eval_signalp41tm[, "signal.peptide"],
                        eval_signalp3nn[, "signal.peptide"],
                        eval_signalp3hmm[, "signal.peptide"],
                        eval_predsi[, "signal.peptide"],
                        eval_phobius[, "signal.peptide"],
                        eval_philius[, "signal.peptide"],
                        eval_signalhsmm[, "prob.sig"],
                        eval_signalhsmm2[, "prob.sig"],
                        eval_signalhsmm3[, "prob.sig"],
                        eval_signalhsmm4[, "prob.sig"])
colnames(all_preds) <- c("real", 
                         "predtat", 
                         "signalp41notm",
                         "signalp41tm",
                         "signalp3nn",
                         "signalp3hmm",
                         "predsi",
                         "phobius",
                         "philius",
                         "signal-hsmm-2010",
                         "signal-hsmm-2004",
                         "signal-hsmm-1997",
                         "signal-hsmm-1990")
HMeasure(all_preds[, "real"], all_preds[, -1])[["metrics"]][, c("AUC", "H")]
auc(c(rep(TRUE, 140), rep(FALSE, 280)), eval_signalhsmm4[, "prob.sig"])


# debug_dat <- lapply(c(150, 200, 3000), function(i)
#   signal_hsmm_train(pos_train[1:i], read.fasta("pub_benchmark.fasta"), aa5))