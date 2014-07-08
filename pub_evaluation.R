library(seqinr)

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


#READ AND PREPROCESS DATA ------------------------------------
#select:(keyword:signal) AND reviewed:yes AND created:[1950 TO 2010] AND taxonomy:"Eukaryota [2759]" AND annotation:(type:signal confidence:experimental)
#3676 proteins
pos_train <- read_uniprot("pub_pos_train.txt", euk = TRUE)

#(keyword:signal) AND reviewed:yes AND created:[2011 TO 2014] AND taxonomy:"Eukaryota [2759]" AND annotation:(type:signal confidence:experimental)
#140 proteins
pos_test <- read_uniprot("pub_pos_test.txt", euk = TRUE)

#NOT (keyword:signal) AND reviewed:yes AND created:[2011 TO 2014] AND taxonomy:"Eukaryota [2759]"
#10634 seqsbut see below
neg_test <- read.fasta("pub_neg_test.fasta")

#eliminate atypical aa and too short seqs
atyp_aa <- which(sapply(neg_test, function(i) any(i %in% c("x", "j", "z", "b"))))
too_short <- which(sapply(neg_test, length) < 50)
neg_test <- neg_test[-unique(c(atyp_aa, too_short))]
#9795 seqs

chosen_neg <- sample(1:length(neg_test))[1:280]

#save test data
test_dat <- c(pos_test, neg_test[chosen_neg])
write.fasta(test_dat, names = names(test_dat), file.out = "pub_benchmark.fasta")

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
#http://www.cbs.dtu.dk/services/SignalP/
eval_signalp3nn <- read_signalp3_nn("eval_signalp3nn.txt")

#BENCHMARK - SignalP3.0 hidden markov models ------------------------------------
#http://www.cbs.dtu.dk/services/SignalP/
eval_signalp3hmm <- read_signalp3_hmm("eval_signalp3hmm.txt")

