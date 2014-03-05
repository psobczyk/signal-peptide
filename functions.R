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


#remove seqs with atypical and/or not identified aas
get_atyp <- function(list_prots) {
  which(colSums(sapply(list_prots, function(prot) 
    sapply(c("B", "U", "X", "Z"), function(atyp_aa)
      atyp_aa %in% prot))) > 0)
}


read_uniprot <- function(connection, euk) {
  
  all_lines <- readLines(connection)
  
  all_seqs <- preliminary_seqs(all_lines, signal = TRUE) 
  
  only_sure <- remove_unsure(all_lines, all_seqs)
  if (euk) {
    only_nuclear <- remove_nonnuclear(all_lines, all_seqs)
    nuclear_and_sure <- intersect(only_nuclear, only_sure)
    sure_seqs <- all_seqs[, nuclear_and_sure]
  } else {
    sure_seqs <- all_seqs[, only_sure]
  }
  
  list_prots <- lapply(1L:ncol(sure_seqs), function(i) {
    start_seq <- sure_seqs[2,i]
    end_seq <- sure_seqs[3,i]
    
    ith_seq <- strsplit(gsub(" ", "", paste0(all_lines[start_seq:end_seq], collapse = "")), "")[[1]]
    
    class(ith_seq) <- "SeqFastaAA"
    attr(ith_seq, "name") <- strsplit(all_lines[sure_seqs[1,i]], "   ")[[1]][2]
    attr(ith_seq, "class") <- "SeqFastaAA"
    
    sig <- as.numeric(na.omit(as.numeric(strsplit(strsplit(all_lines[sure_seqs[4,i]], "SIGNAL       ")[[1]][2], " ")[[1]])))
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
  c(start_n = signal[1], start_h = start_h, start_c = start_c, cs = signal[2])
}


degenerate <- function(seq, aa_group) {
  for (i in 1L:length(aa_group)) {
    seq[seq %in% aa_group[[i]]] <- i
  }
  seq
}



calc_t <- function(list_prots) {
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
  
  t1 <- table(degenerate(n_region, aa5))
  t2 <- table(degenerate(h_region, aa5))
  t3 <- table(degenerate(c_region, aa5))
  t4 <- table(degenerate(rest, aa5))
  
  list(mean_cs = mean(nhc[, 4]), sd_cs = sd(nhc[, 4]), t1 = t1, t2 = t2, t3 = t3, t4 = t4)
}


#helper function starting model
start_model <- function(x, t1, t2, t3){
  #tu musz? zosta? wczytane dane
  numStates <- 3
  n_responses <- length(unique(x))
  #szanse wyprodukowania sygnalow
  state1 <- (t1/sum(t1))[1:n_responses] #n-region
  state2 <- (t2/sum(t2))[1:n_responses]  #h-region
  state3 <- (t3/sum(t3))[1:n_responses]  #c-region
  
  rModels <- list(
    list(
      GLMresponse(formula=x~1, family=multinomial(), pstart=state1, fixed=c(TRUE, TRUE, TRUE, TRUE, TRUE)[1:n_responses])
    ),
    list(
      GLMresponse(formula=x~1, family=multinomial(), pstart=state2, fixed=c(TRUE, TRUE, TRUE, TRUE, TRUE)[1:n_responses])
    ),
    list(
      GLMresponse(formula=x~1, family=multinomial(), pstart=state3, fixed=c(TRUE, TRUE, TRUE, TRUE, TRUE)[1:n_responses])
    )
  )
  
  #transition probs, wyliczone z dlugosci region?w
  transition <- list(transInit(~1, nstates=numStates, family=multinomial("identity"), pstart=c(0.82, 0.18, 0), fixed=c(TRUE, TRUE, TRUE)),
                     transInit(~1, nstates=numStates, family=multinomial("identity"), pstart=c(0, 0.92, 0.08), fixed=c(TRUE, TRUE, TRUE)),
                     transInit(~1, nstates=numStates, family=multinomial("identity"), pstart=c(0, 0, 1), fixed=c(TRUE, TRUE, TRUE)))
  
  #zaczynamy od konkretnego stanu
  instart = c(1, 0, 0)
  inMod <- transInit(~1, ns = 3, ps = instart, family = multinomial("identity"), data=data.frame(1), fixed=c(TRUE, TRUE, TRUE, TRUE))
  
  mod <- makeDepmix(response = rModels, transition = transition,prior=inMod)
  #summary(mod)
  #pars <- c(unlist(getpars(mod)))
  #length(pars)
  free <- c(0, 0, 0, rep(1, 9), rep(c(1, 1, 1), n_responses))
  fitted <- fit(mod, fixed=!free)
  #fitted@response
  #cbind(probka,posterior(fitted))
  #euk_nhc[numer_probki,]
  #summary(fitted)
  fitted
}

#helper function starting model
start_model2 <- function(x, t1, t2, t3, t4){
  #setting params -------
  additional_margin = 10
  pipar <- c(1,0,0,0)
  tpmpar <- matrix(c(0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1,
                     0.01, 0, 0, 0), 4, byrow = TRUE)
  #FIXME!!! --------------
  rdpar <- list(lambda = c(5.360593, 11.56915, 4.554336,additional_margin))
  #FIXME!!!!! ---------------
  odpar <- list(m = matrix(c((t1/sum(t1))[1:4],
                             (t2/sum(t2))[1:4],
                             (t3/sum(t3))[1:4],
                             (t4/sum(t4))[1:4]), 4, byrow = TRUE))
  
  viterbi_path <- hsmm.viterbi(x,od = "mult", rd = "pois",
                                 pi.par = pipar, tpm.par = tpmpar,
                                 od.par = odpar, rd.par = rdpar)
  return(viterbi_path)
}

#calculates concert of the single protein
test_protein <- function(protein, signal, aa_list, t1, t2, t3, t4, memory = FALSE) {
  nhc <- find_nhc(protein, signal)
  
  if (memory) {
    deg_protein <- make_transition_matrix(protein[1L:(nhc[4]+30)], aa_list)
  } else {
    deg_protein <- as.numeric(degenerate(protein[1L:(nhc[4]+30)], aa_list))
  }
  #fitted_model <- start_model(deg_protein, t1, t2, t3)
  #viterbi_path <- fitted_model@posterior[1L:nhc[4], 1]
  viterbi_path <- start_model2(deg_protein, t1, t2, t3, t4)$path[1L:nhc[4]]
  nhc[4] <- nhc[4] + 1
  expected <- c(rep(1, nhc[2] - 1), rep(2, nhc[3] - nhc[2]), rep(3, nhc[4] - nhc[3]))
  c(sum(viterbi_path==expected)/length(viterbi_path),
    vapply(1L:3, function(i) 
      sum(viterbi_path[expected == i] == i)/(nhc[i+1] - nhc[i]), 0))
  
}

#seeks signal peptides in the list of proteins
find_signal <- function(list_prots, ts, alpha = 2, paralell = FALSE) {
  
  #signal peptide must have at least 12 aa - assumption
  pot_sigs <- seq(from = ifelse(round(ts[[1]] - alpha*ts[[2]], 0) > 12, round(ts[[1]] - alpha*ts[[2]], 0), 12), 
                  to = round(ts[[1]] + alpha*ts[[2]], 0), by = 1)
  
  if (paralell) {
    sfLapply(list_prots, function(protein)
      recogn <- sapply(pot_sigs, function(cs) 
        try(test_protein(protein, c(1, cs), aa5, ts[["t1"]], ts[["t2"]], ts[["t3"]], ts[["t4"]]), silent = TRUE)
      )
    )
  } else {
    lapply(list_prots, function(protein)
      recogn <- sapply(pot_sigs, function(cs) 
        try(test_protein(protein, c(1, cs), aa5, ts[["t1"]], ts[["t2"]], ts[["t3"]], ts[["t4"]]), silent = TRUE)
      )
    )       
  }
}

#read result of signalp2; remeber to copy raw result to txt 
read_signalp2 <- function(connection) {
  all_lines <- readLines(connection)
  found_names <- grep(">", all_lines)
  all_names <- all_lines[found_names[seq(2L, length(found_names), by = 2)]]
  
  found_preds <- grep("Signal peptide probability", all_lines)
  res <- as.numeric(unlist(strsplit(all_lines[found_preds], ":"))[seq(2L, length(found_names), by = 2)])
  names(res) <- all_names
  res
}

#read result of signalp2; remeber to copy raw result to txt
read_signalp4 <- function(connection) {
  all_lines <- readLines(connection)
  found_names <- grep(">", all_lines)
  all_names <- all_lines[found_names]
  res <- sapply(strsplit(all_lines[found_names + 6], "   "), function(x) x[7])
  res <- as.numeric(as.factor(res)) - 1
  names(res) <- all_names
  res
}