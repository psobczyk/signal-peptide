#"Error in if (any(y > 1)) stop(\"multinomial response with n>1 not allowed with link='mlogit'\") : \n  missing value where TRUE/FALSE needed\n"


#IMPORTANT
#Nielsen used Swiss-Prot 29 (June 1994) which contains 38303 proteins total
#created:[1950 TO 1995] AND reviewed:yes contains 38440 proteins total, so let's say that we got nearly the same data sets

#Nielsen got 2282 entries eukaryotic entries with signal peptide (after data purification)
#following select grabs 2382 entries, after purification (removing entries with <1, ? and alternate cleavage sites) there are 2299 entries
#select: (keyword:signal) AND reviewed:yes AND created:[1950 TO 1995] AND taxonomy:"Eukaryota [2759]" AND annotation:(type:signal confidence:experimental)

#negative set
#NOT (keyword:signal) AND reviewed:yes AND created:[1950 TO 1995] AND taxonomy:"Eukaryota [2759]"


#Nielsen got 579 entries prokaryotic entries with signal peptide (after data purification)
#following select grabs 603 entries, after purification (removing entries with <1, ? and alternate cleavage sites) there are 597 entries
#(select: (keyword:signal) AND reviewed:yes AND created:[1950 TO 1995] AND annotation:(type:signal confidence:experimental) NOT keyword:"Lipoprotein [KW-0449]" AND (taxonomy:"Bacteria [2]" OR taxonomy:"Archaea [2157]"))


#negative set
#(NOT (keyword:signal) AND reviewed:yes AND created:[1950 TO 1995] AND (taxonomy:"Bacteria [2]" OR taxonomy:"Archaea [2157]"))



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

#old version
# remove_unsure <- function(all_lines, all_seqs) {
#   signals <- all_seqs[4, ]
#   not_probable <- signals[-grep("Probable", all_lines[signals])]
#   not_potential <- not_probable[-grep("Potential", all_lines[not_probable])]
#   cleaved <- not_potential[-grep("Not cleaved.", all_lines[not_potential])]
#   #remove with unknown signal peptide end
#   sure_cleaved <- cleaved[-grep("?", all_lines[cleaved], fixed = TRUE)]
#   #remove with more than 1 signal peptide end
#   sure_cleaved2 <- sure_cleaved[-grep("Or ", all_lines[sure_cleaved], fixed = TRUE)]
#   #remove with signal peptide start <1
#   sure_cleaved3 <- sure_cleaved2[-grep("<1", all_lines[sure_cleaved2], fixed = TRUE)]
#   all_ids <- sort(c(all_seqs[1, ], sure_cleaved3), method = "quick")
#   which(all_ids %in% sure_cleaved3) - 1L:length(sure_cleaved3)
# }

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
  
  list_prots[-get_atyp(list_prots)]
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



start_model <- function(x, t1, t2, t3){
  #tu musz¹ zostaæ wczytane dane
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
  
  #transition probs, wyliczone z dlugosci regionów
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

test_protein <- function(protein, signal, aa_list, t1, t2, t3, memory = FALSE) {
  nhc <- find_nhc(protein, signal)
  
  if (memory) {
    deg_protein <- make_transition_matrix(protein[1L:(nhc[4]+30)], aa_list)
  } else {
    deg_protein <- as.numeric(degenerate(protein[1L:(nhc[4]+30)], aa_list))
  }
  fitted_model <- start_model(deg_protein, t1, t2, t3)
  viterbi_path <- fitted_model@posterior[1L:nhc[4], 1]
  nhc[4] <- nhc[4] + 1
  expected <- c(rep(1, nhc[2] - 1), rep(2, nhc[3] - nhc[2]), rep(3, nhc[4] - nhc[3]))
  c(sum(viterbi_path==expected)/length(viterbi_path),
    vapply(1L:3, function(i) 
      sum(viterbi_path[expected == i] == i)/(nhc[i+1] - nhc[i]), 0))
  
}


find_signal <- function(list_prots, ts, alpha = 2, paralell = FALSE) {
  
  #signal peptide must have at least 12 aa - assumption
  pot_sigs <- seq(from = ifelse(round(ts[[1]] - alpha*ts[[2]], 0) > 12, round(ts[[1]] - alpha*ts[[2]], 0), 12), 
                  to = round(ts[[1]] + alpha*ts[[2]], 0), by = 1)
  
  if (paralell) {
    sfLapply(list_prots, function(protein)
      recogn <- sapply(pot_sigs, function(cs) 
        try(test_protein(protein, c(1, cs), aa5, ts[["t1"]], ts[["t2"]], ts[["t3"]]), silent = TRUE)
      )
    )
  } else {
    lapply(list_prots, function(protein)
      recogn <- sapply(pot_sigs, function(cs) 
        try(test_protein(protein, c(1, cs), aa5, ts[["t1"]], ts[["t2"]], ts[["t3"]]), silent = TRUE)
      )
    )       
  }
}


encoding <- function(x){
  a <- x%/%10
  b <- x%%10
  b+4*(a-1)
}

make_transition_matrix <- function(region, aa_groups){
  b <- degenerate(region, aa_groups)
  tabela <- table(apply(cbind(b[-length(b)],b[-1]), 1, FUN=function(x) encoding(as.numeric(paste(x[1], x[2], sep="")))))
  for(i in 0:3){
    suma <- sum(tabela[(1+i*4):((i+1)*4)])
    for(j in 1:4){
      tabela[i*4+j] <- tabela[i*4+j]/suma
    }
  }
  return(tabela)
}

degenerate_mem <- function(protein, aa_list) {
  deg_aa <- degenerate(protein, aa_list)
  aa_pairs <- apply(expand.grid(1L:length(aa_list), 1:length(aa_list)), 1, function(row) 
    paste0(row, collapse = ""))
  len <- length(deg_aa)
  
  id1_start <- seq(1, len, by = 2)
  id1_end <- seq(2, len, by = 2)
  id2_start <- id1_end
  id2_end <- id1_start[-1]
  
  if (len %% 2 == 1) {
    id1_start <- id1_start[-length(id1_start)]
  } else {
    id2_start <- id2_start[-length(id2_start)]
  }
  
  as.vector(table(factor(c(apply(cbind(id1_start, id1_end), 1, function(row) 
    paste0(deg_aa[row], collapse = "")),
    apply(cbind(id1_start, id1_end), 1, function(row) 
      paste0(deg_aa[row], collapse = ""))), levels = aa_pairs)))
}

transition1 <- make_transition_matrix(n_region, aa5)
transition2 <- make_transition_matrix(h_region, aa5)
transition3 <- make_transition_matrix(c_region, aa5)
transition4 <- make_transition_matrix(rest, aa5)

# values ----------------------------

aa5 = list(positively.charged=c("K", "R", "H"),
           hydrofobic=c("V","I","L","M","F","W","C"),
           polar.uncharged=c("S", "T", "N", "Q"),
           rest=c("D","E","A","P","Y","G"))

# libraries ----------------------------

library(depmixS4)
library(seqinr)
library(snowfall)
library(kernlab)
library(ROCR)


#code

euk <- read_uniprot("euk.txt", euk = TRUE)
euk_not <- read.fasta("euk_not.fasta", seqtype = "AA")
euk_not <- euk_not[-get_atyp(euk_not)]


prok <- read_uniprot("prok.txt", euk = FALSE)

split_prots <- function(list_prot, train = 0.5, valid = 0.3, replace = TRUE) {
  len <- length(list_prot)
  ids <- sample(len, replace = replace)
  train_len <- round(len*train, 0)
  valid_len <- train_len + 1 + round(len*valid, 0)
  train_id <- ids[1L:train_len]
  valid_id <- ids[(train_len + 1):valid_len]
  test_id <- ids[(valid_len + 1):len]
  list(train = ids[1L:train_len],
       valid = ids[(train_len + 1):valid_len],
       test = ids[(valid_len + 1):len])
}

list_to_df <- function(results, target) {
  results <- t(do.call(cbind, lapply(results[sapply(results, class) == "matrix"], rowMeans)))
  cbind(data.frame(results), tar = rep(target, nrow(results)))
}



pos <- split_prots(euk, replace = FALSE)
neg <- split_prots(euk_not, replace = FALSE)
neg[["valid"]] <- neg[["valid"]][1L:length(pos[["valid"]])]
neg[["test"]] <- neg[["test"]][1L:length(pos[["test"]])]


ts <- calc_t(euk[pos[["train"]]])

sfInit( parallel=TRUE, cpus=4 )
sfLibrary(depmixS4)

sfExport("test_protein", "find_nhc", "degenerate", "start_model", 
         "aa5", "ts", "euk", "euk_not", "pos", "neg")
pos_valid <- find_signal(euk[pos[["valid"]]], ts, paralell = TRUE)
neg_valid <- find_signal(euk_not[neg[["valid"]]], ts, paralell = TRUE)
pos_test <- find_signal(euk[pos[["test"]]], ts, paralell = TRUE)
neg_test <- find_signal(euk_not[neg[["test"]]], ts, paralell = TRUE)

sfStop()

valid <- rbind(list_to_df(pos_valid, "pos"), list_to_df(neg_valid, "neg"))
test <- rbind(list_to_df(pos_test, "pos"), list_to_df(neg_test, "neg"))

model <- ksvm(tar ~ ., valid, C = 1, prob.model = TRUE)
list(test = predict(model, test, type = "probabilities")[,1],
     bad_pos = sum(sapply(pos_valid, class) != "matrix"),
     bad_neg = sum(sapply(neg_valid, class) != "matrix"))


res_pc <- princomp(x = res[-5])
plot(res_pc$scores[, 1], res_pc$scores[, 2], cex = 0, xlab = "PCA_1", ylab = "PCA_2", main = "Wyniki PCA")
points(res_pc$scores[res[5] == "neg", 1], res_pc$scores[res[5] == "neg", 2])
points(res_pc$scores[res[5] == "pos", 1], res_pc$scores[res[5] == "pos", 2], col = "red")

dump("res",file="dane_pca.Rdat")
i = 4
titles <- c("Zgodnoœæ calkowita", "Zgodnoœæ - N-region", "Zgodnoœæ - H-region",
            "Zgodnoœæ - C-region")

dens_neg <- density(result_filtered_neg[i, ])
dens_pos <- density(result_filtered_pos[i, ])     
max_y <- ifelse(max(dens_neg[["y"]]) >= max(dens_pos[["y"]]), max(dens_neg[["y"]]), max(dens_pos[["y"]]))
plot(c(0, 1), c(0, max_y), main = titles[i], cex = 0)
lines(dens_neg, lwd = 2)
lines(dens_pos, col = "red", lty = "dashed")
legend("topleft", c("neg", "pos"), lty = c("solid", "dashed"), lwd = c(2, 1), col = c("black", "red"))

dump(c("result_filtered_neg", "result_filtered_pos"), file="densities.Rdat")
pred <- prediction(predict(model, test, type = "probabilities")[,1],
                   abs(as.numeric(test[["tar"]]) - 2))
save("pred",file="pred.Rdat")
plot(performance(pred, "tpr", "fpr"), main = "ROC curve")
text(0.8, 0.2, paste0("AUC = ", round(performance(pred, "auc")@y.values[[1]], 2)))

# for (i in euk_not[ids_n[1:50]]) {
#   print(" ")
#   print(" ")
#   print(" ")
#   print(attr(i,"name"))
#   print(attr(i,"name"))
#   print(attr(i,"name"))
#   print(" ")
#   print(" ")
#   print(" ")
#   find_signal(list(el1=i), ts)
# }
# 
# euk_not["sp|P11088|FILA_MOUSE"]
# 
# problem <- vapply(pot_sigs, function(cs) find_nhc(euk_not["sp|P11088|FILA_MOUSE"], c(1, cs)), c(0, 0, 0, 0))
# 
# for (cs in pot_sigs)
#   find_nhc(euk_not["sp|P11088|FILA_MOUSE"], c(1, cs))
# 
# find_nhc(euk_not["sp|P11088|FILA_MOUSE"], c(1, 20))


# procent_rozpoznania <- NULL
# testowane_bialka <- sample(1:length(analized_sequences),100, replace=FALSE)
# for(numer_probki in testowane_bialka){
#   probka <- as.numeric(degenerate(analized_sequences[[numer_probki]], aa5)[1:(euk_nhc[numer_probki,4]+30)])
#   fitted.model <- start_model(probka, t1, t2, t3)
#   viterbi_path <- fitted.model@posterior[1L:euk_nhc[numer_probki,4],1]
#   expected <- c(rep(1,euk_nhc[numer_probki,2]-1),rep(2,euk_nhc[numer_probki,3]-euk_nhc[numer_probki,2]),rep(3,euk_nhc[numer_probki,4]-euk_nhc[numer_probki,3]+1))
#   procent_rozpoznania <- c(procent_rozpoznania, sum(viterbi_path==expected)/length(viterbi_path))
# }
# 
# #bardzo prosta statystyka opisowa ;)
# mean(procent_rozpoznania)
# 
# wynik <- cbind(probka, fitted.model@posterior)