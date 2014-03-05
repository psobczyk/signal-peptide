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

# values ----------------------------

aa5 = list(positively.charged=c("K", "R", "H"),
           hydrofobic=c("V","I","L","M","F","W","C"),
           polar.uncharged=c("S", "T", "N", "Q"),
           rest=c("D","E","A","P","Y","G"))

transition1 <- make_transition_matrix(n_region, aa5)
transition2 <- make_transition_matrix(h_region, aa5)
transition3 <- make_transition_matrix(c_region, aa5)
transition4 <- make_transition_matrix(rest, aa5)



# libraries ----------------------------

library(depmixS4)
library(dhmm)
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
sfLibrary(dhmm)

sfExport("test_protein", "find_nhc", "degenerate", "start_model", "start_model2", 
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
titles <- c("Zgodno?? calkowita", "Zgodno?? - N-region", "Zgodno?? - H-region",
            "Zgodno?? - C-region")

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