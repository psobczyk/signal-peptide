# constants ----------------------------

aa <- a()[-1]

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

aa6 <- list(`1` = c("A", "C", "Q", "E"), 
            `2` = c("R", "H",  "K"), 
            `3` = c("N", "D", "G", "P", "S", "T"), 
            `4` = c("I", "L",  "M", "F", "W", "Y", "V"))

#generate ort_codes
ort_codes <- sapply(1:20, function(i) {
  res <- rep(0, 20)
  res[i] <- 1
  res
})

rownames(ort_codes) <- aa

# data ----------------

euk <- read_uniprot("euk.txt", euk = TRUE)
euk_not <- read.fasta("euk_not.fasta", seqtype = "AA")
euk_not <- euk_not[-get_atyp(euk_not)]
euk_not <- euk_not[-which(sapply(euk_not, length) > 4000)]
euk_not <- euk_not[-which(sapply(euk_not, length) < 50)]


pos <- split_prots(euk, replace = FALSE)
neg <- split_prots(euk_not, replace = FALSE)
neg[["valid"]] <- neg[["valid"]][1L:length(pos[["valid"]])]
neg[["test"]] <- neg[["test"]][1L:length(pos[["test"]])]

# searching cleave region ---------------------

#########
kmer = 9 
prot_list = euk[pos[["train"]]]
########

cleaves <- get_kmer(prot_list, kmer)

#parts of signal peptide without cleaving region
n_cleaves1 <- get_kmer(prot_list[sapply(prot_list, function(ith_prot) 
  attr(ith_prot, "sig")[2] > 16)], kmer, d = -round(kmer * 1.3, 0))

#mature protein
n_cleaves2 <- get_kmer(prot_list, kmer, d = round(kmer * 1.3, 0))
n_cleaves <- rbind(n_cleaves1, n_cleaves2)
colnames(cleaves) <- paste0("P", 1L:kmer)
dat <- rbind(cbind(cleaves, et = rep("pos", nrow(cleaves))),
             cbind(n_cleaves, et = rep("neg", nrow(n_cleaves))))

edtrain <- cbind(t(apply(dat[, -ncol(dat)], 1, function(i)
  encode_ort(i, ort_codes))), et = dat[, ncol(dat)])

#training ----------------------
sfInit(parallel=TRUE, cpus=6)
ctrl <- trainControl(method = "repeatedcv",
                     repeats = 3,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     allowParallel = TRUE)

svm_fit <- train(et ~ .,
                 data = edtrain,
                 method = "svmLinear",
                 trace = FALSE,
                 maxit = 50,
                 trControl = ctrl,
                 metric = "ROC")
sfStop()



# testing signal peptide ---------------------

ts <- calc_t(euk[pos[["train"]]])

sfInit(parallel=TRUE, cpus=6)
sfLibrary(depmixS4)
sfLibrary(hsmm)

sfExport("encode_ort", "svm_fit", "ort_codes", "test_protein", "find_nhc", "find_cleave", "degenerate", "start_model", "start_model2", 
         "aa5", "ts", "euk", "euk_not", "pos", "neg")
pos_valid <- find_signal(euk[pos[["valid"]]], ts, model = svm_fit$finalModel, paralell = TRUE)
neg_valid <- find_signal(euk_not[neg[["valid"]]], ts, model = svm_fit$finalModel, paralell = TRUE)
pos_test <- find_signal(euk[pos[["test"]]], ts, model = svm_fit$finalModel, paralell = TRUE)
neg_test <- find_signal(euk_not[neg[["test"]]], ts, model = svm_fit$finalModel, paralell = TRUE)

sfStop()

pos_valid <- find_signal(euk[pos[["valid"]]], ts, model = svm_fit$finalModel, paralell = FALSE)
neg_valid <- find_signal(euk_not[neg[["valid"]]][5:6], ts, model = svm_fit$finalModel, paralell = FALSE)
pos_test <- find_signal(euk[pos[["test"]]], ts, model = svm_fit$finalModel, paralell = FALSE)
neg_test <- find_signal(euk_not[neg[["test"]]], ts, model = svm_fit$finalModel, paralell = FALSE)



valid <- rbind(list_to_df(pos_valid, "pos"), list_to_df(neg_valid, "neg"))
test <- rbind(list_to_df(pos_test, "pos"), list_to_df(neg_test, "neg"))

model <- ksvm(tar ~ ., valid, C = 1, prob.model = TRUE)
list(test = predict(model, test, type = "probabilities")[,1],
     bad_pos = sum(sapply(pos_valid, class) != "matrix"),
     bad_neg = sum(sapply(neg_valid, class) != "matrix"))

pred <- prediction(predict(model, test, type = "probabilities")[,1],
                   abs(as.numeric(test[["tar"]]) - 2))
plot(performance(pred, "tpr", "fpr"), main = "ROC curve")
text(0.8, 0.2, paste0("AUC = ", round(performance(pred, "auc")@y.values[[1]], 2)))
