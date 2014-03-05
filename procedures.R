# constants ----------------------------

aa5 = list(positively.charged=c("K", "R", "H"),
           hydrofobic=c("V","I","L","M","F","W","C"),
           polar.uncharged=c("S", "T", "N", "Q"),
           rest=c("D","E","A","P","Y","G"))

# data ----------------

euk <- read_uniprot("euk.txt", euk = TRUE)
euk_not <- read.fasta("euk_not.fasta", seqtype = "AA")
euk_not <- euk_not[-get_atyp(euk_not)]

pos <- split_prots(euk, replace = FALSE)
neg <- split_prots(euk_not, replace = FALSE)
neg[["valid"]] <- neg[["valid"]][1L:length(pos[["valid"]])]
neg[["test"]] <- neg[["test"]][1L:length(pos[["test"]])]

# testing algorithm ---------------------

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

pred <- prediction(predict(model, test, type = "probabilities")[,1],
                   abs(as.numeric(test[["tar"]]) - 2))
plot(performance(pred, "tpr", "fpr"), main = "ROC curve")
text(0.8, 0.2, paste0("AUC = ", round(performance(pred, "auc")@y.values[[1]], 2)))
