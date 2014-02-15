#positive set
#select: (keyword:signal) AND reviewed:yes AND created:[2012 TO 2013] AND taxonomy:"Eukaryota [2759]" AND annotation:(type:signal confidence:experimental)
test_pos <- read_uniprot("test_pos.txt", euk = TRUE)
write.fasta(test_pos, names = names(test_pos), file.out = "test_pos.fasta")

#negative set
#NOT (keyword:signal) AND reviewed:yes AND created:[2012 TO 2013] AND taxonomy:"Eukaryota [2759]"
test_neg <- read.fasta("test_neg_big.fasta")
neg_id <- sample(length(test_neg), 201)
test_neg <- test_neg[neg_id]
test_neg <- test_neg[-which(sapply(test_neg, length) > 4000)]
write.fasta(test_neg, names = names(test_neg), file.out = "test_neg.fasta")

true_labels <- c(rep(1, length(test_pos)), rep(0, length(test_neg)))

sigp2 <- c(read_signalp2("signalp2_pos.txt"), read_signalp2("signalp2_neg.txt"))
sigp4 <- c(read_signalp4("signalp4_pos.txt"), read_signalp4("signalp4_neg.txt"))

pred_sigp2 <- prediction(sigp2, true_labels)
pred_sigp4 <- prediction(sigp4, true_labels)

perf_sigp2 <- performance(pred_sigp2, "tpr", "fpr")
perf_sigp4 <- performance(pred_sigp4, "tpr", "fpr")
perf_piotr <- performance(pred, "tpr", "fpr")
plot(perf_sigp2, main = "ROC curve")
lines(slot(perf_sigp4, "x.values")[[1]], slot(perf_sigp4, "y.values")[[1]], col = "red")
lines(slot(perf_piotr, "x.values")[[1]], slot(perf_piotr, "y.values")[[1]], col = "green", lty = "dashed",
      lwd = 2)
legend("bottomright", c("SignalP 2", "SignalP 4", "Piotr"), 
       lty = c("solid", "solid", "dashed"), 
       lwd = c(1, 1, 2),
       col = c("black", "red", "green"))
