find_nhc <- function(protein) {
  sig <- protein[attr(protein, "sig")[1]:attr(protein, "sig")[2]]
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
  
  noh <- 0
  while(noh == 0) {
    start_h <- start_h + 1
    noh <- noh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), 1, 0)
  }
  c(start_n = attr(protein, "sig")[1], start_h = start_h, start_c = start_c, cs = attr(protein, "sig")[2])
}



with_sig <- vapply(speuk, function(protein) attr(protein, "sig")[1] != 0, FALSE)
all_nhc <- t(vapply(speuk[with_sig], find_nhc, rep(0, 4)))
len_c <- all_nhc[, "cs"] - all_nhc[, "start_c"]
len_h <- all_nhc[, "start_c"] - all_nhc[, "start_h"]
len_n <- all_nhc[, "start_h"] - all_nhc[, "start_n"]
lengths <- matrix(c(len_n, len_h, len_c), ncol = 3)
colnames(lengths) <- c("n", "h", "c")
boxplot(lengths)
