library(seqinr)

aa <- a()[-1]

# ortogonal coding ---------------------------------
#generate ort_codes
ort_codes <- sapply(1:20, function(i) {
  res <- rep(0, 20)
  res[i] <- 1
  res
})

rownames(ort_codes) <- aa

#code ort_codes
encode_ort <- function(seq, codes)
  as.vector(codes[seq, ])


# test encoding
# tmp <- matrix(sample(aa, 25*9, replace = TRUE), nrow = 25)
# t(apply(tmp, 1, function(i)
#   encode_ort(i, ort_codes)))


# k-tuple encoding ---------------------------------

create_features <- function(n, u) {
  #creates a vector of all possible grams
  #u - unigrams
  #n - size of gram (i.e. 2-grams: "AA", "AC", ...)
  grid_list <- lapply(1L:n, function(i) u)
  apply(expand.grid(grid_list), 1, function(x)
    paste(x, collapse = "."))
}

create_grams <- function(seq, n, d) {
  #creates grams from given function
  #n - size of gram (i.e. 2-grams: "AA", "AC", ...)
  #d - distance between two consecutive aa
  ind <- lapply(1:n, function(i) 
    (1 + i - 1):(length(seq) - n + i))
  
  if(d > 0) {
    ind[-1] <- lapply(ind[-1], function(i) i + d)
    not_taken <- ind[[1]][(length(ind[[1]]) - d + 1):length(ind[[1]])]
    ind <- lapply(ind, function(i) i[-not_taken])
  }
  
  pair_matrix <- do.call(cbind, lapply(ind, function(i) seq[i]))
  apply(pair_matrix, 1, function(x) 
    paste(x, collapse="."))
}

count_grams <- function(seq, feature_list, n, d = 0, scale = FALSE) {
  #feature list(list of possible n-grams) is outside, because count_grams is meant to
  #be inside the loop
  if (n > 1) {
    grams <- create_grams(seq, n, d)
  } else {
    grams <- seq
  }
  res <- sapply(feature_list, function(i)
    sum(i == grams))
  if (scale)
    res <- res/(length(seq) - n - 1)
  res
}

get_kmer <- function(list_prot, k, d = 0) {
  #get k-mer of aminoacids within distance d from cleaves place
  if (k %% 2 != 0) {
    k1 <- k2 <- floor(k/2)
  } else {
    k2 <- k/2
    k1 <- k2 - 1
  }
  t(sapply(list_prot, function(ith_prot) 
    ith_prot[(attr(ith_prot, "sig")[2] + d - k1):(attr(ith_prot, "sig")[2] + d + k2)]))
}


calc_unigrams <- function(seqs) {
  #helper function for calculating unigrams
  if (class(seqs) == "matrix") {
    unigram <- apply(seqs, 1, function(i) data.frame(table(factor(i, levels = aa))))
    len <- nrow(seqs)
  } 
  if (class(seqs) == "list") {
    unigram <- lapply(seqs, function(i) data.frame(table(factor(i, levels = aa))))
    len <- length(seqs)
  }
  unigram <- do.call(rbind, unigram)
  agg_unigram <- aggregate(Freq ~ Var1, unigram, sum)
  data.frame(AA = agg_unigram[[1]], freq = agg_unigram[[2]]/len)
}

