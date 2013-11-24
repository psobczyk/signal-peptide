# transitions -----------------------

#calculates transition probability between aminoacids with specific distance, helper function
calc_transit_by <- function(seq, distance, AA_names) {
  trans <- lapply(AA_names, function(aminoacid) {
    ind <- which(seq == aminoacid) + distance
    ind <- ind[ind < length(seq)]
    if (length(ind) == 0) {
      tmp_table <- rep(0, length(AA_names))
      names(tmp_table) <- AA_names
      tmp_table
    } else {
      table(seq[ind])/length(ind)
    }
  })
  
  trans <- do.call(rbind, trans)
  rownames(trans) <- AA_names
  trans
}

#calculates transition probability between aminoacids
calc_transit <- function(seq, AA_names) {
  lapply(1L:5, function(distance) calc_transit_by(seq, distance, AA_names))
}



protein <- spplant[[1]]
s_start <- 1
s_stop <- 40

seq <- protein[s_start:s_stop]

AA <- a()[-1]

seq <- factor(toupper(seq), level = AA)

calc_transit(seq, AA)

#works also for degenerated sequences

aa1 = list(nonpolar.aliphatic = c("G", "A", "P", "V", "L", "I", "M"), 
           aromatic = c("F", "W", "Y"),
           positively.charged=c("K", "R", "H"), 
           negatively.charged=c("D", "E"), 
           polar.uncharged=c("S", "T", "C", "N", "Q"))

deg_seq <- degenerate(as.character(seq), aa1)
deg_seq <- factor(deg_seq, level = as.character(1L:length(aa1)))
calc_transit(deg_seq, as.character(1L:length(aa1)))




