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
tmp = calc_transit(deg_seq, as.character(1L:length(aa1)))[[1]]
bar_names <- as.vector(sapply(1:5, function(x) paste0(1:5, x)))

signals_length <- t(vapply(speuk, function(seq) attr(seq, "sig"), c(0, 0)))
plot(density(signals_length[signals_length[,1] == 1,2]), main = "")

deg_levels <- as.character(1L:length(aa1))
signals1 <- vapply(speuk[signals_length[,1] == 1], function(seq) {
  deg_seq <- degenerate(seq[attr(seq, "sig")[1]:attr(seq, "sig")[2]], aa1)
  deg_seq <- factor(deg_seq, level = deg_levels)
  as.vector(calc_transit_by(deg_seq, 1, deg_levels))
}, rep(0, length(aa1)^2))

nonsignals1 <- vapply(speuk[signals_length[,1] == 0], function(seq) {
  deg_seq <- degenerate(seq[1:68], aa1)
  deg_seq <- factor(deg_seq, level = deg_levels)
  as.vector(calc_transit_by(deg_seq, 1, deg_levels))
}, rep(0, length(aa1)^2))


barplot(rbind(rowMeans(signals1), rowMeans(nonsignals1)), beside = TRUE, col = c("blue", "red"),
        names.arg = bar_names)

signals2 <- vapply(speuk[signals_length[,1] == 1], function(seq) {
  deg_seq <- degenerate(seq[attr(seq, "sig")[1]:attr(seq, "sig")[2]], aa1)
  deg_seq <- factor(deg_seq, level = deg_levels)
  as.vector(calc_transit_by(deg_seq, 2, deg_levels))
}, rep(0, length(aa1)^2))

nonsignals2 <- vapply(speuk[signals_length[,1] == 0], function(seq) {
  deg_seq <- degenerate(seq[1:68], aa1)
  deg_seq <- factor(deg_seq, level = deg_levels)
  as.vector(calc_transit_by(deg_seq, 2, deg_levels))
}, rep(0, length(aa1)^2))

barplot(rbind(rowMeans(signals2), rowMeans(nonsignals2)), beside = TRUE, col = c("blue", "red"), 
        names.arg = bar_names)

signals3 <- vapply(speuk[signals_length[,1] == 1], function(seq) {
  deg_seq <- degenerate(seq[attr(seq, "sig")[1]:attr(seq, "sig")[2]], aa1)
  deg_seq <- factor(deg_seq, level = deg_levels)
  as.vector(calc_transit_by(deg_seq, 3, deg_levels))
}, rep(0, length(aa1)^2))

nonsignals3 <- vapply(speuk[signals_length[,1] == 0], function(seq) {
  deg_seq <- degenerate(seq[1:68], aa1)
  deg_seq <- factor(deg_seq, level = deg_levels)
  as.vector(calc_transit_by(deg_seq, 3, deg_levels))
}, rep(0, length(aa1)^2))

barplot(rbind(rowMeans(signals3), rowMeans(nonsignals3)), beside = TRUE, col = c("blue", "red"), 
        names.arg = bar_names)

signals4 <- vapply(speuk[signals_length[,1] == 1], function(seq) {
  deg_seq <- degenerate(seq[attr(seq, "sig")[1]:attr(seq, "sig")[2]], aa1)
  deg_seq <- factor(deg_seq, level = deg_levels)
  as.vector(calc_transit_by(deg_seq, 4, deg_levels))
}, rep(0, length(aa1)^2))

nonsignals4 <- vapply(speuk[signals_length[,1] == 0], function(seq) {
  deg_seq <- degenerate(seq[1:68], aa1)
  deg_seq <- factor(deg_seq, level = deg_levels)
  as.vector(calc_transit_by(deg_seq, 4, deg_levels))
}, rep(0, length(aa1)^2))

barplot(rbind(rowMeans(signals4), rowMeans(nonsignals4)), beside = TRUE, col = c("blue", "red"), 
        names.arg = bar_names)


