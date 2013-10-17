#próba znalezienia profili białkowych
require(seqinr)
source("czestosc_w_grupach.R")
zip <- function(x,y){
  return(x[1:y])
}
cutSignalSeqs <- mapply(zip, speuk[nonZero], dl)
group=aa1
degeneratedSignalSeqs <- lapply(cutSignalSeqs, function(x) degenerate(x, group))
length(degeneratedSignalSeqs)
#we analyze what is the length of signals seqs
N <- 20
sequencesLongerThanN <- Filter(function(x) length(x)>=N, degeneratedSignalSeqs)
sequencesLongerThanN <- lapply(sequencesLongerThanN, function(x) x[1:N])

m <- matrix(unlist(sequencesLongerThanN), ncol=N,byrow=FALSE)

computed.consensus <- consensus(m, method ="threshold",threshold=0)
computed.consensus

con(unlist(sequencesLongerThanN), length(group))