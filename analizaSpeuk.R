#analiza sekwencji o znanych dlugosciach
require(MASS)

#secretory signal peptide with the exact localization of signal peptide!
speuk <- read.how("speuk.how") #all eukarionts
spgplus <- read.how("spgplus.how")
spgminus <- read.how("spgminus.how")

#now we're fitting gamma dist to empirical lengths
dl <- lapply(speuk, function(x) {y <- attr(x, "sig"); y[2]-y[1]})
nonZero <- dl>0
dl <- unlist(dl[nonZero])
dl <- as.vector(dl)

hist(dl, breaks=40,probability=T)

my.gamma <- list(mean(dl)/var(dl), mean(dl^2)/var(dl))
my.gamma
require(stats4)
x <- seq(0,80)
lines(x, dgamma(x,rate=my.gamma[[1]], shape=my.gamma[[2]]))

#now we're analyzing frequencies of nucleotides in different groups
source("czestosc_w_grupach.R")
zip <- function(x,y){
  return(x[1:y])
}
length(speuk[nonZero])
length(dl)
cutSignalSeqs <- mapply(zip, speuk[nonZero], dl)
pastedSignalSeqs <- unlist(lapply(cutSignalSeqs, function(x) degenerate(x, aa1)))
frequencySignal <- con(pastedSignalSeqs, length(aa1))
cutNonSingalSeqs <- lapply(speuk[!nonZero], function(x) degenerate(x, aa1))
frequencyNonSignal <- con(unlist(cutNonSingalSeqs), length(aa1))
chisq.test(frequencyNonSignal,frequencySignal) #non conclusive - just to give an idea

#now we explore some possible clustering
group <- aa1
freqsSeqs <- lapply(cutSignalSeqs, function(x) con(degenerate(x, group), length(group)))
m <- matrix(unlist(freqsSeqs),ncol=length(group),byrow=FALSE)

wss <- (nrow(m)-1)*sum(apply(m,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(m, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

fit <- kmeans(m, centers=4)
fit$centers

require(pvclust)
pvclust(m)