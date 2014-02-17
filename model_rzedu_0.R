#model markowa rzędu 0 (czestości wystepowania aminikwasów)

source("main.R")

#uruchamiamy ukryte łańcuchy
require(depmixS4)
source("run_model.R")

procent_rozpoznania <- NULL
testowane_bialka <- sample(1:length(analized_sequences),50, replace=FALSE)
for(numer_probki in testowane_bialka){
  probka <- as.numeric(degenerate(analized_sequences[[numer_probki]], aa5)[1:(all_nhc[numer_probki,4]+10)])
  fitted.model <- uruchom_model(probka)
  viterbi_path <- fitted.model@posterior[1:all_nhc[numer_probki,4],1]
  expected <- c(rep(1,all_nhc[numer_probki,2]-1),rep(2,all_nhc[numer_probki,3]-all_nhc[numer_probki,2]),rep(3,all_nhc[numer_probki,4]-all_nhc[numer_probki,3]+1))
  procent_rozpoznania <- c(procent_rozpoznania, sum(viterbi_path==expected)/length(viterbi_path))
}

#bardzo prosta statystyka opisowa ----------------------------
mean(procent_rozpoznania)

# fa ----------------------------
hist(procent_rozpoznania,probability=T,breaks=20)

#wykresy dla hmm ----------------
a <- fitted.model@posterior
par(mar = c(5,5,2,5))
plot(a[,2],type="l", col="red", ylab="marginal probilities")
lines(a[,3], col="blue")
lines(a[,4], col="green")
par(new=TRUE)
barplot(a[,1],,col="pink",xaxt="n",yaxt="n",xlab="",ylab="state",density=4)
axis(side=4,at=c(0,1,2,3), labels=c("","one", "two", "three"),line=F,tick=F)

height <- t(t(c(1,-1,1)))
bardensity <- t(t(c(10,10,0)))
barangle <- t(t(c(45,135,0)))
barplot(height, density = bardensity, angle = barangle) 


s <- strsplit("AWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW",split="")
probka <- as.numeric(degenerate(unlist(s), aa5))
fitted.model <- uruchom_model(probka)
viterbi_path <- fitted.model@posterior
