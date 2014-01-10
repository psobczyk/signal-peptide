#adding memory to model
source("main.R")

encoding <- function(x){
  a <- x%/%10
  b <- x%%10
  b+4*(a-1)
}

make_transition_matrix <- function(region, aa_groups){
  b <- degenerate(region, aa_groups)
  tabela <- table(apply(cbind(b[-length(b)],b[-1]), 1, FUN=function(x) encoding(as.numeric(paste(x[1], x[2], sep="")))))
  for(i in 0:3){
    suma <- sum(tabela[(1+i*4):((i+1)*4)])
    for(j in 1:4){
      tabela[i*4+j] <- tabela[i*4+j]/suma
    }
  }
  return(tabela)
}

transition1 <- make_transition_matrix(n_region, aa5)
transition2 <- make_transition_matrix(h_region, aa5)
transition3 <- make_transition_matrix(c_region, aa5)
transition4 <- make_transition_matrix(reszta, aa5)


procent_rozpoznania2 <- NULL
n <- 10
testowane_bialka <- sample(1:length(analized_sequences),n, replace=FALSE)
for(numer_probki in testowane_bialka){
  probka <- (degenerate(analized_sequences[[numer_probki]], aa5)[1:(all_nhc[numer_probki,4]+30)])
  probka <- apply(cbind(probka[-length(probka)],probka[-1]), 1, FUN=function(x) encoding(as.numeric((paste(x[1], x[2], sep="")))))
  fitted.model <- uruchom_model2(probka)
  viterbi_path <- fitted.model@posterior[1:all_nhc[numer_probki,4],1]
  expected <- c(rep(1,all_nhc[numer_probki,2]-1),rep(2,all_nhc[numer_probki,3]-all_nhc[numer_probki,2]),rep(3,all_nhc[numer_probki,4]-all_nhc[numer_probki,3]+1))
  procent_rozpoznania2 <- c(procent_rozpoznania2, sum(viterbi_path==expected)/length(viterbi_path))
}

#bardzo prosta statystyka opisowa ;)
mean(procent_rozpoznania2)

rbind(expected, viterbi_path)

cbind(probka, fitted.model@posterior)