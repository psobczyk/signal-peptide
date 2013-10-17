#program analizujący sekwencje sygnałowe
source("czestosc_w_grupach.R")

family <- c( "ctp", "mtp", "notp") #spplant
groups <- c("aa1", "aa2", "aa3", "aa4")

aux <- function(x,y){
  dists <- distribution(eval(as.name(x)),eval(as.name(y)))
  chisq.test(dists[[1]], dists[[2]])$p.value
} 
result <- NULL
result[1] <- aux("spplant", "aa1")
result[2] <- aux("spplant", "aa2")
result[3] <- aux("spplant", "aa3")
result[4] <- aux("spplant", "aa4")
result[5] <- aux("ctp", "aa1")
result[6] <- aux("ctp", "aa2")
result[7] <- aux("ctp", "aa3")
result[8] <- aux("ctp", "aa4")
result[9] <- aux("mtp", "aa1")
result[10] <- aux("mtp", "aa2")
result[11] <- aux("mtp", "aa3")
result[12] <- aux("mtp", "aa4")
result[13] <- aux("notp", "aa1")
result[14] <- aux("notp", "aa2")
result[15] <- aux("notp", "aa3")
result[16] <- aux("notp", "aa4")
