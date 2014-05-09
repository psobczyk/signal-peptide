#viterbi

viterbi <- function(probka, pipar, tp, od){
  
  viterbi <- matrix(nrow=length(probka), ncol=4)
  psi <- matrix(nrow=length(probka), ncol=4)
  for(j in 1:4){
    viterbi[1,j]=log(pipar[j])+log(od[j,probka[1]])
    psi[1,j] = 1
  }  
  for(i in 2:length(probka)){
    for(j in 1:4){
      max=-1000
      max.i = 0
      for(k in 1:4){
        if(viterbi[i-1,k]+log(tp[k,j]) + log(od[j,probka[i]])>max){
          max = viterbi[i-1,k]+log(tp[k,j]) + log(od[j,probka[i]])
          max.i = k
        }
      }
      viterbi[i,j]=max
      psi[i,j]=max.i
    }  
  }
  path = NULL
  path[length(probka)] = which.max(viterbi[length(probka),])
  for(i in (length(probka)-1):1){
    path[i] = psi[i, path[i+1]]
  }
  return(list(path =path,
              viterbi = viterbi,
              psi = psi))
}

