#' Plotting signal peptide hmm results
#' 
#' Plots the Viterbi path and cut site
#' 
#' 

plot.Viterbi.results <- function(viterbi_path, protein){
  length = length(viterbi_path)
  c.site <- ifelse(any(viterbi_path==3), max(which(viterbi_path==3)), length)
  plot(x=1:length, type="n", rep(1, length), xaxt='n', yaxt='n', ylim=c(0,1), xlim=c(1, length), ylab="", xlab="Aminoacids")
  axis(side=1, at=1:length, labels=protein[1:length])
  for(j in 1:length){
    if(viterbi_path[j]==1)
      rect(j-0.5,0,j+0.5, 1, col="darkorchid1")
    if(viterbi_path[j]==2)
      rect(j-0.5, 0, j+0.5, 1, col="honeydew1")
    if(viterbi_path[j]==3)
      rect(j-0.5, 0,j+0.5, 1, col="aliceblue")
    if(viterbi_path[j]==4)
      rect(j-0.5,0,j+0.5, 1, col="beige")
  }
  abline(v=c.site+0.5, lwd =3, col="red")
  text(x=c.site+0.5, y=1, labels="cut", cex=2)
  legend(x="bottomright", legend=c("n-region", "h-region", "c-region", "cleavage site", "rest of protein"),
         fill=c("darkorchid1", "honeydew1", "red", "aliceblue", "beige"))
}
