#' Plotting signal peptide hmm results
#' 
#' Plots the Viterbi path and cut site
#' 
#' 

plot.Viterbi.results <- function(viterbi_path, protein){
  length = length(viterbi_path)
  circle.y <- length/4
  c.site <- ifelse(any(viterbi_path==3), max(which(viterbi_path==3)), length)
  par(mar=c(2,1,4,1)) 
  plot(x=1:length, type="n", ylim=c(0,circle.y*2), xlim=c(1, length), ylab="", axes=F, main="Predicted regions for aminoacids in peptide")
  #axis(side=1, at=1:length, labels=protein[1:length])
  for(j in 1:length){
    if(viterbi_path[j]==1)
      circle(j,circle.y, 0.5, col="slateblue1")
    if(viterbi_path[j]==2)
      circle(j,circle.y, 0.5, col="pink2")
    if(viterbi_path[j]==3)
      circle(j,circle.y, 0.5, col="ivory2")
    if(viterbi_path[j]==4)
      circle(j,circle.y, 0.5, col="seagreen1")
    #rect(j-0.5,circle.y-0.5, j+0.5, circle.y+0.5, col="seagreen1")
  }
  abline(v=c.site+0.5, lwd =3, col="red")
  text(x=c.site+0.5, y=circle.y+3, labels="cut", cex=2,pos=2)
  text(x=c.site, y=circle.y-3, pos=1, cex=2, labels="Aminoacids")
  legend(x="bottomright", legend=c("n-region", "h-region", "c-region", "cleavage site", "mature protein"),
         fill=c("slateblue1", "pink2", "ivory2", "red",  "seagreen1"))
  text(x=1:length, y=circle.y-1.5, labels=protein[1:length])
}





