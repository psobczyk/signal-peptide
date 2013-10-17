#packages -----------------------------------------------
chooseCRANmirror(ind = 1)
install.packages("seqinr")
require("seqinr", character.only = TRUE)


# aminoacid groups --------------------------------------
aa1 = list(nonpolar.aliphatic = c("G", "A", "P", "V", "L", "I", "M"), 
           aromatic = c("F", "W", "Y"),
           positively.charged=c("K", "R", "H"), 
           negatively.charged=c("D", "E"), 
           polar.uncharged=c("S", "T", "C", "N", "Q"))
#\cite{timberlake92}

aa2 = list(nonpolar = c("G", "A", "P", "V", "L", "I", "M", "F"), 
           positively.charged=c("K", "R", "H"), 
           negatively.charged=c("D", "E"), 
           polar.uncharged=c("S", "T", "C", "N", "Q", "Y", "W"))
#\cite{patron05}

aa3 = list(monoamino.monocarboxylic = c("G", "A", "P", "V", "L", "I", "M", 
                                        "F", "W", "S", "T", "C", "N", "Q", "Y"), 
           monoamino.dicarboxylic = c("D", "E"), 
           diamino.monocarboxylic = c("K", "R", "H"))
#\cite{devlin92}

aa4 = list(aliphatic = c("G", "A", "V", "L", "I"), sulfur = c("M", "C"), 
           aromatic = c("F", "W", "Y"), 
           neutral = c("S", "T", "N", "Q"), 
           positively.charged=c("K", "R", "H"), 
           negatively.charged=c("D", "E"), 
           imino.acid="P")
#\cite{koolman96}

# functions ---------------------------------------------------------

read.how <- function(file) {
  lines <- readLines(file)
  ind <- c(grep(" ", lines), length(lines) + 1)
  result <- lapply(1L:(length(ind) - 1), function(i) {
    two_seqs <- lines[(ind[i] + 1):(ind[i + 1] - 1)]
    split_seq <- length(two_seqs)/2
    aa <- strsplit(paste0(two_seqs[1L:split_seq], collapse = ""), "")[[1]]
    class(aa) <- "SeqFastaAA"
    attr(aa, "name") <- lines[ind[i]]
    sig <- which(unlist(strsplit(two_seqs[(split_seq + 1):(2*split_seq)], "")) == "S")
    if(length(sig) == 0) {
      attr(aa, "sig") <- c(0, 0)
    } else {
      attr(aa, "sig") <- c(min(sig), max(sig))
    } 
    aa
  }) 
  names(result) <- lines[ind[-length(ind)]]
  result
}

degenerate <- function(seq, aa_group) {
  for (i in 1L:length(aa_group)) {
    seq[seq %in% aa_group[[i]]] <- i
  }
  seq 
}


download.file("http://www.cbs.dtu.dk/services/SignalP/data/euk.2.how", "speuk.how")
download.file("http://www.cbs.dtu.dk/services/SignalP/data/gram+.2.how", "spgplus.how")
download.file("http://www.cbs.dtu.dk/services/SignalP/data/gram-.2.how", "spgminus.how")
download.file("http://www.cbs.dtu.dk/services/TargetP/datasets/planta.SP.269.rr.fasta", "spplant.fasta")
download.file("http://www.cbs.dtu.dk/services/TargetP/datasets/planta.cTP.141.rr.fasta", "ctp.fasta")
download.file("http://www.cbs.dtu.dk/services/TargetP/datasets/planta.mTP.368.rr.fasta", "mtp.fasta")
download.file("http://www.cbs.dtu.dk/services/TargetP/datasets/planta.cyt.87.rr.fasta", "cyt1.fasta")
download.file("http://www.cbs.dtu.dk/services/TargetP/datasets/planta.cyt.108.rr.fasta", "cyt2.fasta")



#secretory signal peptide with the exact localization of signal peptide!
speuk <- read.how("speuk.how") #all eukarionts
spgplus <- read.how("spgplus.how")
spgminus <- read.how("spgminus.how")

#secretory signal peptide without the exact localization of signal peptide!
spplant <- read.fasta("spplant.fasta")

#chloroplast transit peptide
ctp <- read.fasta("ctp.fasta")

#mitochondrion transit peptide
mtp <- read.fasta("mtp.fasta")

#cytoplasmatic proteins, no transit peptide
notp <- c(read.fasta("cyt1.fasta"), read.fasta("cyt2.fasta"))

#degenerate list of sequences
#lapply(speuk[1:10], function(x) degenerate(x, aa1))

#degenerate single sequence
#degenerate(speuk[[50]], aa1)