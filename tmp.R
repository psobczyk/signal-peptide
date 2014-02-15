setwd("D:/pobieranie2/uniprot_sprot.dat")

#number of lines: 54392766
#select: (keyword:signal) AND reviewed:yes AND created:[1950 TO 1994] AND reviewed:yes AND taxonomy:"Eukaryota [2759]"

all_lines <- readLines("dane.txt")
prot_ids <- grep("\\<ID   ", all_lines)
seqs_start <- grep("\\<SQ   ", all_lines) + 1
seqs_end <- grep("//", all_lines)
webs <- grep("http://", all_lines)


#removes proteins with probable or potential signal peptides
#removes proteins without cleavage site for signalase
remove_unsure <- function(all_lines, prot_ids) {
  signals <- grep("SIGNAL", all_lines)
  true_signals <- signals[grep("FT", all_lines[signals])]
  not_probable <- true_signals[-grep("Probable", all_lines[true_signals])]
  not_potential <- not_probable[-grep("Potential", all_lines[not_probable])]
  cleaved <- not_potential[-grep("Not cleaved.", all_lines[not_potential])]
  all_ids <- sort(c(prot_ids, cleaved), method = "quick")
  all_ids[which(all_ids %in% cleaved) - 1]
}

tmp <- remove_unsure(all_lines, prot_ids)


tmp <- grep("OG   ", all_lines)
#proteins not encoded in nucleus
#notnuc_info <- which(pbsapply(all_lines, function(line) substr(line, 0, 2) == "OG"))
notnuc_info <- grep("OG   ", all_lines)
all_ids <- sort(c(prot_ids, notnuc_info), method = "quick")
prot_ids2 <- all_ids[which(all_ids %in% notnuc_info) - 1]

test_prot <- pblapply(prot_ids2, function(i) {
  name <- strsplit(all_lines[i], "   ")[[1]][2]
  while(substr(all_lines[i], 0, 2) != "SQ")
    i <- i +1
  start_seq <- i + 1
  
  while(substr(all_lines[i], 0, 2) != "//")
    i <- i +1
  end_seq <- i - 1
  seq <- strsplit(gsub(" ", "", paste0(all_lines[start_seq:end_seq], collapse = "")), "")[[1]]
  class(seq) <- "SeqFastaAA"
  attr(seq, "name") <- name
  seq
})