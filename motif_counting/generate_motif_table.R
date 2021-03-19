get_motifs <- function() {
  library(seqinr)
  
  RxLR <- "^\\w{1,56}R\\wLR"
  RxLR_EER <- "^\\w{1,56}R\\wLR\\w{1,40}[D|E][D|E][K|R]"
  
  QGH_xLR <- "^\\w{1,56}[Q|G|H]\\wLR"
  QGH_xLR_EER <- "^\\w{1,56}[Q|G|H]\\wLR\\w{1,40}[D|E][D|E][K|R]"
  
  RxL_QKG <- "^\\w{1,56}R\\wL[Q|K|G]"
  RxL_QKG_EER <- "^\\w{1,56}R\\wL[Q|K|G]\\w{1,40}[D|E][D|E][K|R]"
  
  EER <- "^\\w{1,97}[D|E][D|E][K|R]"
  
  # old: reg.pat<-"^\\w{20}\\w{1,80}R\\wLR"
  # extramotif: [RKHGQ][X][LMYFIV][RNK]
  
  motifs <- c(RxLR, RxLR_EER, QGH_xLR, QGH_xLR_EER, RxL_QKG, RxL_QKG_EER, EER)
  names <- c("RxLR", "RxLR_EER", "QGH_xLR", "QGH_xLR_EER", "RxL_QKG", "RxL_QKG_EER",
             "EER")
  
 # motifs <- c(RxLR, RxLR_EER, GxLR)
 # names <- c("RxLR", "RxLREER", "GxLR")
  
  return(list("motifs"=motifs, "names"=names))
} # get_motifs function

process_secretome <- function(motifs, names, fasta.file) {
  
  ORF <- read.fasta(fasta.file)
  SeqIDs <- names(ORF)
  
  filename <- rep(strsplit(fasta.file, "/")[[1]][9], length(SeqIDs))
  results_DF <- data.frame(filename, SeqIDs,  stringsAsFactors = FALSE)
  
  #initialize data frame
  colnames <- c("filename", "seqID")
  for (name in names) {
    colnames <- c(colnames, name)
  }
  
  for (i in 3:length(colnames)) {
    hit_List <- rep("N", length(SeqIDs))
    results_DF <- data.frame(results_DF, hit_List, check.names = FALSE, 
                             stringsAsFactors = FALSE)
  }
  
  colnames(results_DF) <- colnames
  #end initializing data frame
  
  
  #start processing ORF:

  for (i in 1:length(motifs)) {
    lenREG <- 0
    
    res <- 0
    res <- try(REGEX <- regex.search(ORF, motif = "custom", motifs[i]))
    
    if (substr(res[1], 1, 5) == "Error")
      lenREG <- 0
    else 
      lenREG <- length(REGEX)
    
    for (ID in 1:lenREG) {
      sequence <- attr(REGEX[ID], "name")
      results_DF[results_DF$seqID==sequence, names[i]] <- "Y"
    }

  } #loop through each motif, append information to avg_* vectors
  
  return(results_DF)
}

main <- function(secretome_path) {
  info <- get_motifs()
  genomeName <- strsplit(secretome_path, "/")[[1]][8]
  file_names <- c(rep(genomeName, length(info$names)))
  
  normal_DF <- process_secretome(info$motifs, info$names,secretome_path)
  
  #outFile <- paste("EffectR_results_", genomeName, ".txt", sep="")
  outFile <- "/share/rwmwork/mjnur/effectR/motif_redo_201810/results_20190123/20210319_EffectR_results_NR_cleaved.txt"
  cat("writing to outfile, " , outFile, "\n")
  
  isBremia <- grepl("lac", secretome_path)

  if (isBremia)
    write.table(normal_DF, outFile, row.names = FALSE, col.names = TRUE, sep=",", append=FALSE)
  else
    write.table(normal_DF, outFile, row.names = FALSE, col.names = FALSE, sep=",", append=TRUE)
}


args <- commandArgs(trailingOnly = TRUE)
secretome_path <- args[1]

main(secretome_path)
