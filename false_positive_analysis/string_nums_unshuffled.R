library(effectR) ; library(seqinr)

## initialize DataFrame RESULTS_DF and filenames
genomePath = "/share/rwmwork/mjnur/effectR/motif_redo_201810/genomes/NR_cleaved_secretomes/"

names = "/share/rwmwork/mjnur/effectR/motif_redo_201810/genomes/NR_cleaved_secretomes/10_shuffles/filenames"
shuffle_files <- list.files(names, pattern="*.fasta*")
shuffle_files <- substring(shuffle_files, 2)

zeros = rep(0, length(shuffle_files))

results_DF <- data.frame(shuffle_files, zeros,zeros, zeros, stringsAsFactors = FALSE)
colnames(results_DF) = c("genome_file", "RXLR_hits", "EER_hits", "RXLR_and_EER_hits")

RxLR <- "^\\w{1,56}R\\wLR"

QGH_xLR <- "^\\w{1,56}[Q|G|H]\\wLR"

RxL_QKG <- "^\\w{1,56}R\\wL[Q|K|G]"

EER <- "^\\w{1,97}[D|E][D|E][K|R]"

## start calculating the hits in each shuffle file
  
  for (file in shuffle_files) {
    orf <- read.fasta(paste0(genomePath,file))
    
    rxlr_hits = regex.search(orf, motif="custom", "^\\w{1,56}R\\wLR")
    #rxlr_hits = regex.search(orf, motif="custom", "^\\w{1,56}[R|K|H|G|Q]\\w{1,2}[L|M|Y|F|I|V][R|N|K]")

    qghxlr_hits = regex.search(orf, motif="custom", "^\\w{1,56}[Q|G|H]\\wLR")
    rxlqkg_hits = regex.search(orf, motif="custom", "^\\w{1,56}R\\wL[Q|K|G]")
    
    eer_hits = regex.search(orf, motif="custom", "^\\w{1,97}[D|E][D|E][K|R]")
    

    rxlrs_hits = rxlr_hits
    #rxlrs_hits = union(rxlr_hits, qghxlr_hits)
    #rxlrs_hits = union(rxlrs_hits, rxlqkg_hits)



    rxlr_eer_hits = intersect(rxlrs_hits, eer_hits)
    
    results_DF[results_DF$`genome_file`==file, "RXLR_hits"] = length(rxlrs_hits)
    results_DF[results_DF$`genome_file`==file, "EER_hits"] = length(eer_hits)

    results_DF[results_DF$`genome_file`==file, "RXLR_and_EER_hits"] = length(rxlr_eer_hits)
  } #inner for
  

  
  
write.table(results_DF, "20200529_unshuffled_String_strictRXLRlenientEER.csv", quote=FALSE, sep=",", row.names=FALSE)


