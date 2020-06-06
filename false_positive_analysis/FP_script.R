library(effectR) ; library(seqinr)

## initialize DataFrame RESULTS_DF and filenames
	
  #use when working remotely
  #genomePath <- "/Users/munir/mnt/CLUSTER/effectR/motif_redo_201810/genomes/NR_cleaved_secretomes/10_shuffles/"
	
  
  #use when working with cluster
  setwd("/share/rwmwork/mjnur/effectR/motif_redo_201810/false_positive_analysis_201811")
  genomePath = "/share/rwmwork/mjnur/effectR/motif_redo_201810/genomes/NR_cleaved_secretomes/10_shuffles/"
	#end use when working with cluster
  
  shuffle_files <- list.files( paste0(genomePath, "filenames"), pattern="*.fasta")
  shuffle_files <- substring(shuffle_files, 2)
	zeros = rep(0, length(shuffle_files))
	
	results_DF <- data.frame(shuffle_files, zeros, zeros,zeros,zeros, zeros, zeros, stringsAsFactors = FALSE)
	colnames(results_DF) = c("10_shuffle_file", "RXLR_hits", "EER_hits", "RXLR_and_EER_hits", "RXLR_SD", "EER_SD", "RXLR_and_EER_SD")

	RxLR <- "^\\w{1,56}R\\wLR"
	
	QGH_xLR <- "^\\w{1,56}[Q|G|H]\\wLR"
	
	RxL_QKG <- "^\\w{1,56}R\\wL[Q|K|G]"
	
	EER <- "^\\w{1,97}[D|E][D|E][K|R]"

## start calculating the hits in each shuffle file
 	for (file in shuffle_files) {
 	  
 	  write(file, "progress.txt", append = TRUE)

 	  rxlr_results = c()
 	  eer_results = c()
	  rxlr_and_eer_results = c()
 	  
 	  for (i in 1:10) {
 	    write(i, "progress.txt", append = TRUE)
 	    
 	    orf <- read.fasta(paste0(genomePath,"shuffle_",file,i))
 	    
 	    rxlr_hits = regex.search(orf, motif="custom", "^\\w{1,56}R\\wLR")
	    #rxlr_hits = regex.search(orf, motif="custom", "^\\w{1,56}[R|K|H|G|Q]\\w{1,2}[L|M|Y|F|I|V][R|N|K]")
 	    qghxlr_hits = regex.search(orf, motif="custom", "^\\w{1,56}[Q|G|H]\\wLR")
 	    #rxlqkg_hits = regex.search(orf, motif="custom", "^\\w{1,56}R\\wL[Q|K|G]")
 	    
 	    eer_hits = regex.search(orf, motif="custom", "^\\w{1,97}[D|E][D|E][K|R]")
 	    
	    rxlrs_hits = rxlr_hits
 	    rxlrs_hits = union(rxlr_hits, qghxlr_hits)
 	    #rxlrs_hits = union(rxlrs_hits, rxlqkg_hits)

	    rxlr_and_eer_hits = intersect(rxlrs_hits, eer_hits)
 	    
 	    rxlr_results = c(rxlr_results, length(rxlrs_hits))
 	    eer_results = c(eer_results, length(eer_hits))
	    rxlr_and_eer_results = c(rxlr_and_eer_results, length(rxlr_and_eer_hits))
 	  } #inner for
 	  
 	  results_DF[results_DF$`10_shuffle_file`==file, "RXLR_hits"] = sum(rxlr_results)
 	  results_DF[results_DF$`10_shuffle_file`==file, "RXLR_SD"] = sd(rxlr_results)
 	  
 	  results_DF[results_DF$`10_shuffle_file`==file, "EER_hits"] = sum(eer_results)
 	  results_DF[results_DF$`10_shuffle_file`==file, "EER_SD"] = sd(eer_results)

	  results_DF[results_DF$`10_shuffle_file`==file, "RXLR_and_EER_hits"] = sum(rxlr_and_eer_results)
	  results_DF[results_DF$`10_shuffle_file`==file, "RXLR_and_EER_SD"] = sd(rxlr_and_eer_results)

 	  
 	  
 	  
 	} #outer for
 	
 	write.table(results_DF, "FP_Results_20200513_kylemotifs.csv", quote=FALSE, sep=",", row.names=FALSE)
 	
 	
 	
