library(effectR) ; library(seqinr)

## initialize DataFrame RESULTS_DF and filenames
	
  #use when working remotely
  #genomePath <- "/Users/munir/mnt/CLUSTER/effectR/motif_redo_201810/genomes/NR_cleaved_secretomes/10_shuffles/"
	
  
  #use when working with cluster
  setwd("/share/rwmwork/mjnur/effectR/motif_redo_201810/false_positive_analysis_201811")
  genomePath = "/share/rwmwork/mjnur/effectR/motif_redo_201810/genomes/NR_cleaved_secretomes/10_shuffles/"
	#end use when working with cluster
  
  shuffle_files <- list.files( paste0(genomePath, "filenames"), pattern="*.fasta")
  shuffle_files <- substring(shuffle_files,2)
	zeros = rep(0, length(shuffle_files))
	
	results_DF <- data.frame(shuffle_files, zeros, zeros,zeros, zeros, zeros, zeros, stringsAsFactors = FALSE)
	colnames(results_DF) = c("10_shuffle_file", "HMM_RXLR_hits", "HMM_WY_hits", "SD_HMM_RXLR", "SD_HMM_WY", "unshuffled_HMM_WYs", "unshuffled_HMM_RXLRs")


## start calculating the hits in each shuffle file
 	for (file in shuffle_files) {
 	  
 	  write(file, "progress.txt", append = TRUE)

 	  wy_results = c()
 	  rxlr_results = c()
 	  
 	  for (i in 1:10) {
 	    write(i, "progress.txt", append = TRUE)
 	    
 	    orf <- paste0(genomePath,"shuffle_",file,i)

	    rxlr_cmd = paste('hmmsearch --incT 0 --incdomT 0 --noali whisson_et_al_rxlr_eer_cropped.hmm', orf, " | grep 'over threshold' | awk ' {print $5} ' ")


    	    wy_cmd = paste('hmmsearch --incT 0 --incdomT 0 --noali WY_fold.hmm', orf, " | grep 'over threshold' | awk ' {print $5} ' ")

	   
 	    
 	    
 	    wy_hits = as.numeric( system(wy_cmd, intern=TRUE) )
	    rxlr_hits = as.numeric( system(rxlr_cmd, intern=TRUE) )

 	    rxlr_results = c(rxlr_results, rxlr_hits)
	    wy_results = c(wy_results, wy_hits)
 	  } #inner for
 	  
 	  results_DF[results_DF$`10_shuffle_file`==file, "HMM_RXLR_hits"] = sum(rxlr_results)
 	  results_DF[results_DF$`10_shuffle_file`==file, "SD_HMM_RXLR"] = sd(rxlr_results)
 	  
 	  results_DF[results_DF$`10_shuffle_file`==file, "HMM_WY_hits"] = sum(wy_results)
 	  results_DF[results_DF$`10_shuffle_file`==file, "SD_HMM_WY"] = sd(wy_results)
 	  
 	  
 	  
 	} #outer for
 	
 	write.table(results_DF, "FP_Results_HMMs_20190324.csv", quote=FALSE, sep=",", row.names=FALSE)
 	
 	
 	
