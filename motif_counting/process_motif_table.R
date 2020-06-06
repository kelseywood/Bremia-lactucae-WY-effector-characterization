#DF <- read.csv("/share/rwmwork/mjnur/effectR/
#motif_redo_201810/results_20181015/EffectR_results_NR_cleaved_secretomes.txt", stringsAsFactors = F)
#setwd("~/mnt/CLUSTER/effectR/motif_redo_201810/results_20190123/")
setwd("/share/rwmwork/mjnur/effectR/motif_redo_201810/results_20190123")
DF <- read.csv("EffectR_results_NR_cleaved.txt", stringsAsFactors = F)

#import WY and HMM_RXLR
WYs <- read.table("all_WY_cleaved_IDs.csv")
hmm_RXLRs <- read.table("all_whissonHMM_IDs.csv")

WY <- rep('N', nrow(DF))
hmm_RXLR <- rep('N', nrow(DF))

DF <- data.frame(DF, WY, stringsAsFactors = F)

for (gene in WYs$V1) {
  DF[DF$seqID==gene, "WY"] <- "Y"
}
for (gene in hmm_RXLRs$V1) {
  DF[DF$seqID==gene, "hmm_RXLR"] <- "Y"
}

#end importing WY and HMM_RXLR information

filenames <- as.character(unique(unlist(DF$filename)))

#DF_RxLR contains RxLR, GxLR, QxLR, RxLK (and all 4 HMM) results

DF_RxLR <- DF[which(DF$RxLR=='Y' | DF$QGH_xLR=='Y' | DF$hmm_RXLR=='Y'), ]

DF_RxLR_EER <- DF[which(DF$RxLR_EER=='Y' | DF$QGH_xLR_EER=='Y'), ]

DF_EER <- DF[which(DF$EER=='Y'), ]

#DF_WY contains WY results
DF_WY <- DF[which(DF$WY=='Y'), ]

zeros <- rep(0, length(filenames))
resultsDF <- data.frame(filenames, zeros, zeros, zeros, zeros, zeros, zeros, zeros,
                        zeros, zeros, zeros, zeros,zeros,zeros,zeros,zeros,zeros, zeros, zeros)
                        
names(resultsDF) <- c("filename", "num_NR_secreted_proteins", 
                      "RXLR", "WY", "EER", #singles
                      "RXLREER", "RXLR_and_EER","RXLR_EER_noWY","RXLR_noEER", 
                      "RXLR_noEER_noWY",  "noRXLR_EER_noWY",
                      "WY_noRXLR", "WY_and_RXLR", #doubles, EER after and anywhere
                      "WY_and_RXLR_and_EER", "WY_and_RXLR_noEER","WY_noRXLR_and_EER", "WY_noRXLR_noEER",
                      "WY_and_RXLREER", "WY_noRXLREER")



for (genome in filenames) {
  resultsDF[resultsDF$filename==genome, ]$num_NR_secreted_proteins <- 
    nrow(DF[which(DF$filename==genome), ])
  
  genome_DF <- DF[which(DF$filename==genome), ]
  genome_DF_RxLR <- DF_RxLR[which(DF_RxLR$filename==genome), ]
  genome_DF_EER <- DF_EER[which(DF_EER$filename==genome), ]
  genome_DF_WY <- DF_WY[which(DF_WY$filename==genome), ]
  genome_DF_RxLREER <- DF_RxLR_EER[which(DF_RxLR_EER$filename==genome), ]
  
  
  #singles
  resultsDF[resultsDF$filename==genome, ]$RXLR <- nrow(genome_DF_RxLR)
  resultsDF[resultsDF$filename==genome, ]$EER <- nrow(genome_DF_EER)
  resultsDF[resultsDF$filename==genome, ]$WY <- nrow(genome_DF_WY)
  resultsDF[resultsDF$filename==genome, ]$RXLREER <- nrow(genome_DF_RxLREER)
  
  #doubles
  #setdiff(x,y) returns items in x that are NOT in y
  
  resultsDF[resultsDF$filename==genome, ]$RXLR_and_EER <- 
    nrow(genome_DF_RxLR[which(genome_DF_RxLR$EER =='Y'), ])

  resultsDF[resultsDF$filename==genome, ]$RXLR_noEER <- 
    length(setdiff(genome_DF_RxLR$seqID, genome_DF_EER$seqID))

  resultsDF[resultsDF$filename==genome, ]$WY_and_RXLREER <- 
    nrow(genome_DF_RxLREER[which(genome_DF_RxLREER$WY =='Y'), ])
  
  resultsDF[resultsDF$filename==genome, ]$WY_and_RXLR <- 
    nrow(genome_DF_RxLR[which(genome_DF_RxLR$WY =='Y'), ])
  
  resultsDF[resultsDF$filename==genome, ]$WY_noRXLR <- 
    length(setdiff(genome_DF_WY$seqID, genome_DF_RxLR$seqID))
  
  resultsDF[resultsDF$filename==genome, ]$WY_noRXLREER <-
    length(setdiff(genome_DF_WY$seqID, genome_DF_RxLREER$seqID))
  
  #triples
  #"WY_and_RXLR_and_EER", "WY_and_RXLR_noEER","WY_noRXLR_and_EER", "WY_noRXLR_noEER",
  WYnoRXLR = setdiff(genome_DF_WY$seqID, genome_DF_RxLR$seqID)
  WYandRXLR = intersect(genome_DF_WY$seqID, genome_DF_RxLR$seqID)
  
  #WYnoRXLR:
  resultsDF[resultsDF$filename==genome, ]$WY_noRXLR_noEER <-
    length(setdiff(WYnoRXLR, genome_DF_EER$seqID))
  
   resultsDF[resultsDF$filename==genome, ]$WY_noRXLR_and_EER <-
     length(intersect(WYnoRXLR, genome_DF_EER$seqID))
  
  #WYandRXLR:
  resultsDF[resultsDF$filename==genome, ]$WY_and_RXLR_noEER <-
    length(setdiff(WYandRXLR, genome_DF_EER$seqID))
  
  resultsDF[resultsDF$filename==genome, ]$WY_and_RXLR_and_EER <-
    length(intersect(WYandRXLR, genome_DF_EER$seqID))
    
 # "RXLR_noEER_noWY", "RXLR_EER_noWY", "noRXLR_EER_noWY",
  # setdiff(x,y) returns items in x that are NOT in y
  resultsDF[resultsDF$filename==genome, ]$RXLR_noEER_noWY <-
    nrow(genome_DF_RxLR[which(genome_DF_RxLR$WY=='N' & genome_DF_RxLR$EER=='N'), ])
  
    
  resultsDF[resultsDF$filename==genome, ]$RXLR_EER_noWY <-
    nrow(genome_DF_RxLR[which(genome_DF_RxLR$WY=='N' & genome_DF_RxLR$EER=='Y'), ])

    
  EERnoWY = genome_DF_EER[which(genome_DF_EER$WY=='N'), ]
  resultsDF[resultsDF$filename==genome, ]$noRXLR_EER_noWY <-
    length(setdiff(EERnoWY$seqID, genome_DF_RxLR$seqID))
}

write.csv(resultsDF, file="2020603_Overlap_Table_strictRXLRnoHMM.csv", quote=F, row.names=F)




#upset plots
library(UpSetR)
for (genome in filenames) {
  
  genome_DF <- DF[which(DF$filename==genome), ]
  genome_DF_RxLR <- DF_RxLR[which(DF_RxLR$filename==genome), ]
  genome_DF_EER <- DF_EER[which(DF_EER$filename==genome), ]
  genome_DF_WY <- DF_WY[which(DF_WY$filename==genome), ]
  genome_DF_RxLREER <- DF_RxLR_EER[which(DF_RxLR_EER$filename==genome), ]
  
  RxLR_IDs <- genome_DF_RxLR$seqID
  EER_IDs <- genome_DF_EER$seqID
  RxLR_EER_IDs <- genome_DF_RxLREER$seqID
  WY_IDs <- genome_DF_WY$seqID
  
  RxLR_and_EER <- intersect(RxLR_IDs, EER_IDs)
  
  listInput <- list(WY = WY_IDs, RXLR = RxLR_IDs, EER = EER_IDs)
  myIntersections <- list(list("WY", "RXLR"), list("WY"), list("WY","RXLR", "EER"))
  #pdf(file=paste0("UPSET_", genome, ".pdf"))
  upset(fromList(listInput), order.by="freq", 
        sets.bar.color = "#56B4E9", intersections = myIntersections)
  #dev.off()
}


# load libraries
require(ggplot2)
require(reshape2)
library(tidyverse)

# melt data from wide to long
copy <- resultsDF
copy = copy %>% select(filename, WY_noRXLR_noEER, 
                       WY_noRXLR_and_EER, WY_and_RXLR_noEER, 
                       WY_and_RXLR_and_EER)

dat_m <- melt(copy, id.vars = "filename")

# plot

ggplot(dat_m, aes(x=filename, y=value, fill = variable)) + 
  geom_bar(stat = "identity")
