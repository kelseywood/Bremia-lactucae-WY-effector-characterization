get_averages <- function(PONDR_results_file) {
  #function to parse PONDR results file, get averages for each AA position
  
  #print(PONDR_results_file)
  
  WY_results = read.table(PONDR_results_file, comment.char="-")
  names(WY_results) = c("Position", "AA", "Score", "Disorder")
  #names(WY_results) = c("Position", "AA", "Score")
  
  
  maxPos <- max(WY_results$Position)
  averages <- c(rep( -1, maxPos))
  
  for (pos in 1:maxPos) {
    #get all values at location pos for each AA sequence
    pos_items <- subset(WY_results, Position == pos)
    #get means
    averages[pos] <- mean(pos_items$Score)
  }
  
  return(averages)
} #returns positional averages array for a PONDR results file


plotter <- function(rawTitle, prettyTitle) {
 
  #put in beginning of plotter:
  #read in genomes.txt
  #read in old names.txt
  #RXLR_path equals ls(RXLR_resPath).str.contains(RXLR_and_EER)
  #RXLR_path equals ls(RXLR_resPath).str.contains(RXLR_and_EER)
  RXLR_EERs_file = list.files(path = "/share/rwmwork/mjnur/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811/Results_20181102/Results_20181102_RXLR_and_EER",
                                 pattern = paste0('*', rawTitle, '*'), include.dirs=TRUE)
  RXLR_EERs_file = paste0("/share/rwmwork/mjnur/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811/Results_20181102/Results_20181102_RXLR_and_EER/", 
                          RXLR_EERs_file[1])
  
  WYs_withRXLR_file = list.files(path = "/share/rwmwork/mjnur/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811/Results_20181102/Results_20181102_WY_and_RXLR",
                                 pattern = paste0('*', rawTitle, '*'), include.dirs=TRUE)
  WYs_withRXLR_file = paste0("/share/rwmwork/mjnur/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811/Results_20181102/Results_20181102_WY_and_RXLR/",
                           WYs_withRXLR_file[1])
  
  WYs_noRXLR_file = list.files(path = "/share/rwmwork/mjnur/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811/Results_20181102/Results_20181102_WY_noRXLR",
                                 pattern = paste0('*', rawTitle, '*'), include.dirs=TRUE)
  WYs_noRXLR_file = paste0("/share/rwmwork/mjnur/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811/Results_20181102/Results_20181102_WY_noRXLR/",
                            WYs_noRXLR_file[1])
  
  secretomes_file = list.files(path = "/share/rwmwork/mjnur/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811/Results_20181102/Results_20181102_Secretomes",
                                 pattern = paste0('*', rawTitle, '*'), include.dirs=TRUE)
  secretomes_file = paste0("/share/rwmwork/mjnur/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811/Results_20181102/Results_20181102_Secretomes/",
                           secretomes_file[1])
  
  WYs_noRXLR = get_averages((WYs_noRXLR_file))
  WYs_withRXLR = get_averages((WYs_withRXLR_file))
  secretome = get_averages(secretomes_file)
  #secretome = get_averages((WYs_noRXLR_file))
  RXLR_EERS = get_averages((RXLR_EERs_file))
  
  
  
  
  pos = seq(1:150)
  disorder = data.frame(WYs_noRXLR[1:150], 
                      RXLR_EERS[1:150], WYs_withRXLR[1:150], secretome[1:150], pos)
  colnames(disorder) = c("WYnoRXLR", "RXLR_EERs", "WYwithRXLR","Secretome", "pos")
  disorder2 = melt(disorder, id.var="pos")

  title = "genome"

  ##save as pdf 6 by 8 inches :) meow nom nom nom
  x1 = ggplot(data=disorder, aes(x=pos)) +
  geom_line(aes(y=WYnoRXLR, colour = "WY (no RXLR)"), size=1.2) + 
  geom_line(aes(y=WYwithRXLR, colour = "WY (with RXLR)"), size=1.2) + 
  geom_line(aes(y=RXLR_EERs, colour = "RXLR+EER"), size=1.2) +
  geom_line(aes(y=Secretome, colour = "Secretome"), size=1.2) +
  ggtitle(label=prettyTitle) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) + 
  #geom_hline(yintercept=0.5, linetype="dashed", 
  #           color = "black", size=1) + 
  ylab(label=NULL) +
  xlab(label=NULL) +
  theme(plot.title = element_text(color="black", face="italic"),
	legend.title=element_blank()) +
   # theme(legend.position = "none",
          #plot.title = element_text(color="black", size=21, face="bold"),
          #axis.title.x = element_text(color="black", size=21),
          #axis.title.y = element_text(color="black", size=21),
          #axis.text.y = element_text(size=14),
          #axis.text.x = element_text(size=14)) +
  scale_colour_manual(breaks = c("WY (no RXLR)", "WY (with RXLR)", "RXLR+EER", "Secretome"),
                      values = c("#7030a0", "#4372c4", "#fec002", "#444444")) 
  return(x1)
}



library("reshape2")
library("ggplot2")
library("ggpubr")
library("grid")
library("gridExtra")

#setwd("/Users/munir/mnt/CLUSTER/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811")

setwd("/share/rwmwork/mjnur/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811")
rawNames = readLines("rawGenomes.txt")
prettyNames = readLines("genomes.txt")

disorderGraphs = c()
#for (genome in 1:length(rawNames)) {
for (genome in 1:length(rawNames)) {
  cat(rawNames[genome], " = ", prettyNames[genome], "\n")
  newPlot = plotter(rawNames[genome], prettyNames[genome])
  disorderGraphs[[genome]] = newPlot
}

#put in beginning of plotter:
  #read in genomes.txt
  #read in old names.txt
  #RXLR_path equals ls(RXLR_resPath).str.contains(RXLR_and_EER)
  #RXLR_path equals ls(RXLR_resPath).str.contains(RXLR_and_EER)




#bremia = get_averages("/Users/munir/mnt/CLUSTER/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811/Results_20181102/Results_20181102_WY_noRXLR/PONDR_WYnoRXLR_NR_secreted_B_lac.protein.fasta.txt")


ggarrange(disorderGraphs[[2]], disorderGraphs[[3]], disorderGraphs[[4]], 
          disorderGraphs[[5]], disorderGraphs[[6]], disorderGraphs[[7]], disorderGraphs[[8]],
          disorderGraphs[[9]], disorderGraphs[[10]], disorderGraphs[[11]], disorderGraphs[[12]], disorderGraphs[[13]],
          ncol=4, nrow=3, common.legend=TRUE, legend="bottom")

ggsave("20190403_rest_of_genomes_legendRight_4by3.pdf", plot = last_plot(), device = "pdf",
       scale = 1, width = 8, height = 6, units = "in",
       dpi = 300, limitsize = TRUE)

