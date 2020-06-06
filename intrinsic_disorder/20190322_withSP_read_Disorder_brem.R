get_averages <- function(PONDR_results_file) {
  #function to parse PONDR results file, get averages for each AA position
  
  WY_results = read.table(PONDR_results_file, comment.char="-")
  names(WY_results) = c("Position", "AA", "Score", "Disorder")
  #names(WY_results) = c("Position", "AA", "Score")
  
  
  maxPos <- max(WY_results$Position)
  averages <- c(rep( -1, maxPos))
  mins <- c(rep( -1, maxPos))
  maxs <- c(rep( -1, maxPos))
  
  
  for (pos in 1:maxPos) {
    #get all values at location pos for each AA sequence
    pos_items <- subset(WY_results, Position == pos)
    #get means
    averages[pos] <- mean(pos_items$Score)
    std <- function(x) sd(x) / sqrt(length(x))
    
    mins[pos] <- averages[pos] - std(pos_items$Score)
    maxs[pos] <- averages[pos] + std(pos_items$Score)
    
  }
  
  return(list(averages[1:150], mins[1:150], maxs[1:150]))
} #returns positional averages array for a PONDR results file

#make sure to highlight and run the above code so get_averages() function is saved in environment

#Blac_WYs <- get_averages("Results_pondr_Blac_cleavedWYs_Jul29.txt")
setwd("/Users/munir/mnt/CLUSTER/effectR/motif_redo_201810/intrinsic_Disorder_analysis_201811")
Blac_WYs_withRXLR <- get_averages("Results_20181102/Results_20181102_WY_and_RXLR/PONDR_WY_and_RXLR_NR_secreted_B_lac.protein.fasta.txt")
Blac_WYs_noRXLR <- get_averages("Results_20181102/Results_20181102_WY_noRXLR/PONDR_WYnoRXLR_NR_secreted_B_lac.protein.fasta.txt")
Blac_RXLR_EERS <- get_averages("Results_20181102/Results_20181102_RXLR_and_EER/PONDR_RXLR_and_EER_NR_secreted_B_lac.protein.fasta.txt")
secretome <- get_averages("Results_20181102/Results_20181102_Secretomes/PONDR_Secretome_NR_secreted_B_lac.protein.fasta.txt")
#secretome <- get_averages("Results_20181102/Results_20181102_RXLR_and_EER/PONDR_RXLR_and_EER_NR_secreted_B_lac.protein.fasta.txt")


library("reshape2")
library("ggplot2")

pos = seq(1:150)
disorder = data.frame(Blac_WYs_noRXLR[[1]], Blac_RXLR_EERS[[1]], Blac_WYs_withRXLR[[1]], secretome[[1]], pos)
colnames(disorder) = c("WYnoRXLR", "RXLR_EERs", "WYwithRXLR","Secretome", "pos")
disorder2 = melt(disorder, id.var="pos")



##save as pdf 6 by 8 inches :) meow nom nom nom
ggplot(data=disorder, aes(x=seq(1:150))) +
  geom_line(aes(x=pos, y=WYnoRXLR, colour = "WY (no RXLR)"), size=1.8) +
  geom_ribbon(alpha=0.3, aes(ymin=Blac_WYs_noRXLR[[2]], ymax=Blac_WYs_noRXLR[[3]]), fill="#fec002" ) +
  
  geom_line(aes(x=pos, y=WYwithRXLR, colour = "WY (with RXLR)"), size=1.8) +
  geom_ribbon(alpha=0.3, aes(ymin=Blac_WYs_withRXLR[[2]], ymax=Blac_WYs_withRXLR[[3]]), fill="#444444" ) +
  
  geom_line(aes(x=pos, y=RXLR_EERs, colour = "RXLR+EER"), size=1.8) +
  geom_ribbon(alpha=0.3, aes(ymin=Blac_RXLR_EERS[[2]], ymax=Blac_RXLR_EERS[[3]]), fill="#7030a0" ) +
  
  geom_line(aes(x=pos, y=Secretome, colour = "Secretome"), size=1.8) +
  geom_ribbon(alpha=0.3, aes(ymin=secretome[[2]], ymax=secretome[[3]]), fill="#4372c4" ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 

  ylab(label="VSL2 disorder score") +
  xlab(label="Amino acid position") +
  theme(legend.position = c(1, 1), legend.justification = c(1,1), 
        legend.title=element_blank(), legend.text=element_text(size=13),
        legend.key.height = unit(.79, "cm"), legend.key.width = unit(.48, "cm"), 
        plot.title = element_text(color="black", size=21, face="bold"),
        axis.title.x = element_text(color="black", size=21),
        axis.title.y = element_text(color="black", size=21),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14)) +
  scale_colour_manual(breaks = c("WY (no RXLR)", "WY (with RXLR)", "RXLR+EER", "Secretome"),
                      values = c("#7030a0", "#4372c4", "#fec002", "#444444")) 
#rxlr+eer, secretome, wynorxlr, wywithrxlr
#values = c("#7030a0", "#009E73", "#666666", "#fec002")) 

 # scale_colour_manual(breaks = c("WYs (no RXLR)", "WYs (with RXLR)", "RXLR + EER", "Secretome"),
  #                    values = c("#009E73", "#666666", "#D55E00", "#0072B2")) 

ggsave("20190324_bremiaDisorder_newColors.pdf", plot=last_plot(), width=8, height=6, units="in")
