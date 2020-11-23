get_averages <- function(PONDR_results_file) {
  #function to parse PONDR results file, get averages for each AA position
  
  WY_results = read.table(PONDR_results_file, comment.char="-")
  names(WY_results) = c("Position", "AA", "Score", "Disorder")

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
  
  #return(list(averages[1:150], mins[1:150], maxs[1:150]))
  return(averages)
} #returns positional averages array for a PONDR results file

example_output <- get_averages("example.out")


library("reshape2")
library("ggplot2")

# make disorder dataframe with file1's output
# [1:150] to force the list to be 150 sequences long, even if the end ones are empty
disorder = data.frame(example_output[1:150], seq(1:150))
colnames(disorder) = c("gene_file1", "pos")

# plot data
ggplot(data=disorder, aes(x=seq(1:150))) +
  geom_line(aes(x=pos, y=gene_file1, colour = "gene_file1"), size=1.8)

# example to plot 2 groups of data:
   # example_output <- get_averages("example.out")
   # example_output2 <- get_averages("example2.out")
   # disorder = data.frame(example_output[1:150], example_output2[1:150], seq(1:150))
   # colnames(disorder) = c("gene_file1", "gene_file2", "pos")
   # disorder2 = melt(disorder, id.var="pos")
   # 
   # ggplot(data=disorder, aes(x=seq(1:150))) +
   #    geom_line(aes(x=pos, y=gene_file1, colour = "gene_file1"), size=1.8) +
   #    geom_line(aes(x=pos, y=gene_file2, colour = "gene_file2"), size=1.8)


# saving plot:
# ggsave("20190324_bremiaDisorder_newColors.pdf", plot=last_plot(), width=8, height=6, units="in")