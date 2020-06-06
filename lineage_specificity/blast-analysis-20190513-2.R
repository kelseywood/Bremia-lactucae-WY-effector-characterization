# the Bremia WYs were blasted against the 13 other genomes using blastp
# blastp results file  uses default e-value cutoff blastp used, 
#    but in R I only use hits with e-value < 0.01 (line 58)

setwd("/Users/munir/Desktop/20190531_blast_test")

#### get prefixes of fasta sequences (of database to blast against) to be able to know which genome a blastp hit came from
# before blasting, each genome you're blasting against must have a genome-wide prefix for gene IDs.
# 2 methods of doing this in the command line (inplace and creating new file)"

## (in command line) add prefix to all gene IDs inplace to fasta file:
# perl -pi -e "s/^>/>phosphate-/g" your.fasta

## (in command line) alternatively, can do this to generate a new fasta file with prefixes added:
#     perl -p -e "s/^>/>phosphate-/g" your.fasta > with_prefixes.fasta

#### edit filename_to_prefix_to_genome.txt to include prefixes used to subset blast results
# column 1: prefix of the genome
# column 2: genome fasta file
# column 3: output genome title you want in plot
prefixes <- as.character(read.csv(pipe("cut -f1 -d',' filename_to_prefix_to_genomeName.csv"), header=F)$V1)
genome_names <- as.character(read.csv(pipe("cut -f3- -d',' filename_to_prefix_to_genomeName.csv"), header=F)$V1)

#### read in data produced by blastp (make sure to use -out "6 std qcovs for blast run)
# data was produced with:
#     blastp -db DB_July2018 -query secreted_wy_proteins.fasta -outfmt "6 std qcovs" -out results_file.tab 
data <- read.table("20190510_secretedWYs_vs_all_addedGenomes.tab",sep="\t",stringsAsFactors =F)

# rename columns to blast -outfmt 6 std qcovs standard
names(data) <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

genomes <- list()
uniq_genomes <- list()

for (i in 1:length(prefixes)) {
  new_frame <- data[grep(prefixes[i], data$sseqid),]
  
  cat("processing: ", prefixes[i], " ( # of hits: ", nrow(new_frame) , ")\n")
  
  # sort by bitscore, arrange by query sequence in bremia (for viewing purposes)
  sorted_by_bitscore <- new_frame[ order(new_frame[,'qseqid'], new_frame[,'bitscore'], decreasing=T), ]
  
  # only keep highest bitscore hit for every query sequence
  # i think this should keep the highest quality of alignment, before looking at % identity
  best_hit <- sorted_by_bitscore[ !duplicated(sorted_by_bitscore$qseqid), ]
  
  # add species name so we know where these hits came from
  best_hit$species <- genome_names[i]

  # append the best hit for each WY in bremia. so in total, only as many hits as #WY sequences
  genomes[[i]] <- best_hit
}

# genomes[[]] is a list of dataframes, so bind them all into one dataframe
# we have created a column for species above in the for loop, so we know what where hits came from
all_oomycetes <- do.call("rbind", genomes)

# select hits where the e-value is at least smaller than 0.01 to get rid of garbage hits
selected_data <- subset(all_oomycetes, evalue<0.01)


# now time to plot:
library("ggplot2")
library("gtools")


ggplot(selected_data, aes(factor(species, levels=rev(unique(species))), pident),trim=FALSE) +
  ## violin or box plot:
  
  #geom_violin(aes(fill=factor(species, levels=(unique(species))))) + 
  #geom_boxplot(aes(fill=factor(species, levels=(unique(species))))) + 
  geom_boxplot(outlier.shape=NA, aes(fill=factor(species, levels=(unique(species))))) +

  ## include all individual points:
  geom_dotplot(binaxis='y', dotsize=0.45, binwidth=1, stackdir="centerwhole",
               method="histodot") +

  theme(axis.text.y = element_text(angle=0, hjust=1, face="italic", size=12)) +

  theme(legend.position = "none") + # removes legend
  coord_flip(ylim=c(18.5,100.5)) + 
  ylab("Best BLAST hit % identity") + 
  xlab("Species") + 
  theme(axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14))

#### save as PDF 6 x 8 inches
# ggsave("20190510_plots/20190513_secreted_WYs_vs_all_removed_empty.pdf", plot=last_plot(), width=8, height=6, units="in")
