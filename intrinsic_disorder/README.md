## Running PONDR VSL2 on multiple genes

- (PONDR VSL2 documenation here: http://www.pondr.com)

- For this analysis, we needed to run PONDR VSL2 on many genomes. to do this, we wrote a bash script that effectively:
	- takes in a fasta file as input with multiple sequences
	- runs PONDR VSL2 on the sequences one-by-one
	- outputs a txt file of the positional disorder scores, that can be parsed + graphed using R

- usage: 
