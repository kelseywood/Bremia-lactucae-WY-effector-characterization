## Running PONDR VSL2 on multiple genes

- (PONDR VSL2 documenation here: http://www.pondr.com)

- For this analysis, we needed to run PONDR VSL2 on many genomes. to do this, we wrote a bash script that effectively:
	- takes in a fasta file as input with multiple sequences
	- runs PONDR VSL2 jar file on the sequences one-by-one (requires VSL2.jar to be in the same directory)
	- outputs a txt file of the positional disorder scores, that can be parsed + graphed using R

- usage: 

	1. run VSL2 PONDR file to get the raw output into a textfile 
	
		example usage: bash run_pondr.sh (Input fasta file) (Output text file)
		
		```bash
		bash run_pondr.sh example.fasta example.out
		```
		
	2. read in the textfile with R and plot results
	
		there are example R files available in this directory 
		
		the most barebones one is **read_disorder.R**, which reads in example.out and plots a line graph of positional averages of all the sequences' disorder series, and limits it to the first 150 AA positions.
