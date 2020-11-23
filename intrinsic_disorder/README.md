## Running PONDR VSL2 on multiple genes

- (PONDR VSL2 documenation here: http://www.pondr.com)

- For this analysis, we needed to run PONDR VSL2 on many genomes. to do this, we wrote a bash script that effectively:
	- takes in a fasta file as input with multiple sequences
	- runs PONDR VSL2 jar file on the sequences one-by-one (requires VSL2.jar to be in the same directory)
	- outputs a txt file of the positional disorder scores, that can be parsed + graphed using R

- usage: 

	1. run VSL2 PONDR file to get the raw output into a textfile. 
	
	**Note**: In the background, this **run_pondr.sh** script reformats the fasta file to be a single line fasta file (each sequence string is in a single line and doesn't span multiple lines), and the script also removes the '*' characters. these are requirements by VSL2, and you can see the resulting formatted fasta file in the same directory with the '_formatted' suffix if you want to.
	
	example usage: bash run_pondr.sh (Input fasta file) (Output text file)
		
		bash run_pondr.sh example.fasta example.out
		
	2. read+parse the textfile with R and plot results
	
	there are example R files available in this directory, the most barebones one is **read_disorder.R**, which reads+parses in example.out and plots a line graph of positional averages of all sequences' disorder series, and limits it to the first 150 AA positions. you can run this in RStudio, line by line and see the data turn into a nice figure!
