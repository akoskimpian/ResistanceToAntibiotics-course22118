'''-Pseudocode-

Import libraries

Read in resistance fasta

Resistance K-mer dict = {key: kmer seq, value: genes it is in} somehow store number of unique kmers for every gene, in order to calculate coverage

Create empty results dict: store how many times we have seen a res kmer: important because need a depth around 50 to validate finding 
Parse through fastq-s, get sequence Fastq-s are big -> a way to quickly load them and deal with them? 
Generate reverse complement of sequence 
For kmer in fastq: 
	if kmer in res_kmer_dict: 
		add kmer to results dict, value is count, +1 to count

If kmer was seen less than 10 times, delete it Calculate the coverage of each resistance gene, 
how many unique kmers of it were found / total unique kmers of it. We need 100% if i understand correctly 
("So when you check the k-mers of the sample against the k-mers in the resistance database then the entire gene must have been covered"'''
