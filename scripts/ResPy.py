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
("So when you check the k-mers of the sample against the k-mers in the resistance database then the entire gene must have been covered"
'''

import gzip
import sys

class ResistanceFinder:
    def __init__(self, k=19):
        self.k = k
        self.kmer_db = {}  # kmer -> list of (gene_index, pos)
        self.genes = []    # List of (name, sequence, depth)
        
        # Create reverse complement table
        self.reverse_table = str.maketrans("ACGT", "TGCA")

    def load_fasta(self, fasta_path):
        """Parsing the fasta with the resistance genes and building the k-mer indexes"""
        print(f"Indexing FASTA: {fasta_path}")
        with open(fasta_path, 'r') as f:
            name, seq = None, []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if name:
                        self._process_gene(name, "".join(seq))
                    name = line[1:]
                    seq = []
                else:
                    seq.append(line)
            if name:
                self._process_gene(name, "".join(seq))

    def _process_gene(self, name, sequence):
	'''Extracting info about the resistance genes'''
        gene_idx = len(self.genes)
        sequence = sequence.upper()
        # Tracking depth: a number for each basepair, if a kmer is found in that position, +1 to that number
        depth = [0] * len(sequence)
        self.genes.append({'name': name, 'seq': sequence, 'depth': depth})

        # Slide window to find all kmer-s and put them in kmer_db
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i+self.k]
            if kmer not in self.kmer_db:
                self.kmer_db[kmer] = []
            self.kmer_db[kmer].append((gene_idx, i))

    def _get_rev_comp(self, seq):
        return seq.translate(self.trans_table)[::-1]

