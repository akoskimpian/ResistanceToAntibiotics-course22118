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
from array import array


class ResistanceFinder:
    def __init__(self, k=19):
        self.k = k
        self.kmer_db = {}  # kmer -> list of (gene_idx, pos)
        self.genes = []    # List of dictionaries
        self.trans_table = str.maketrans("ACGT", "TGCA")


    def load_fasta(self, fasta_path):
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
        gene_idx = len(self.genes)
        sequence = sequence.upper()
        # Difference array: length + 1 to handle the boundary of 'pos + k'
        # Using 'i' for signed integers for negative values
        diff_arr = array('i', [0] * (len(sequence) + 1))
        self.genes.append({'name': name, 'seq': sequence, 'diff': diff_arr})


        # Slide window to populate kmer_db
        db = self.kmer_db
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i+self.k]
            if kmer not in db:
                db[kmer] = []
            db[kmer].append((gene_idx, i))


    def _get_rev_comp(self, seq):
        return seq.translate(self.trans_table)[::-1]


    def process_fastq(self, fastq_path):
        print(f"Processing FASTQ: {fastq_path}")
        # Localize k and map_read function for speed
        k = self.k
        db = self.kmer_db
        genes = self.genes
       
        with gzip.open(fastq_path, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    read = line.strip().upper()
                    # Inline mapping - new step
                    for sequence_variant in (read, self._get_rev_comp(read)):
                        for j in range(len(sequence_variant) - k + 1):
                            kmer = sequence_variant[j:j+k]
                            if kmer in db:
                                for gene_idx, pos in db[kmer]:
                                    d = genes[gene_idx]['diff']
                                    d[pos] += 1
                                    d[pos + k] -= 1


    def report(self, min_coverage=0.95, min_depth=10):
        print(f"\n{'Gene Name':<50} | {'Coverage':<10} | {'Avg Depth':<10}")
        print("-" * 75)
       
        results = []
        for gene in self.genes:
            diff = gene['diff']
            length = len(gene['seq'])
           
            # Reconstruct actual depth from difference array - prefi sum method, new step
            actual_depth = [0] * length
            current_val = 0
            covered_bases = 0
            total_sum = 0
           
            for i in range(length):
                current_val += diff[i]
                actual_depth[i] = current_val
                if current_val > 0:
                    covered_bases += 1
                    total_sum += current_val
           
            coverage = covered_bases / length
            avg_depth = total_sum / length
           
            if coverage >= min_coverage and avg_depth >= min_depth:
                results.append({
                    'name': gene['name'],
                    'coverage': coverage,
                    'depth': avg_depth
                })


        results.sort(key=lambda x: (x['coverage'], x['depth']), reverse=True)
        for r in results:
            print(f"{r['name'][:50]:<50} | {r['coverage']:>9.2%} | {r['depth']:>10.2f}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        # Defaults for your specific files
        default_fasta = "resistance_genes.fsa.txt"
        default_fastq = "Unknown3_raw_reads_1.txt.gz"
       
        scanner = ResistanceFinder(k=19)
        scanner.load_fasta(default_fasta)
        scanner.process_fastq(default_fastq)
        scanner.report()
    else:
        scanner = ResistanceFinder(k=19)
        scanner.load_fasta(sys.argv[1])
        scanner.process_fastq(sys.argv[2])
        scanner.report()


