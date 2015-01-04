import sys

from Bio import SeqIO
from my_lib import pattern_matching_list

with open(sys.argv[1]) as input_file:
    fasta_seqs = list(SeqIO.parse(input_file, 'fasta'))

genome = str(fasta_seqs[0].seq)
pattern = str(fasta_seqs[1].seq)

matching_list = pattern_matching_list(genome, pattern)

if len(sys.argv) == 3:
    with open(sys.argv[2], 'w') as output_file:
        output_file.write(' '.join(map(str, matching_list)) + '\n')
else:
    print ' '.join(map(str, matching_list))
