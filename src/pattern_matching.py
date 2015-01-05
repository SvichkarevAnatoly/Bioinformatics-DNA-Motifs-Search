import sys

from Bio import SeqIO

from my_lib import *

with open(sys.argv[1]) as input_file:
    fasta_seqs = list(SeqIO.parse(input_file, 'fasta'))

genome = str(fasta_seqs[0].seq)
pattern = str(fasta_seqs[1].seq)

matching_list = pattern_matching_list(genome, pattern)
matching_list_str = matching_list_to_string(matching_list)

if len(sys.argv) == 3:
    with open(sys.argv[2], 'w') as output_file:
        output_file.write(matching_list_str)
else:
    print matching_list_str
