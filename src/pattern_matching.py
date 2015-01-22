import re
import sys

from Bio import SeqIO
import Bio.motifs as motifs
from Bio.Alphabet import IUPAC
import MOODS

from my_lib import searching_result_to_str


dna_alf = sorted(list(IUPAC.unambiguous_dna.letters))

with open(sys.argv[1]) as pwm_transfac_file:
    pwm_records = motifs.parse(pwm_transfac_file, "TRANSFAC")

with open(sys.argv[2]) as fasta_file:
    fasta_seqs = list(SeqIO.parse(fasta_file, 'fasta'))

tf_name = sys.argv[3]
pwm_matrix = [pwm for pwm in pwm_records if pwm['ID'] == tf_name][0]
consensus = pwm_matrix.consensus
pwm_matrix = pwm_matrix.counts
matrix = [pwm_matrix[n] for n in dna_alf]

output_file = open(sys.argv[3] + "_out.txt", 'w')

threshold = 0.7 * MOODS.max_score(matrix)
for i, fasta_seq in enumerate(fasta_seqs):
    sequence = str(fasta_seq.seq)
    interval = re.split("=| ", fasta_seq.description)[2]
    results = MOODS.search(sequence, [matrix], threshold, convert_log_odds=False, both_strands=True,
                           pseudocount=0, threshold_from_p=False)

    reversed_sequence = fasta_seq.seq[::-1]
    reversed_results = MOODS.search(reversed_sequence, [matrix], threshold, convert_log_odds=False, both_strands=True,
                                    pseudocount=0, threshold_from_p=False)

    sequence_length = len(sequence)
    result_str = searching_result_to_str(interval, results, reversed_results, sequence_length)
    output_file.write(str(result_str))

output_file.close()
