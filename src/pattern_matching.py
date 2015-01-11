import sys

from my_lib import searching_result_to_str

from Bio import SeqIO
import Bio.motifs as motifs
from Bio.Alphabet import IUPAC
import MOODS

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

threshold = 0.7 * MOODS.max_score(matrix)
for i, fasta_seq in enumerate(fasta_seqs):
    sequence = str(fasta_seq.seq)
    results = MOODS.search(sequence, [matrix], threshold, convert_log_odds=False,
                           pseudocount=0, threshold_from_p=False)
    print searching_result_to_str(sequence, consensus, results, MOODS.max_score(matrix))

# if len(sys.argv) == 3:
# with open(sys.argv[2], 'w') as output_file:
# output_file.write(matching_list_str)
# else:
#     print matching_list_str
