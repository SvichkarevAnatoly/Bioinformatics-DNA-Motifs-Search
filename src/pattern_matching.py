import sys

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
print pwm_matrix.consensus
pwm_matrix = pwm_matrix.counts
matrix = [pwm_matrix[n] for n in dna_alf]

for fasta_seq in fasta_seqs:
    sequence = str(fasta_seq.seq)
    print sequence
    results = MOODS.search(sequence, [matrix], 1500, convert_log_odds=False,
                           threshold_from_p=False)
    print results
    break

print pwm_matrix.consensus
print sequence[56:56+len(matrix[0])]

# if len(sys.argv) == 3:
#     with open(sys.argv[2], 'w') as output_file:
#         output_file.write(matching_list_str)
# else:
#     print matching_list_str
