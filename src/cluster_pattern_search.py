import os
import re

from Bio import SeqIO
import Bio.motifs as motifs
from Bio.Alphabet import IUPAC
import MOODS


dna_alf = sorted(list(IUPAC.unambiguous_dna.letters))

pwms_file_name = "../data/PWMs_TRANSFAC.txt"
cluster_fasta_file_name = "../data/clusters/mm10_clusters.fa"

with open(pwms_file_name) as pwm_transfac_file:
    pwm_records = motifs.parse(pwm_transfac_file, "TRANSFAC")

with open(cluster_fasta_file_name) as cluster_fasta_file:
    cluster_seqs = list(SeqIO.parse(cluster_fasta_file, 'fasta'))

matrix_id_list = []
matrix_list = []
threshold_list = []
for pwm in pwm_records:
    matrix_id_list.append(pwm['ID'])
    matrix = [pwm.counts[n] for n in dna_alf]
    threshold = 0.7 * MOODS.max_score(matrix)
    threshold_list.append(threshold)
    matrix_list.append(matrix)

cluster_file_path, cluster_file_basename = os.path.split(cluster_fasta_file_name)
cluster_file_basename_list = os.path.splitext(cluster_file_basename)
cluster_match_file_name = cluster_file_basename_list[0] + '_match.txt'
cluster_match_full_file_name = os.path.join(cluster_file_path, cluster_match_file_name)


def result_to_str(interval, result_list, sequence_length):
    result_str = ">" + interval + '\n'
    for ind, result in enumerate(result_list):
        pos_str = matrix_id_list[ind] + ' '
        for (pos, score) in result:
            if pos < 0:
                pos_str += '-'
                pos += sequence_length
            pos_str += str(pos) + ' '
        result_str += pos_str + '\n'
    return result_str


with open(cluster_match_full_file_name, 'w') as cluster_match_file:
    for i, fasta_seq in enumerate(cluster_seqs):
        sequence = str(fasta_seq.seq)
        cluster_interval = re.split("=| ", fasta_seq.description)[2]
        search_result_list = MOODS.search(sequence, matrix_list, threshold_list, convert_log_odds=False,
                                          both_strands=True, pseudocount=0, threshold_from_p=False)

        sequence_len = len(sequence)
        cluster_match_file.write(result_to_str(cluster_interval, search_result_list, sequence_len))