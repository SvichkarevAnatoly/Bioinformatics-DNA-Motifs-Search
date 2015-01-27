import os
import re

from Bio import SeqIO
import Bio.motifs as motifs
from Bio.Alphabet import IUPAC
import MOODS


dna_alf = sorted(list(IUPAC.unambiguous_dna.letters))

pwms_file_name = "../data/PWMs_TRANSFAC.txt"
cluster_fasta_file_name = "../data/clusters/mm10_clusters.fa"
cluster_bed_file_name = "../data/clusters/mm10_clusters_ucsc.bed"

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
cluster_excel_file_name = cluster_file_basename_list[0] + '_excel.txt'
cluster_excel_full_file_name = os.path.join(cluster_file_path, cluster_excel_file_name)


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


def result_to_excel_str(interval, result_list, sequence_length):
    result_str = interval
    for ind, result in enumerate(result_list):
        if not result:
            result_str += " #"
        else:
            result_str += ' '
            for i, (pos, score) in enumerate(result):
                if pos < 0:
                    result_str += '-'
                    pos += sequence_length
                result_str += str(pos)
                if len(result) > i + 1:
                    result_str += ';'
    result_str += '\n'
    return result_str


result_list_list = []
sequence_len_list = []
cluster_interval_str_list = []
for i, fasta_seq in enumerate(cluster_seqs):
    sequence = str(fasta_seq.seq)
    search_result_list = MOODS.search(sequence, matrix_list, threshold_list, convert_log_odds=False,
                                      both_strands=True, pseudocount=0, threshold_from_p=False)
    result_list_list.append(search_result_list)
    sequence_len_list.append(len(sequence))
    cluster_interval_str = re.split("=| ", fasta_seq.description)[2]
    cluster_interval_str_list.append(cluster_interval_str)

with open(cluster_match_full_file_name, 'w') as cluster_match_file:
    for i, result in enumerate(result_list_list):
        cluster_match_file.write(result_to_str(cluster_interval_str_list[i], result, sequence_len_list[i]))

cluster_result_dict = {}
for i, result in enumerate(result_list_list):
    cluster_result_dict[cluster_interval_str_list[i]] = result_to_excel_str(cluster_interval_str_list[i], result,
                                                                            sequence_len_list[i])

with open(cluster_bed_file_name, 'r') as cluster_bed_file:
    with open(cluster_excel_full_file_name, 'w') as cluster_excel_file:
        for interval_line in cluster_bed_file:
            (chr_str, start, end) = interval_line.strip().split()
            start_int = int(start)
            interval_str = chr_str + ':' + str(start_int + 1) + '-' + end
            cluster_excel_file.write(cluster_result_dict[interval_str])