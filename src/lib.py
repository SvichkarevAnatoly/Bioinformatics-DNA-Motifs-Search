import os
import re
from Bio.Alphabet import IUPAC
import MOODS


def parse_interval_line(interval_line):
    chr_name, start, end = re.split(":|-", interval_line.strip())
    start = int(start)
    end = int(end)
    return [chr_name, start, end]


def interval_param_list_to_str(interval_param):
    return interval_param[0] + ':' + str(interval_param[1]) + '-' + str(interval_param[2])


def interval_param_to_str(chr_name, start, end):
    return chr_name + ':' + str(start) + '-' + str(end)


def interval_length(interval_param):
    return interval_param[2] - interval_param[1]


def interval_extend(interval_param, new_length):
    start = interval_param[1]
    end = interval_param[2]
    center = (end + start) / 2
    new_start = center - new_length / 2
    new_end = center + new_length / 2
    return [interval_param[0], new_start, new_end]


# TODO: how to test input file
def read_excel_motif_matrix_list_from_file(input_file):
    dna_alp = dna_alphabet()
    id_list = []
    matrix_list = []
    matrix = {n: [] for n in dna_alp}
    for line in input_file:
        param_list = line.strip().split()
        if len(param_list) == 1:
            id_list.append(param_list[0][2:])
            matrix_list.append(matrix)
            matrix = {n: [] for n in dna_alp}
        else:
            if param_list[1].isdigit():
                for i, nucleotide in enumerate(dna_alp):
                    matrix[nucleotide].append(int(param_list[i + 1]))
    matrix_list.pop(0)
    matrix_list.append(matrix)
    return [id_list, matrix_list, matrix]


def dna_alphabet():
    return sorted(list(IUPAC.unambiguous_dna.letters))


def create_output_file_name(input_file_name):
    path, basename = os.path.split(input_file_name)
    basename_part_list = os.path.splitext(basename)
    return os.path.join(path, basename_part_list[0] + '_out' + basename_part_list[1])


# TODO: refactored
# TODO: make for list tf_names
def create_matrices_from_pwms(pwm_records, tf_names=None):
    if tf_names is not None:
        pwm_records = [pwm for pwm in pwm_records if pwm['ID'] == tf_names]
        tf_names = [tf_names]
    else:
        tf_names = [pwm['ID'] for pwm in pwm_records]
    matrices = [record.counts for record in pwm_records]
    alp = dna_alphabet()
    matrices = [[matrix[n] for n in alp] for matrix in matrices]
    return matrices, tf_names


def search_motif(sequence, matrices, threshold_factor=0.7, both_strands=False):
    threshold_func = lambda m: threshold_factor * MOODS.max_score(m)
    threshold = map(threshold_func, matrices)
    return MOODS.search(sequence, matrices, threshold,
                        convert_log_odds=False, both_strands=both_strands,
                        pseudocount=0, threshold_from_p=False)


def searching_result_to_str(interval, results, rev_results, sequence_length):
    result_str = ""
    if len(results[0]) == 0 and len(rev_results[0]) == 0:
        return result_str
    else:
        result_str += ">" + interval + '\n'
        for result in results:
            pos_str = ""
            for (pos, score) in result:
                if pos < 0:
                    pos_str += '-'
                    pos += sequence_length
                pos_str += str(pos) + ' '
            result_str += pos_str + '\n'
        for result in rev_results:
            pos_str = ""
            for (pos, score) in result:
                if pos < 0:
                    pos_str += '-'
                    pos += sequence_length
                pos_str += str(pos) + ' '
            result_str += pos_str + '\n'
    return result_str
