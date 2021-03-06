import os
import re

from Bio.Alphabet import IUPAC
from Bio import Alphabet
from Bio.Seq import Seq
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


def get_pwm_id_names(pwm_records):
    return [pwm['ID'].upper() for pwm in pwm_records]


def filter_pwms_in_tfs(pwms, tf_names):
    pwm_ids = [pwm_name.upper() for pwm_name in get_pwm_id_names(pwms)]
    tfs_set = set(tf_names)
    intersection = tfs_set & set(pwm_ids)
    if tfs_set == intersection:
        return [pwm for pwm in pwms if pwm['ID'].upper() in tf_names]
    else:
        # TODO: find exception class
        raise Exception("I know python!")


def create_matrices_from_pwms(pwms, tf_names):
    if tf_names is not None:
        pwms = [pwm for pwm in pwms if pwm['ID'].upper() in tf_names]
    matrices = [record.counts for record in pwms]
    alp = dna_alphabet()
    matrices = [[matrix[n] for n in alp] for matrix in matrices]
    return matrices


def search_motif(sequence, matrices, threshold_factor=0.7, both_strands=False):
    threshold_func = lambda m: threshold_factor * MOODS.max_score(m)
    threshold = map(threshold_func, matrices)
    return MOODS.search(sequence, matrices, threshold,
                        convert_log_odds=False, both_strands=both_strands,
                        pseudocount=0, threshold_from_p=False,
                        algorithm="naive")


def get_join_position_str(positions, seq_length):
    positions_str = []
    for pos in positions:
        if pos >= 0:
            positions_str.append(str(pos))
        else:
            pos += seq_length
            positions_str.append(str(pos) + "(-)")
    return ';'.join(positions_str)


def local_pos(pos, seq_length):
    if pos < 0:
        pos += seq_length
    return pos


def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


def complement(sequence):
    return str(Seq(sequence).complement())


def biopython_seqs(sequences):
    return [Seq(seq, IUPAC.unambiguous_dna) for seq in sequences]


def check_seq_correct(seqs):
    pwm_length = len(seqs[0])
    for seq in seqs:
        if not Alphabet._verify_alphabet(seq):
            raise ValueError("Seq contains letter not from Alphabet")
        if len(seq) != pwm_length:
            raise ValueError("Seqs length differ")