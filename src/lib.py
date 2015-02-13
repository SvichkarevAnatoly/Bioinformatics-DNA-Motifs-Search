import re
from Bio.Alphabet import IUPAC


def parse_interval_line(interval_line):
    chr_name, start, end = re.split(":|-", interval_line.strip())
    start = int(start)
    end = int(end)
    return [chr_name, start, end]


def interval_param_to_str(interval_param):
    return interval_param[0] + ':' + str(interval_param[1]) + '-' + str(interval_param[2])


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