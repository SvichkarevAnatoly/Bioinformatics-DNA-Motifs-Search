from argparse import Namespace
import cStringIO
import os
from Bio import SeqIO
from Bio import motifs
import errno


def create_fasta(sequence):
    fasta_str = '\n'.join([
        ">seq",
        sequence
    ]) + '\n'
    fasta_handler = cStringIO.StringIO()
    fasta_handler.write(fasta_str)
    fasta_handler.seek(0)
    fasta = list(SeqIO.parse(fasta_handler, "fasta"))
    fasta_handler.close()
    return fasta


def create_pwm(pwm_str):
    pwm_handler = cStringIO.StringIO()
    pwm_handler.write(pwm_str)
    pwm_handler.seek(0)
    pwm_records = motifs.parse(pwm_handler, "TRANSFAC")
    pwm_handler.close()
    return pwm_records


def to_ind(nucleotide):
    return {'A': 0,
            'C': 1,
            'G': 2,
            'T': 3}[nucleotide]


def generate_simple_pwm_str(motif_name, sequence):
    pwm_str = "ID  " + motif_name + "\n" \
              "P0  A C G T\n"
    for i, nucleotide in enumerate(sequence):
        counters = [0] * 4
        counters[to_ind(nucleotide)] = 9
        pwm_str += str(i+1) + "   " + " ".join(map(str, counters)) + '\n'
    pwm_str += "//\n"
    return pwm_str


def generate_pwm_str(motif_name, pwm_matrix):
    pwm_str = "ID  " + motif_name + "\n" \
              "P0  A C G T\n"
    for i, counters in enumerate(pwm_matrix):
        pwm_str += str(i+1) + "   " + " ".join(map(str, counters)) + '\n'
    pwm_str += "//\n"
    return pwm_str


def create_args(sequence, pwm_str, is_excel=True, reverse_complement=False, threshold=0.7, tf_names=None):
    args = Namespace()
    args.pwm = create_pwm(pwm_str)
    args.fasta = create_fasta(sequence)
    args.output = cStringIO.StringIO()
    args.tf = tf_names
    args.reverse_complement = reverse_complement
    args.excel = is_excel
    args.threshold = threshold
    return args


def read_output_file(output):
    output.seek(0)
    return output.read()


def get_score(sequence, matrix):
    score = 0
    for i, nucleotide in enumerate(sequence):
        nucleotide_ind = to_ind(nucleotide)
        score += matrix[nucleotide_ind][i]
    return score


def silent_remove(file_name):
    try:
        os.remove(file_name)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured

