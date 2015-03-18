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


def generate_pwm_str(motif_name, sequence):
    pwm_str = "ID " + motif_name + "\n" \
              "P0  A C G T\n"
    for i, nucleotide in enumerate(sequence):
        counters = [0] * 4
        counters[to_ind(nucleotide)] = 9
        pwm_str += str(i+1) + "   " + " ".join(map(str, counters)) + '\n'
    pwm_str += "//\n"
    return pwm_str


def create_args(sequence, pwm_str=None):
    args = Namespace()
    args.pwm = create_pwm(pwm_str)
    args.fasta = create_fasta(sequence)
    args.output = cStringIO.StringIO()
    args.tf = None
    args.reverse_complement = False
    args.excel = False
    args.threshold = 0.7
    return args


def read_output_file(output):
    output.seek(0)
    return output.read()


def silent_remove(file_name):
    try:
        os.remove(file_name)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured

