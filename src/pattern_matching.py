import argparse
from Bio import SeqIO
import Bio.motifs as motifs
import lib


def create_parser():
    parser = argparse.ArgumentParser(description="Matching position weight matrices (PWM) against DNA sequences")
    parser.add_argument("fasta", type=argparse.FileType('r'), help="fasta file with DNA sequences")
    parser.add_argument("pwm", type=argparse.FileType('r'), help="file with position weight matrices (PWM)")
    parser.add_argument("-o", "--output", dest="matching", type=argparse.FileType('w'), metavar='matching',
                        help="output file with matching results")

    parser.add_argument("-tf", "--factor", dest="tf", type=int, metavar='tf',
                        help="transcription factor name in pwm file. "
                        "If not specified, matching with all tf in pwm file.")
    return parser


def process(args):
    # TODO: write triggers in parser for writing input data
    seq_list = list(SeqIO.parse(args.fasta, "fasta"))
    pwm_record_list = motifs.parse(args.pwm, "TRANSFAC")
    if args.tf is not None:
        # TODO: if tf not exist in pwm_record_list
        matrix = lib.create_matrix_from_pwms(args.tf, pwm_record_list)

    return None


def save(result, args):
    # TODO:
    pass


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)


if __name__ == "__main__":
    main()