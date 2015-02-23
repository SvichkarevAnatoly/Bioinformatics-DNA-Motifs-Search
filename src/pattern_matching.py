import argparse
from Bio import SeqIO
import Bio.motifs as motifs


def create_parser():
    parser = argparse.ArgumentParser(description="Matching position weight matrices (PWM) against DNA sequences")
    parser.add_argument("fasta", type=argparse.FileType('r'), help="fasta file with DNA sequences")
    parser.add_argument("pwm", type=argparse.FileType('r'), help="file with position weight matrices (PWM)")
    parser.add_argument("-o", "--output", dest="matching", type=argparse.FileType('w'), metavar='matching',
                        help="output file with matching results")
    return parser


def process(args):
    # TODO: write triggers in parser for writing input data
    seq_list = list(SeqIO.parse(args.fasta, "fasta"))
    pwm_records = motifs.parse(args.pwm, "TRANSFAC")

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