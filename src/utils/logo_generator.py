import argparse
import os
from signal import signal, SIGPIPE, SIG_DFL
import sys

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import motifs


sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import lib


class ReadSeqAction(argparse.Action):
    def __call__(self, parser, args, seqs_handler, option_string=None):
        seqs = []
        for seq_line in seqs_handler:
            seq_string = seq_line.strip()
            seq = Seq(seq_string, IUPAC.unambiguous_dna)
            seqs.append(seq)
        seqs_handler.close()
        setattr(args, self.dest, seqs)


class ReadPwmAction(argparse.Action):
    def __call__(self, parser, args, pwm_handler, option_string=None):
        pwm = motifs.parse(pwm_handler, "TRANSFAC")
        pwm_handler.close()
        setattr(args, self.dest, pwm)


def create_parser():
    parser = argparse.ArgumentParser(
        description="Create sequence logo from sequences")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--seqs", type=argparse.FileType('r'),
                       action=ReadSeqAction,
                       help="file with DNA sequences same length."
                            " One sequence in on line")
    group.add_argument("-p", "--pwm", type=argparse.FileType('r'),
                       action=ReadPwmAction,
                       help="file with PWM.")
    parser.add_argument("-o", "--output", required=True,
                        dest="output", type=argparse.FileType('w'),
                        help="output file with logo in vector SVG format")

    return parser


def process(args):
    if args.seqs is not None:
        lib.check_seq_correct(args.seqs)
        motif = motifs.create(args.seqs)
    else:
        motif = args.pwm[0]
    return motif


def save(result, args):
    result.weblogo(args.output.name,
                   format='SVG',
                   color_scheme='color_classic')


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)
    args.output.close()


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()

