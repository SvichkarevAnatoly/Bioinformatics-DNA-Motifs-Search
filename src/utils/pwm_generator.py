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


def create_parser():
    parser = argparse.ArgumentParser(
        description="Create position weight matrices (PWM) from DNA sequences")
    parser.add_argument("seqs", type=argparse.FileType('r'),
                        action=ReadSeqAction,
                        help="file with DNA sequences same length."
                             " One sequence in on line")
    parser.add_argument("-o", "--output", nargs='?', dest="output",
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="output file with PWM")
    parser.add_argument("-m", "--motif", nargs='?', dest="motif",
                        type=str, default=None,
                        help="pwm motif name. Default no name.")

    return parser


def process(args):
    lib.check_seq_correct(args.seqs)
    motif = motifs.create(args.seqs)
    return motif


def save(result, args):
    if args.motif is not None:
        args.output.write("ID  " + args.motif + '\n')
    args.output.write(result.format("transfac"))


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)
    args.output.close()


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()
