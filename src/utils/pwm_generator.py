import argparse
from signal import signal, SIGPIPE, SIG_DFL
import sys

from Bio.Alphabet import IUPAC
from Bio import Alphabet
from Bio.Seq import Seq


class ReadSeqAction(argparse.Action):
    def __call__(self, parser, args, seqs_handler, option_string=None):
        seqs = []
        for seq_line in seqs_handler:
            seq_string = seq_line.strip()
            seq = Seq(seq_string, IUPAC.unambiguous_dna)
            Alphabet._verify_alphabet(seq)
            # TODO: check same length
            seqs.append(seq)
        seqs_handler.close()
        setattr(args, self.dest, seqs)


def create_parser():
    parser = argparse.ArgumentParser(description="TODO")
    parser.add_argument("seqs", type=argparse.FileType('r'),
                        action=ReadSeqAction,
                        help="TODO")
    parser.add_argument("-o", "--output", nargs='?', dest="output",
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="TODO")
    return parser


def check(seqs):
    pwm_length = seqs[0]
    for seq in seqs:
        Alphabet._verify_alphabet(seq)
        if len(seq) != pwm_length:
            raise ValueError("seqs length differ")


def process(args):
    check(args.seqs)


def save(result, args):
    pass


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)
    args.output.close()


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()
