import argparse
from signal import signal, SIGPIPE, SIG_DFL

from Bio.Alphabet import IUPAC
from Bio import Alphabet
from Bio.Seq import Seq
from Bio import motifs


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
        description="Create sequence logo from sequences")
    parser.add_argument("seqs", type=argparse.FileType('r'),
                        action=ReadSeqAction,
                        help="file with DNA sequences same length."
                             " One sequence in on line")
    parser.add_argument("-o", "--output", required=True,
                        dest="output", type=argparse.FileType('w'),
                        help="output file with logo")

    return parser


def check(seqs):
    pwm_length = len(seqs[0])
    for seq in seqs:
        if not Alphabet._verify_alphabet(seq):
            raise ValueError("Seq contains letter not from Alphabet")
        if len(seq) != pwm_length:
            raise ValueError("Seqs length differ")


def process(args):
    check(args.seqs)
    motif = motifs.create(args.seqs)
    return motif


def save(result, args):
    result.weblogo(args.output.name, format='SVG')


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)
    args.output.close()


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()

