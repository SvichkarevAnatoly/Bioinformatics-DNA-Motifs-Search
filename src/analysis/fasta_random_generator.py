import argparse
import os
from signal import signal, SIGPIPE, SIG_DFL
import sys

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import motifs

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import lib


def create_parser():
    parser = argparse.ArgumentParser(
        description="Create random DNA sequences in fasta format")
    parser.add_argument("-o", "--output", nargs='?', dest="output",
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="output file with fasta")
    parser.add_argument("-l", "--length", dest="length", type=int,
                        default=10, help="length every fasta seq")
    parser.add_argument("-n", "--number", dest="number", type=int,
                        default=1, help="number fasta seqs")
    return parser


def process(args):
    return None


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
