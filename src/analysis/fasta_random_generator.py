import argparse
import random
from signal import signal, SIGPIPE, SIG_DFL
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


def to_nucleotide(number):
    return {0: 'A',
            1: 'C',
            2: 'G',
            3: 'T'}[number]


def save(args):
    for i in range(args.number):
        seq = ""
        for j in range(args.length):
            seq += to_nucleotide(random.randrange(4))

        id = "seq" + str(i)
        seq = Seq(seq)
        seq_record = SeqRecord(seq, id, description="")

        SeqIO.write(seq_record, args.output, "fasta")


def main():
    parser = create_parser()
    args = parser.parse_args()
    save(args)
    args.output.close()


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()
