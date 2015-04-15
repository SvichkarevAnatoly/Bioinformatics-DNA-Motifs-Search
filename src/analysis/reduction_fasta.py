import argparse
from signal import signal, SIGPIPE, SIG_DFL
import sys

from Bio import SeqIO
from bbcflib.track import track


class ReadFastaAction(argparse.Action):
    def __call__(self, parser, args, fasta_handler, option_string=None):
        seqs = list(SeqIO.parse(fasta_handler, "fasta"))
        fasta_handler.close()
        setattr(args, self.dest, seqs)


class ReadBedAction(argparse.Action):
    def __call__(self, parser, args, bed_handler, option_string=None):
        all_bed_fields = ['chrom', 'chromStart', 'chromEnd',
                          'c1', 'c2', 'c3', 'c4', 'c5', 'c6',
                          'peakName']
        select_fields = ['chrom', 'chromStart', 'chromEnd', 'peakName']

        with track(bed_handler.name, format='txt',
                   separator='\t', fields=all_bed_fields) as t:
            bed_peaks = t.read(fields=select_fields)

        setattr(args, self.dest, bed_peaks)


def create_parser():
    parser = argparse.ArgumentParser(
        description="Order fasta according to bed file, change header to peak name")
    parser.add_argument("fasta", type=argparse.FileType('r'),
                        action=ReadFastaAction,
                        help="fasta file with DNA sequences")
    parser.add_argument("bed", type=argparse.FileType('r'),
                        action=ReadBedAction,
                        help="bed file with peaks")
    parser.add_argument("-o", "--output", nargs='?', dest="output",
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="output file with fasta")
    return parser


def save(args):
    for peak in args.bed:
        start = str(int(peak[1]) + 1)
        range_interval = peak[0] + ':' + start + '-' + peak[2]
        fasta_record = None
        for fasta in args.fasta:
            if range_interval in fasta.description:
                fasta_record = fasta
                break

        if fasta_record is None:
            raise Exception("fasta not found")

        fasta_record.id = peak[3]
        fasta_record.description = ""
        SeqIO.write(fasta_record, args.output, "fasta")


def main():
    parser = create_parser()
    args = parser.parse_args()
    save(args)
    args.output.close()


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()
