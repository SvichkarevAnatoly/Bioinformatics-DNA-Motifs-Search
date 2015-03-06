import argparse
import sys
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO
import Bio.motifs as motifs

import lib


class SeqSearchResults(object):
    def __init__(self, seq_name):
        self.seq_name = seq_name
        self.tf_dict = None
        self.tfs = None

    def create_tf_dict(self, tf_names):
        self.tfs = tf_names
        self.tf_dict = {tf: DirectionMatchingTF(tf) for tf in tf_names}

    def fill_directed_matching(self, matching):
        for i, tf in enumerate(self.tfs):
            self.tf_dict[tf].directed = matching[i]

    def fill_reversed_matching(self, matching):
        for i, tf in enumerate(self.tfs):
            self.tf_dict[tf].reversed = matching[i]


class DirectionMatchingTF(object):
    def __init__(self, tf_name):
        self.tf = tf_name


def create_parser():
    parser = argparse.ArgumentParser(description="Matching position weight matrices (PWM) against DNA sequences")
    parser.add_argument("fasta", type=argparse.FileType('r'), help="fasta file with DNA sequences")
    parser.add_argument("pwm", type=argparse.FileType('r'), help="file with position weight matrices (PWM)")
    parser.add_argument("-o", "--output", nargs='?', dest="output",
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="output file with matching results. "
                             "If not specified, write output to stdout.")

    # TODO: make for list tf_names
    parser.add_argument("-tf", "--factor", dest="tf", type=str,
                        help="transcription factor name in pwm file. "
                             "If not specified, matching with all tf in pwm file.")
    parser.add_argument("-th", "--threshold", dest="threshold", type=float, default=0.7,
                        help="The parameter threshold split for better control on what parts of the scoring are used. "
                             "If not specified, threshold=0.7.")
    parser.add_argument("-r", "--reversed", action="store_true", default=False,
                        help="For searching in both direction. "
                             "If not specified, search only in direct.")
    parser.add_argument("-e", "--excel", action="store_true", default=False,
                        help="For saving results in easy paste to excel format. "
                             "If not specified, saving results in compact format.")
    return parser


def process(args):
    # TODO: write triggers in parser for writing input data
    seqs = list(SeqIO.parse(args.fasta, "fasta"))
    args.fasta.close()

    pwm_record_list = motifs.parse(args.pwm, "TRANSFAC")
    args.pwm.close()

    matrices, tf_names = lib.create_matrices_from_pwms(pwm_record_list, args.tf)

    results = []
    for seq in seqs:
        sequence = str(seq.seq)
        matching = lib.search_motif(sequence, matrices, args.threshold, args.reversed)

        seq_result = SeqSearchResults(seq.description)
        seq_result.create_tf_dict(tf_names)
        seq_result.fill_directed_matching(matching)

        if args.reversed:
            reversed_sequence = seq.seq[::-1]
            reversed_matching = lib.search_motif(reversed_sequence, matrices, args.threshold, args.reversed)
            seq_result.fill_reversed_matching(reversed_matching)

        results.append(seq_result)
    return results


def save(result, args):
    for seq_result in result:
        args.output.write('>' + seq_result.seq_name + '\n')
        for tf in seq_result.tfs:
            args.output.write(tf + ' ')
            directed_matching = seq_result.tf_dict[tf].directed
            positions = [pos_tuple[0] for pos_tuple in directed_matching]
            args.output.write(';'.join(map(str, positions)) + '\n')
    args.output.close()


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()