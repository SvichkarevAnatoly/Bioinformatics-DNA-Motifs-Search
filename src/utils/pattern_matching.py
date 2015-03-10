import argparse
import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO
import Bio.motifs as motifs

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
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

    def fill_backward_matching(self, matching):
        for i, tf in enumerate(self.tfs):
            self.tf_dict[tf].backward = matching[i]


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
                             "Default stdout.")

    parser.add_argument("-tf", "--factor", nargs='+', dest="tf", type=str,
                        help="transcription factor name in pwm file. "
                             "Default matching with all tf in pwm file.")
    parser.add_argument("-th", "--threshold", dest="threshold", type=float, default=0.7,
                        help="The parameter threshold split for better control on what parts of the scoring are used. "
                             "Default 0.7.")
    parser.add_argument("-rc", "--reverse-complement", dest="reverse_complement", action="store_true", default=False,
                        help="For searching in both strands. "
                             "Default False.")
    parser.add_argument("-b", "--backward", dest="backward", action="store_true", default=False,
                        help="For searching in both direction. "
                             "Default only in direct.")
    parser.add_argument("-e", "--excel", dest="excel", action="store_true", default=False,
                        help="For saving results in easy paste to excel format. "
                             "Default human readable format.")
    return parser


def process(args):
    # TODO: write triggers in parser for writing input data
    seqs = list(SeqIO.parse(args.fasta, "fasta"))
    args.fasta.close()

    pwm_records = motifs.parse(args.pwm, "TRANSFAC")
    args.pwm.close()

    if args.tf is None:
        args.tf = lib.get_pwm_ids(pwm_records)
    pwms = lib.filter_pwms_in_tfs(pwm_records, args.tf)
    matrices = lib.create_matrices_from_pwms(pwms, args.tf)

    results = []
    for seq in seqs:
        sequence = str(seq.seq)
        matching = lib.search_motif(sequence, matrices, args.threshold, args.reverse_complement)

        seq_result = SeqSearchResults(seq.description)
        seq_result.create_tf_dict(args.tf)
        seq_result.fill_directed_matching(matching)

        if args.backward:
            backward_sequence = seq.seq[::-1]
            backward_matching = lib.search_motif(backward_sequence, matrices, args.threshold, args.reverse_complement)
            seq_result.fill_backward_matching(backward_matching)

        results.append(seq_result)
    return results


def save_excel(result, args):
    for seq_result in result:
        args.output.write('[' + seq_result.seq_name + ']')
        for tf in seq_result.tfs:
            matching = seq_result.tf_dict[tf]
            positions = [pos for pos, score in matching.directed]
            if positions:
                positions_str = ' ' + ';'.join(map(str, positions))
            else:
                positions_str = " #"
            args.output.write(positions_str)
            if hasattr(matching, 'backward'):
                args.output.write('|')
                backward_positions = [pos_tuple[0] for pos_tuple in matching.backward]
                if backward_positions:
                    backward_positions_str = ';'.join(map(str, backward_positions))
                else:
                    backward_positions_str = "#"
                args.output.write(backward_positions_str)
        args.output.write('\n')


def save_human_readable(result, args):
    for seq_result in result:
        args.output.write('>' + seq_result.seq_name + '\n')
        for tf in seq_result.tfs:
            args.output.write(tf + ' ')
            matching = seq_result.tf_dict[tf]
            positions = [pos_tuple[0] for pos_tuple in matching.directed]
            args.output.write(';'.join(map(str, positions)) + '\n')

            if hasattr(matching, 'backward'):
                args.output.write('backward ')
                backward_positions = [pos_tuple[0] for pos_tuple in matching.backward]
                args.output.write(';'.join(map(str, backward_positions)) + '\n')


def save(result, args):
    if args.excel:
        save_excel(result, args)
    else:
        save_human_readable(result, args)


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)
    args.output.close()


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()