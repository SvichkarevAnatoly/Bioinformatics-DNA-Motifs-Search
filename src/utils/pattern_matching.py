import argparse
import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO
import Bio.motifs as motifs

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import lib


class SeqSearchResults(object):
    def __init__(self, seq_name, sequence, tf_names, tf_lengths):
        self.seq_name = seq_name
        self.sequence = sequence
        self.tfs = tf_names
        self.tf_length_dict = {tf: length for tf, length in zip(tf_names, tf_lengths)}
        self.tf_dict = {tf: [] for tf in tf_names}

    def fill_matches(self, matches):
        for i, tf in enumerate(self.tfs):
            self.tf_dict[tf] = matches[i]

    def best_match(self, tf_name):
        tf_len = self.tf_length_dict[tf_name]
        matches = self.tf_dict[tf_name]
        half_seq_len = len(self.sequence) / 2
        best_match = (0, sys.float_info.min)
        for match in matches:
            if match[1] > best_match[1]:
                best_match = match
            elif match[1] == best_match[1]:
                match_pos = match[0] if match[0] >= 0 else match[0] + tf_len
                best_match_pos = best_match[0] if best_match[0] >= 0 else best_match[0] + tf_len
                match_dist_to_center = abs(half_seq_len - (match_pos + tf_len / 2))
                best_match_dist_to_center = abs(half_seq_len - (best_match_pos + tf_len / 2))
                if match_dist_to_center < best_match_dist_to_center:
                    best_match = match
        return best_match

    def match_subseq(self, pos, tf_name, delta):
        tf_len = self.tf_length_dict[tf_name]
        if pos >= 0:
            left = max(pos - delta, 0)
            right = min(pos + tf_len + delta, len(self.sequence))
            subseq = self.sequence[left: right]
        else:
            pos += len(self.sequence)
            rc_seq = lib.reverse_complement(self.sequence)
            left = max(pos - delta, 0)
            right = min(pos + tf_len + delta, len(self.sequence))
            subseq = rc_seq[left: right]
        return subseq


class ReadFastaAction(argparse.Action):
    def __call__(self, parser, args, fasta_handler, option_string=None):
        seqs = list(SeqIO.parse(fasta_handler, "fasta"))
        fasta_handler.close()
        setattr(args, self.dest, seqs)


class ReadPWMAction(argparse.Action):
    def __call__(self, parser, args, pwm_handler, option_string=None):
        pwm_records = motifs.parse(pwm_handler, "TRANSFAC")
        pwm_handler.close()
        setattr(args, self.dest, pwm_records)


class UpperCaseAction(argparse.Action):
    def __call__(self, parser, args, tf_names, option_string=None):
        tf_names = [val.upper() for val in tf_names]
        setattr(args, self.dest, tf_names)


def create_parser():
    parser = argparse.ArgumentParser(description="Matching position weight matrices (PWM) against DNA sequences")
    parser.add_argument("fasta", type=argparse.FileType('r'),
                        action=ReadFastaAction,
                        help="fasta file with DNA sequences")
    parser.add_argument("pwm", type=argparse.FileType('r'),
                        action=ReadPWMAction,
                        help="file with position weight matrices (PWM)")
    parser.add_argument("-o", "--output", nargs='?', dest="output",
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="output file with matching results. "
                             "Default stdout.")
    parser.add_argument("-tf", "--factor", nargs='+', dest="tf", type=str,
                        action=UpperCaseAction,
                        help="transcription factor name in pwm file. "
                             "Default matching with all tf in pwm file.")
    parser.add_argument("-th", "--threshold", dest="threshold", type=float, default=0.7,
                        help="The parameter threshold split for better control on what parts of the scoring are used. "
                             "Default 0.7.")
    parser.add_argument("-rc", "--reverse-complement", dest="reverse_complement", action="store_true", default=False,
                        help="Scans against reverse complement sequence in addition to "
                             "the input sequence. Hits on reverse complement are reported "
                             "at position [position - sequence_length], which is always "
                             "negative. The actual hit site for any hit is always "
                             "seq[pos, pos + matrix_length]. "
                             "Default False.")
    parser.add_argument("-e", "--excel", dest="excel", action="store_true", default=False,
                        help="For saving results in easy paste to excel format. "
                             "Default human readable format.")
    return parser


def process(args):
    if args.tf is None:
        args.tf = lib.get_pwm_ids(args.pwm)

    pwms = lib.filter_pwms_in_tfs(args.pwm, args.tf)
    matrices = lib.create_matrices_from_pwms(pwms, args.tf)
    tf_lengths = [len(m[0]) for m in matrices]

    results = []
    for seq in args.fasta:
        sequence = str(seq.seq)
        matches = lib.search_motif(sequence, matrices, args.threshold, args.reverse_complement)

        seq_result = SeqSearchResults(seq.description, sequence, args.tf, tf_lengths)
        seq_result.fill_matches(matches)

        results.append(seq_result)
    return results


def save_excel(result, args):
    for seq_result in result:
        seq_length = len(seq_result.sequence)
        args.output.write('[' + seq_result.seq_name + ']')
        for tf in seq_result.tfs:
            matches_tf = seq_result.tf_dict[tf]
            positions = [pos for pos, score in matches_tf]
            positions_str = lib.get_join_position_str(positions, seq_length)

            best_match = seq_result.best_match(tf)
            best_subseq = seq_result.match_subseq(best_match[0], tf, 5)

            args.output.write(' ' + positions_str + ' ' + best_subseq)
        args.output.write('\n')


def save_human_readable(result, args):
    for seq_result in result:
        seq_length = len(seq_result.sequence)
        args.output.write('>' + seq_result.seq_name + '\n')
        for tf in seq_result.tfs:
            args.output.write(tf + ' ')
            matches_tf = seq_result.tf_dict[tf]
            positions = [pos_tuple[0] for pos_tuple in matches_tf]
            positions_str = lib.get_join_position_str(positions, seq_length)
            args.output.write(positions_str + '\n')


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