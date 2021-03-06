import argparse
import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO
import Bio.motifs as motifs
import MOODS
from bbcflib.track import track

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import lib


class SeqSearchResults(object):
    def __init__(self, seq_name, sequence, peak, tf_names, tf_lengths):
        self.seq_name = seq_name
        self.sequence = sequence
        self.peak = peak
        self.tfs = tf_names
        self.tf_length_dict = {tf: length for tf, length in zip(tf_names, tf_lengths)}
        self.tf_dict = {tf: [] for tf in tf_names}
        self.tf_max_scores_dict = {tf: 0 for tf in tf_names}

    def fill_matches(self, matches, max_scores):
        for i, tf in enumerate(self.tfs):
            self.tf_dict[tf] = matches[i]
            self.tf_max_scores_dict[tf] = max_scores[i]

    def best_matches(self, tf_name):
        tf_len = self.tf_length_dict[tf_name]
        matches = self.tf_dict[tf_name]
        seq_len = len(self.sequence)
        best_matches = []
        best_match = (0, sys.float_info.min)
        for match in matches:
            if match[1] > best_match[1]:
                best_match = match
                best_matches = [best_match]
            elif match[1] == best_match[1]:
                match_pos = match[0] if match[0] >= 0 else match[0] + seq_len
                best_match_pos = best_match[0] if best_match[0] >= 0 else best_match[0] + seq_len

                match_dist_to_border = min(match_pos, seq_len - (match_pos + tf_len))
                best_match_to_border = min(best_match_pos, seq_len - (best_match_pos + tf_len))

                if match_dist_to_border > best_match_to_border:
                    best_match = match
                    best_matches = [best_match]
                elif match_dist_to_border == best_match_to_border:
                    best_matches.append(match)
        return best_matches

    def match_subseq(self, pos, tf_name, delta):
        tf_len = self.tf_length_dict[tf_name]
        if pos >= 0:
            left = max(pos - delta, 0)
            right = min(pos + tf_len + delta, len(self.sequence))
            subseq = self.sequence[left: right]
        else:
            pos += len(self.sequence)
            rc_seq = lib.complement(self.sequence)
            left = max(pos - delta, 0)
            right = min(pos + tf_len + delta, len(rc_seq))
            subseq = rc_seq[left: right]
            subseq = subseq[::-1]
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


class ReadBedAction(argparse.Action):
    def __call__(self, parser, args, bed_handler, option_string=None):
        all_bed_fields = ['chrom', 'chromStart', 'chromEnd',
                          'peakName',
                          'c1']
        select_fields = ['chrom', 'chromStart', 'chromEnd', 'peakName']

        with track(bed_handler.name, format='txt',
                   separator='\t', fields=all_bed_fields) as t:
            bed_peaks = t.read(fields=select_fields)

        setattr(args, self.dest, bed_peaks)


class UpperCaseAction(argparse.Action):
    def __call__(self, parser, args, tf_names, option_string=None):
        tf_names = [val.upper() for val in tf_names]
        setattr(args, self.dest, tf_names)


def create_parser():
    parser = argparse.ArgumentParser(description="Matching position frequency matrices (PFM) against DNA sequences")
    parser.add_argument("fasta", type=argparse.FileType('r'),
                        action=ReadFastaAction,
                        help="fasta file with DNA sequences")
    parser.add_argument("pwm", type=argparse.FileType('r'),
                        action=ReadPWMAction,
                        help="file with position weight matrices (PWM)")
    parser.add_argument("bed", type=argparse.FileType('r'),
                        action=ReadBedAction,
                        help="bed file with peaks")

    parser.add_argument("-o", "--output", nargs='?', dest="output",
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="output bed file with extended matching results. "
                             "Default stdout.")
    parser.add_argument("-oc", "--outcompact", nargs='?', dest="output_compact",
                        type=argparse.FileType('w'), default=None,
                        help="output bed file with compact matching results. "
                             "Default not output.")

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
                             "at position [position - sequence_length] "
                             "in complement of input sequence, which is always "
                             "negative. The actual hit site for any hit is always "
                             "seq[pos, pos + matrix_length]. "
                             "Default False.")
    parser.add_argument("-c", "--constriction", dest="constriction", type=int, default=100,
                        help="Distance at which peak-caller MACS (and MACS2) "
                             "expands the boundaries of localization. "
                             "Default 100.")
    return parser


def process(args):
    if args.tf is None:
        args.tf = lib.get_pwm_id_names(args.pwm)

    pwms = lib.filter_pwms_in_tfs(args.pwm, args.tf)
    matrices = lib.create_matrices_from_pwms(pwms, args.tf)
    tf_lengths = [len(m[0]) for m in matrices]
    max_scores = [MOODS.max_score(m) for m in matrices]

    results = []
    for seq in args.fasta:
        peak = args.bed.next()
        if peak[3] != seq.description:
            raise Exception("seq header in fasta not equal to peak name in bed")

        sequence = str(seq.seq)
        matches = lib.search_motif(sequence, matrices, args.threshold, args.reverse_complement)

        seq_result = SeqSearchResults(seq.description, sequence, peak[0:3], args.tf, tf_lengths)
        seq_result.fill_matches(matches, max_scores)

        results.append(seq_result)
    return results


def check_filtered(match, best_matches, seq_length, tf_len, constriction):
    if match in best_matches:
        if seq_length <= 2*constriction:
            return 1
        else:
            pos = lib.local_pos(match[0], seq_length)
            if pos >= constriction and (pos + tf_len < seq_length - constriction):
                return 1
            else:
                return 0
    else:
        return 0


def write_compact(args, match_info, tf_length):
    compact_match_info = list(match_info[:6])
    start_peak = int(compact_match_info[1])
    start_bs = int(match_info[6])
    compact_match_info[1] = start_peak + start_bs
    compact_match_info[2] = compact_match_info[1] + tf_length
    compact_match_line = '\t'.join(map(str, compact_match_info))
    args.output_compact.write(compact_match_line + '\n')


def save(result, args):
    for seq_result in result:
        seq_length = len(seq_result.sequence)
        general_info = list(seq_result.peak)
        general_info.append(seq_result.seq_name + "_BS_")
        bs_index = 1
        for tf in seq_result.tfs:
            matches_tf = seq_result.tf_dict[tf]
            if not matches_tf:
                continue
            tf_length = seq_result.tf_length_dict[tf]
            max_score_tf = seq_result.tf_max_scores_dict[tf]
            best_matches = seq_result.best_matches(tf)
            for match in matches_tf:
                match_info = list(general_info)
                match_info[3] += str(bs_index)
                bs_index += 1

                score = match[1] / max_score_tf
                match_info.append("%.5f" % score)
                pos = match[0]
                if pos >= 0:
                    match_info.append('+')
                else:
                    match_info.append('-')

                local_pos = lib.local_pos(pos, seq_length)
                match_info.append(local_pos)
                match_info.append(local_pos + tf_length)

                filtered = check_filtered(match, best_matches, seq_length, tf_length, args.constriction)
                if args.output_compact is not None:
                    if filtered == 1:
                        write_compact(args, match_info, tf_length)

                match_info.append(filtered)

                predicted_site_seq = seq_result.match_subseq(match[0], tf, 0)
                match_info.append(predicted_site_seq)

                match_line = '\t'.join(map(str, match_info))
                args.output.write(match_line + '\n')


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)
    args.output.close()


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()