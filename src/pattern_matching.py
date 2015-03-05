import argparse

from Bio import SeqIO
import Bio.motifs as motifs
import sys

import lib


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

        result_cortege = (seq.description, zip(tf_names, matching))
        if args.reversed:
            reversed_sequence = seq.seq[::-1]
            reversed_matching = lib.search_motif(reversed_sequence, matrices, args.threshold, args.reversed)
            # TODO: test
            result_cortege = (result_cortege, reversed_matching)

        results.append(result_cortege)
    return results


def save(result, args):
    for seq_result in result:
        seq_name = seq_result[0]
        args.output.write('>' + seq_name + '\n')
        for tf_result in seq_result[1]:
            tf_name = tf_result[0]
            args.output.write(tf_name + ' ')
            positions = [pos_tuple[0] for pos_tuple in tf_result[1]]
            args.output.write(';'.join(map(str, positions)) + '\n')
    args.output.close()


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)


if __name__ == "__main__":
    main()