import argparse
import os
import re
import sys
from signal import signal, SIGPIPE, SIG_DFL

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import lib


def create_parser():
    parser = argparse.ArgumentParser(
        description="Converting matching results in plain text excel format in bed file orders")
    parser.add_argument("excel", type=argparse.FileType('r'),
                        help="text file with excel matching")
    parser.add_argument("-b", "--bedfile", nargs='?',
                        type=argparse.FileType('r'),
                        help="bed file with intervals. "
                             "Order identifiers to order in bed file. "
                             "Default not order.")
    parser.add_argument("-o", "--output", nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="output file with formatted matching results. "
                             "Default stdout.")
    return parser


def process(args):
    interval_pattern = re.compile(r"chr\w+:\d+-\d+")
    matching_pattern = re.compile(r"\]")

    excel_lines = args.excel.readlines()
    args.excel.close()

    result_lines = []
    if args.bed is not None:
        formatted_lines_dict = {}
        for line in excel_lines:
            interval_str = interval_pattern.search(line).group()
            interval = lib.parse_interval_line(interval_str)
            matching = matching_pattern.split(line)[1]

            interval[1] -= 1
            bed_interval = lib.interval_param_list_to_str(interval)
            formatted_lines_dict[bed_interval] = matching

        bed_lines = args.bed.read().splitlines()
        args.bed.close()

        for bed_line in bed_lines:
            result_lines.append(bed_line + formatted_lines_dict[bed_line])
    else:
        for line in excel_lines:
            interval = interval_pattern.search(line).group()
            matching = matching_pattern.split(line)[1]
            result_lines.append(interval + matching)

    return result_lines


def save(result, args):
    args.output.writelines(result)


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)
    args.output.close()


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()
