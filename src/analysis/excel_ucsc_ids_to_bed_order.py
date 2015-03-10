import argparse
import os
import re
import sys
from signal import signal, SIGPIPE, SIG_DFL

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))


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

    lines = args.excel.readlines()
    result_lines = []
    for line in lines:
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
