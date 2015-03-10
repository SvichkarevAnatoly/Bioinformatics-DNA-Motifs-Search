import argparse
import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))


def create_parser():
    parser = argparse.ArgumentParser(
        description="Converting matching results in plain text excel format in bed file orders")
    parser.add_argument("excel_file", type=argparse.FileType('r'),
                        help="text file with excel matching")
    parser.add_argument("-o", "--output", nargs='?', dest="output",
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="output file with formatted matching results. "
                             "Default stdout.")
    return parser


def process(args):
    results = []
    return results


# TODO
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
