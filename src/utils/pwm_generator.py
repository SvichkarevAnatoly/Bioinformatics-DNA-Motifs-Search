import argparse
from signal import signal, SIGPIPE, SIG_DFL
import sys


def create_parser():
    parser = argparse.ArgumentParser(description="TODO")
    parser.add_argument("seqs", type=argparse.FileType('r'),
                        help="TODO")
    parser.add_argument("-o", "--output", nargs='?', dest="output",
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="TODO")
    return parser


def process(args):
    pass


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
