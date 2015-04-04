import argparse
import os
from signal import signal, SIGPIPE, SIG_DFL
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import lib


def create_parser():
    parser = argparse.ArgumentParser(
        description="Central extension each interval of the specified file to the same length")
    parser.add_argument("bedfile", type=argparse.FileType('r'), help="file with bed format intervals")
    parser.add_argument("-l", "--length", dest="length", type=int,
                        help="common extended length. "
                             "If not specified, is extended to the maximum length of the interval in the input file.")
    parser.add_argument("-o", "--output", nargs='?', dest="output",
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="output file with extended bed format intervals. "
                             "If not specified, write output to stdout.")
    return parser


def interval_center_extender(intervals, length):
    intervals = map(lib.parse_interval_line, intervals)

    if length is not None:
        interval_length = length
    else:
        interval_length = max(map(lib.interval_length, intervals))

    return map(lambda interval: lib.interval_extend(interval, interval_length), intervals)


def process(args):
    intervals = args.bedfile.readlines()
    args.bedfile.close()
    return interval_center_extender(intervals, args.length)


def save(result, args):
    for interval_param in result:
        args.output.write("%s\n" % lib.interval_param_list_to_str(interval_param))


def main():
    parser = create_parser()
    args = parser.parse_args()
    result = process(args)
    save(result, args)
    args.output.close()


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    main()