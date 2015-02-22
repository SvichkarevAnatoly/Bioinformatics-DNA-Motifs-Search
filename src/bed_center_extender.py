import argparse
import lib


def parse_args():
    parser = argparse.ArgumentParser(description="Resulting intervals in the input file to the same length")
    parser.add_argument("bedfile", type=argparse.FileType('r'), help="file with bed format intervals")
    parser.add_argument("-l", dest="length", type=int, metavar='length', help="common extended length")
    # TODO: remove outfile if args error
    parser.add_argument("-o", dest="outfile", type=argparse.FileType('w'), metavar='outfile',
                        help="output file with extended bed format intervals")
    return parser.parse_args()


def interval_center_extender(bed_file, length):
    interval_list = bed_file.readlines()
    interval_list = map(lib.parse_interval_line, interval_list)

    if length is not None:
        interval_length = length
    else:
        interval_length = max(map(lib.interval_length, interval_list))

    return map(lambda interval: lib.interval_extend(interval, interval_length), interval_list)


if __name__ == "__main__":
    args = parse_args()
    interval_center_extender(args.bedfile, args.length)
    # TODO: write to output file
