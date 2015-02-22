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


def bed_center_extender(args_list):
    bed_file = args_list[0]
    length = args_list[1]
    outfile = args_list[2]

    bed_line_list = bed_file.readlines()
    # TODO: try rewrite
    bed_param_line_list = map(lib.parse_interval_line, bed_line_list)

    if length is not None:
        interval_length = length
    else:
        interval_length = max(map(lib.interval_length, bed_param_line_list))

    return map(lib.interval_extend, bed_param_line_list, [interval_length] * len(bed_param_line_list))


if __name__ == "__main__":
    args = parse_args()
    args_list1 = [args.bedfile, args.length, args.outfile]
    bed_center_extender(args_list1)
    # TODO: write to output file
