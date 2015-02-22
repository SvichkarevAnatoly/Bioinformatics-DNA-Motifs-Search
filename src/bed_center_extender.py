import argparse
import lib


def create_parser():
    parser = argparse.ArgumentParser(description="Resulting intervals in the input file to the same length")
    parser.add_argument("bedfile", type=argparse.FileType('r'), help="file with bed format intervals")
    parser.add_argument("-l", dest="length", type=int, metavar='length', help="common extended length")
    # TODO: remove outfile if args error
    parser.add_argument("-o", dest="outfile", type=argparse.FileType('w'), metavar='outfile',
                        help="output file with extended bed format intervals")
    return parser


def interval_center_extender(bed_file, length):
    interval_list = bed_file.readlines()
    interval_list = map(lib.parse_interval_line, interval_list)

    if length is not None:
        interval_length = length
    else:
        interval_length = max(map(lib.interval_length, interval_list))

    return map(lambda interval: lib.interval_extend(interval, interval_length), interval_list)


def write_bed_file(interval_param_list, outfile):
    for interval_param in interval_param_list:
        outfile.write("%s\n" % lib.interval_param_list_to_str(interval_param))


def workflow(args):
    bedfile = args.bedfile
    extended_interval_list = interval_center_extender(bedfile, args.length)
    output_file = args.outfile
    if output_file is None:
        output_file = lib.create_output_file_name(bedfile.name)
        output_file = open(output_file, 'w')
    write_bed_file(extended_interval_list, output_file)
    output_file.close()
    bedfile.close()


def main():
    parser = create_parser()
    args = parser.parse_args()
    workflow(args)


if __name__ == "__main__":
    main()