import argparse

parser = argparse.ArgumentParser(description="Resulting intervals in the input file to the same length")
parser.add_argument("bedfile", type=argparse.FileType('r'), help="file with bed format intervals")
parser.add_argument("-l", dest="length", type=int, metavar='length', help="common extended length")
# TODO: remove outfile if args error
parser.add_argument("-o", dest="outfile", type=argparse.FileType('w'), metavar='outfile',
                    help="output file with extended bed format intervals")

args = parser.parse_args()
print '\n'.join(map(str, [args.bedfile, args.length, args.outfile]))


# TODO: read bed file with intervals different length


# TODO: compute central range 500 length


# TODO: write to output file